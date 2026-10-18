"""Implements CONTRIBUTING.md rule #2: provenance is data, not a side effect.

HTTPCache provides transparent on-disk caching for all Layer 1 network calls
(gnomAD, GWAS Catalog, GTEx Portal, Ensembl). This prevents:
- Hitting API rate limits on repeated notebook runs or CI reruns.
- Transient failures poisoning permanent cache entries.

Caching behaviour:
- Successful (200) responses are cached indefinitely (up to cache_ttl_days).
- GraphQL responses whose top-level payload contains an ``errors`` key are
  *never* cached -- they represent a server-side problem that may resolve
  on a retry, not a stable result.
- Persistent HTTP 404 responses are cached so the same missing resource is
  not re-fetched on every call.
- Transient server errors (429, 500, 502, 503, 504) are never cached.
- All entries carry a ``fetched_at`` ISO timestamp; entries older than
  ``cache_ttl_days`` are treated as misses and re-fetched.
"""
from __future__ import annotations

import hashlib
import json
import os
import time
from datetime import datetime, timedelta, timezone
from typing import Any, Optional, Tuple

import requests


class HTTPCache:
    """Transparent on-disk response cache keyed by SHA-256 of the request signature."""

    def __init__(
        self,
        cache_dir: str,
        enabled: bool = True,
        ttl_days: int = 7,
        refresh: bool = False,
    ) -> None:
        self.cache_dir = cache_dir
        self.enabled = enabled
        self.ttl_days = ttl_days
        self.refresh = refresh
        self.hits = 0
        self.misses = 0
        self.oldest_entry_utc: Optional[str] = None
        if self.enabled and self.cache_dir:
            os.makedirs(self.cache_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _cache_key(self, method: str, url: str, params, json_body) -> str:
        payload = {
            "m": method.upper(),
            "u": url,
            "p": sorted(params.items()) if params else None,
            "b": json_body,
        }
        raw = json.dumps(payload, sort_keys=True, default=str).encode("utf-8")
        return hashlib.sha256(raw).hexdigest()

    def _path(self, key: str) -> str:
        return os.path.join(self.cache_dir, f"{key}.json")

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def get(
        self, method: str, url: str, params, json_body
    ) -> Tuple[Optional[dict], bool]:
        """Return (entry_dict, cache_hit).

        ``entry_dict`` contains ``data`` and ``error`` keys matching what was
        stored by :meth:`put`.  Returns ``(None, False)`` on a miss or when
        the cache is disabled.
        """
        if not self.enabled or self.refresh or not self.cache_dir:
            return None, False
        cpath = self._path(self._cache_key(method, url, params, json_body))
        if not os.path.exists(cpath):
            return None, False
        try:
            with open(cpath, "r", encoding="utf-8") as fh:
                entry = json.load(fh)
            fetched_at_str = entry.get("fetched_at")
            if fetched_at_str:
                age = datetime.now(timezone.utc) - datetime.fromisoformat(fetched_at_str)
                if age > timedelta(days=self.ttl_days):
                    return None, False  # expired
                if (
                    self.oldest_entry_utc is None
                    or fetched_at_str < self.oldest_entry_utc
                ):
                    self.oldest_entry_utc = fetched_at_str
            self.hits += 1
            return entry, True
        except Exception:
            return None, False

    def put(
        self,
        method: str,
        url: str,
        params,
        json_body,
        data: Any,
        error: Optional[str],
    ) -> None:
        """Persist a response to disk.

        Rejects:
        - any error other than ``"http_404"`` (transient failures must not be
          cached so the next call gets a fresh attempt).
        - GraphQL payloads that contain a top-level ``errors`` key (cache
          poisoning prevention: the server returned HTTP 200 but the query
          itself failed).
        """
        if not self.enabled or not self.cache_dir:
            return
        if error is not None and error != "http_404":
            return
        if isinstance(data, dict) and data.get("errors"):
            return
        cpath = self._path(self._cache_key(method, url, params, json_body))
        now_str = datetime.now(timezone.utc).isoformat()
        try:
            with open(cpath, "w", encoding="utf-8") as fh:
                json.dump({"data": data, "error": error, "fetched_at": now_str}, fh)
            if self.oldest_entry_utc is None or now_str < self.oldest_entry_utc:
                self.oldest_entry_utc = now_str
        except Exception:
            pass


# ---------------------------------------------------------------------------
# Fetch factory
# ---------------------------------------------------------------------------


def make_fetch(cache: HTTPCache, timeout: int = 30, max_retries: int = 3, backoff: float = 2.0):
    """Return a ``fetch(method, url, params, json_body)`` callable.

    The returned callable returns ``(payload, error, from_cache)``:

    - ``payload`` -- deserialized JSON body, or ``None`` on failure.
    - ``error``   -- ``None`` on success; a short string code on failure
      (e.g. ``"http_503"``, ``"invalid_json"``).
    - ``from_cache`` -- ``True`` iff the response came from the disk cache;
      callers that must rate-limit live API calls (e.g. gnomAD) should skip
      their inter-request sleep when this is ``True``.
    """

    def fetch(method: str, url: str, params=None, json_body=None):
        # 1. Cache lookup
        entry, hit = cache.get(method, url, params, json_body)
        if hit:
            return entry.get("data"), entry.get("error"), True

        cache.misses += 1
        last_err = "unknown"
        for attempt in range(max_retries + 1):
            try:
                resp = requests.request(
                    method, url, params=params, json=json_body, timeout=timeout
                )
            except requests.exceptions.RequestException as exc:
                last_err = type(exc).__name__
            else:
                if resp.status_code == 200:
                    try:
                        body = resp.json()
                    except ValueError:
                        return None, "invalid_json", False
                    cache.put(method, url, params, json_body, body, None)
                    return body, None, False
                if resp.status_code == 404:
                    cache.put(method, url, params, json_body, None, "http_404")
                    return None, "http_404", False
                if resp.status_code in (429, 500, 502, 503, 504):
                    last_err = f"http_{resp.status_code}"
                else:
                    return None, f"http_{resp.status_code}", False
            if attempt < max_retries:
                time.sleep(backoff**attempt)
        return None, last_err, False

    return fetch
