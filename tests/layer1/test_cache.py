import json
import os
import tempfile
from datetime import datetime, timedelta, timezone

from regvar.layer1.cache import HTTPCache, make_fetch


def test_cache_storage_and_retrieval():
    with tempfile.TemporaryDirectory() as tmpdir:
        cache = HTTPCache(cache_dir=tmpdir, enabled=True, ttl_days=7)
        cache.put("GET", "https://api.example.com/test", {"a": 1}, None, {"result": "ok"}, None)
        entry, hit = cache.get("GET", "https://api.example.com/test", {"a": 1}, None)
        assert hit is True
        assert entry["data"] == {"result": "ok"}
        assert entry["error"] is None
        assert cache.hits == 1


def test_cache_rejects_graphql_errors():
    """Cache poisoning prevention: HTTP 200 with GraphQL 'errors' key must not be persisted."""
    with tempfile.TemporaryDirectory() as tmpdir:
        cache = HTTPCache(cache_dir=tmpdir, enabled=True, ttl_days=7)
        payload = {"errors": [{"message": "Rate limit exceeded"}], "data": None}
        cache.put("POST", "https://api.example.com/graphql", None, {"query": "q"}, payload, None)
        entry, hit = cache.get("POST", "https://api.example.com/graphql", None, {"query": "q"})
        assert hit is False
        assert entry is None
        # Verify 0 files written
        assert len(os.listdir(tmpdir)) == 0


def test_cache_handles_404():
    """Persistent 404 should be cached and returned as a hit on second call."""
    with tempfile.TemporaryDirectory() as tmpdir:
        cache = HTTPCache(cache_dir=tmpdir, enabled=True, ttl_days=7)
        cache.put("GET", "https://api.example.com/missing", None, None, None, "http_404")
        entry, hit = cache.get("GET", "https://api.example.com/missing", None, None)
        assert hit is True
        assert entry["error"] == "http_404"


def test_cache_ttl_expiration():
    with tempfile.TemporaryDirectory() as tmpdir:
        cache = HTTPCache(cache_dir=tmpdir, enabled=True, ttl_days=7)
        cache.put("GET", "https://api.example.com/old", None, None, {"val": 1}, None)
        ckey = cache._cache_key("GET", "https://api.example.com/old", None, None)
        cpath = os.path.join(tmpdir, f"{ckey}.json")
        # Artificially age the entry to 8 days ago
        expired_time = (datetime.now(timezone.utc) - timedelta(days=8)).isoformat()
        with open(cpath, "w", encoding="utf-8") as fh:
            json.dump({"data": {"val": 1}, "error": None, "fetched_at": expired_time}, fh)

        entry, hit = cache.get("GET", "https://api.example.com/old", None, None)
        assert hit is False


def test_cache_disabled_writes_no_files():
    with tempfile.TemporaryDirectory() as tmpdir:
        cache = HTTPCache(cache_dir=tmpdir, enabled=False)
        cache.put("GET", "https://api.example.com/test", None, None, {"data": 1}, None)
        assert len(os.listdir(tmpdir)) == 0
        entry, hit = cache.get("GET", "https://api.example.com/test", None, None)
        assert hit is False
