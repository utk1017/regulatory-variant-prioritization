#!/bin/bash
set -e
echo "Downloading bigBedToBed from UCSC..."
wget -q http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed -O scripts/bigBedToBed
chmod +x scripts/bigBedToBed
echo "Done."
