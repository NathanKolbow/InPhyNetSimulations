#!/bin/bash
set -euo pipefail

ARCHIVE="data.zip"
SOURCE="data"

zip -r "$ARCHIVE" "$SOURCE" -x "$SOURCE/**/temp-data/*"