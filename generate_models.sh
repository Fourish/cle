#!/usr/bin/env bash
set -euo pipefail

BASE_DIRECTORY=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

python -B -u \
    "$BASE_DIRECTORY/src/0003-binding-site-selection/generate_models.py" \
    --config "$BASE_DIRECTORY/src/profiles/pipeline_8xql.json" \
    --output-dir "$BASE_DIRECTORY/work/models" "$@"

