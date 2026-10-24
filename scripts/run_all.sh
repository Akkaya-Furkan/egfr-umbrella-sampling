#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
PYTHON="${PYTHON:-python}"
"$PYTHON" scripts/download_data.py
"$PYTHON" scripts/prepare_egfr.py
systemd-inhibit --what=sleep:idle --why='P23 umbrella model production' \
  bash -c 'systemd-inhibit --list > results/inhibitors.txt; exec "$1" scripts/run_model.py' _ "$PYTHON" \
  > results/model_run.log 2>&1
"$PYTHON" scripts/analyze_model.py > results/analysis_run.log 2>&1
"$PYTHON" -m pytest -q --junitxml=results/test-results.xml
"$PYTHON" scripts/verify_results.py
