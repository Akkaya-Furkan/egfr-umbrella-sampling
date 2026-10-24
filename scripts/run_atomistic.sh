#!/usr/bin/env bash
# Start with: nohup systemd-inhibit --what=sleep:idle --why='P23 atomistic pilot' \
#   bash scripts/run_atomistic.sh > atomistic-run.log 2>&1 < /dev/null &
set -euo pipefail
cd "$(dirname "$0")/.."
export OPENMM_PYTHON="${OPENMM_PYTHON:-$HOME/.conda/envs/openmm_env/bin/python}"
export AMBERHOME="${AMBERHOME:-$HOME/.conda/envs/md_env}"
export OMP_NUM_THREADS=4 OPENMM_CPU_THREADS=4
analysis_python="${ANALYSIS_PYTHON:-.venv/bin/python}"
# One deadline includes downloading, preparation and MD; reserve analysis time.
deadline=$(( $(date +%s) + 9900 ))
if [[ -e results/atomistic/simulation.json ]]; then
  echo 'Existing atomistic run found. Archive it before starting a new run.' >&2
  exit 1
fi
mkdir -p data/atomistic results/atomistic
"$analysis_python" scripts/download_data.py
curl -fL --max-time 60 'https://models.rcsb.org/v1/1m17/ligand?auth_asym_id=A&auth_seq_id=999&encoding=sdf' -o data/atomistic/AQ4.sdf
timeout --signal=TERM --kill-after=30s "$((deadline-$(date +%s)))s" bash scripts/prepare_atomistic.sh > results/atomistic/preparation.log 2>&1
"$AMBERHOME/bin/python" scripts/audit_atomistic.py
md_status=0
timeout --signal=TERM --kill-after=30s "$((deadline-$(date +%s)))s" \
  "$OPENMM_PYTHON" -u scripts/run_atomistic.py --deadline-unix "$deadline" \
  --production-ps "${PRODUCTION_PS:-400}" > results/atomistic/md.log 2>&1 || md_status=$?
"$analysis_python" scripts/analyze_atomistic.py
"$OPENMM_PYTHON" scripts/verify_atomistic.py
exit "$md_status"
