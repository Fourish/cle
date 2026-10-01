#!/usr/bin/env bash
set -euo pipefail

fail() {
    printf '%s\n' "$*" >&2
    exit 2
}

[[ $# -eq 2 ]] || fail "Expected exactly two arguments: system ID and stage"
[[ "$1" =~ ^[Ee]([0-9]+)$ ]] || fail "System ID must have the form E07"
printf -v ID 'E%02d' "$((10#${BASH_REMATCH[1]}))"
step="$2"

PRECISION="single"
DEVICE="0"
UTILS_DIRECTORY="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
BASE_DIRECTORY="$(cd -- "${UTILS_DIRECTORY}/../.." && pwd -P)"
JOB_OPENMM_DIRECTORY="${BASE_DIRECTORY}/data/${ID}/charmm-gui-result/openmm"
JOB_RESULT_DIRECTORY="${BASE_DIRECTORY}/result/${ID}"
OUTPUT_BASE="${JOB_RESULT_DIRECTORY}/step${step}"

restart_file=""
case "$step" in
    6.1)
        input_file="${JOB_OPENMM_DIRECTORY}/step6.1_equilibration.inp"
        step_args=(-b "${JOB_OPENMM_DIRECTORY}/sysinfo.dat" --minimize-platform CPU)
        ;;
    6.[2-6])
        input_file="${JOB_OPENMM_DIRECTORY}/step${step}_equilibration.inp"
        previous_step="6.$((${step#6.} - 1))"
        restart_file="${JOB_RESULT_DIRECTORY}/step${previous_step}.rst"
        step_args=(-irst "$restart_file")
        ;;
    7)
        input_file="${UTILS_DIRECTORY}/step7.inp"
        restart_file="${JOB_RESULT_DIRECTORY}/step6.6.rst"
        step_args=(-irst "$restart_file" --checkpoint-interval 2500000)
        ;;
    *)
        fail "Unsupported stage: $step"
        ;;
esac

[[ -d "$JOB_OPENMM_DIRECTORY" ]] || fail "Missing system input directory: $JOB_OPENMM_DIRECTORY"
[[ -f "$input_file" ]] || fail "Missing stage input: $input_file"
[[ -z "$restart_file" || -f "$restart_file" ]] || fail "Missing previous-stage state: $restart_file"
for extension in rst chk dcd pdb out; do
    [[ ! -e "${OUTPUT_BASE}.${extension}" ]] || fail "Output already exists: ${OUTPUT_BASE}.${extension}"
done

mkdir -p -- "$JOB_RESULT_DIRECTORY"
printf '%s\t%s\n' "$step" "$(date '+%Y%m%d%H%M%S')" > "${JOB_RESULT_DIRECTORY}/CURRENT"
cd -- "$JOB_OPENMM_DIRECTORY"
set -o noclobber

status=0
PYTHONPATH="${JOB_OPENMM_DIRECTORY}${PYTHONPATH:+:${PYTHONPATH}}" \
python -B -u "${UTILS_DIRECTORY}/openmm_run_rich.py" \
    --platform CUDA \
    -i "$input_file" \
    -t "${JOB_OPENMM_DIRECTORY}/toppar.str" \
    -p "${JOB_OPENMM_DIRECTORY}/step5_input.psf" \
    -c "${JOB_OPENMM_DIRECTORY}/step5_input.crd" \
    "${step_args[@]}" \
    -orst "${OUTPUT_BASE}.rst" \
    -ochk "${OUTPUT_BASE}.chk" \
    -odcd "${OUTPUT_BASE}.dcd" \
    -opdb "${OUTPUT_BASE}.pdb" \
    --precision "$PRECISION" \
    --device-index "$DEVICE" \
    > "${OUTPUT_BASE}.out" || status=$?

printf '%s\t%s\texit=%s\n' "$step" "$(date '+%Y%m%d%H%M%S')" "$status" >> "${JOB_RESULT_DIRECTORY}/HISTORY"
exit "$status"
