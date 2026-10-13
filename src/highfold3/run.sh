#!/usr/bin/env bash
set -euo pipefail

fail() {
    printf '%s\n' "$*" >&2
    exit 2
}

usage() {
    printf '%s\n' \
        "Usage: bash src/highfold3/run.sh JOB data|inference|all" \
        "Input: data/highfold3/inputs/JOB.json" \
        "Output: result/highfold3/JOB/STAGE/; log: STAGE.out" \
        "Run in the prepared highfold3 environment."
}

if [[ $# -eq 1 && ( "$1" == "--help" || "$1" == "-h" ) ]]; then
    usage
    exit 0
fi
[[ $# -eq 2 ]] || fail "Expected two arguments: job name and stage (data, inference, all)"
[[ "$1" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || fail "Job name may contain letters, digits, underscores, dots and hyphens"
ID="$1"
step="$2"

UTILS_DIRECTORY="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
BASE_DIRECTORY="$(cd -- "${UTILS_DIRECTORY}/../.." && pwd -P)"
JOB_INPUT_FILE="${BASE_DIRECTORY}/data/highfold3/inputs/${ID}.json"
JOB_RESULT_DIRECTORY="${BASE_DIRECTORY}/result/highfold3/${ID}"
OUTPUT_DIRECTORY="${JOB_RESULT_DIRECTORY}/${step}"
OUTPUT_LOG="${JOB_RESULT_DIRECTORY}/${step}.out"

MODEL_DIRECTORY="${HIGHFOLD3_MODEL_DIR:-/root/fourish/alphafold3/model}"
DATABASE_DIRECTORY="${HIGHFOLD3_DB_DIR:-/root/fourish/alphafold3/database}"
DEVICE="${HIGHFOLD3_DEVICE:-0}"
CPU_THREADS="${HIGHFOLD3_CPUS:-4}"
BUCKETS="${HIGHFOLD3_BUCKETS:-256,1024}"
ATTENTION="${HIGHFOLD3_ATTENTION:-triton}"

input_file="$JOB_INPUT_FILE"
case "$step" in
    data)
        stage_args=(--run_data_pipeline=true --run_inference=false)
        ;;
    inference)
        stage_args=(--run_data_pipeline=false --run_inference=true)
        if [[ -d "${JOB_RESULT_DIRECTORY}/data" ]]; then
            shopt -s nullglob
            prepared_inputs=("${JOB_RESULT_DIRECTORY}/data/"*/*_data.json)
            [[ ${#prepared_inputs[@]} -eq 1 ]] || fail "Expected one prepared JSON from the data stage: ${JOB_RESULT_DIRECTORY}/data"
            input_file="${prepared_inputs[0]}"
        fi
        ;;
    all)
        stage_args=(--run_data_pipeline=true --run_inference=true)
        ;;
    *)
        fail "Unsupported stage: $step (expected data, inference or all)"
        ;;
esac

[[ -f "$input_file" ]] || fail "Missing prediction input: $input_file"
[[ "$step" == data || -d "$MODEL_DIRECTORY" ]] || fail "Missing model directory; set HIGHFOLD3_MODEL_DIR: $MODEL_DIRECTORY"
[[ "$step" == inference || -d "$DATABASE_DIRECTORY" ]] || fail "Missing database directory; set HIGHFOLD3_DB_DIR: $DATABASE_DIRECTORY"
[[ ! -e "$OUTPUT_DIRECTORY" ]] || fail "Output already exists: $OUTPUT_DIRECTORY"
[[ ! -e "$OUTPUT_LOG" ]] || fail "Output already exists: $OUTPUT_LOG"

mkdir -p -- "$JOB_RESULT_DIRECTORY"
mkdir -- "$OUTPUT_DIRECTORY"
printf '%s\t%s\n' "$step" "$(date '+%Y%m%d%H%M%S')" > "${JOB_RESULT_DIRECTORY}/CURRENT"
printf 'Input: %s\nOutput: %s\nLog: %s\n' "$input_file" "$OUTPUT_DIRECTORY" "$OUTPUT_LOG"
set -o noclobber

status=0
python -B -u "${UTILS_DIRECTORY}/run_alphafold.py" \
    --json_path "$input_file" \
    --output_dir "$OUTPUT_DIRECTORY" \
    --model_dir "$MODEL_DIRECTORY" \
    --db_dir "$DATABASE_DIRECTORY" \
    --head_to_tail=true \
    --gpu_device "$DEVICE" \
    --buckets "$BUCKETS" \
    --flash_attention_implementation "$ATTENTION" \
    --num_recycles 10 \
    --num_diffusion_samples 5 \
    --jackhmmer_n_cpu "$CPU_THREADS" \
    --nhmmer_n_cpu "$CPU_THREADS" \
    --max_template_date 2026-10-12 \
    "${stage_args[@]}" \
    > "$OUTPUT_LOG" 2>&1 || status=$?

printf '%s\t%s\texit=%s\n' "$step" "$(date '+%Y%m%d%H%M%S')" "$status" >> "${JOB_RESULT_DIRECTORY}/HISTORY"
if [[ "$status" -ne 0 ]]; then
    printf 'HighFold3 failed (exit %s); see %s\n' "$status" "$OUTPUT_LOG" >&2
fi
exit "$status"
