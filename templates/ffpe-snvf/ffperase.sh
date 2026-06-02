#!/usr/bin/env bash
set -euo pipefail

echo "python path:    $(which python3)"
echo "sklearn version: $(python3 -W ignore -c 'import sklearn; print(sklearn.__version__)')"  # must be 0.24.1

# ---- Workflow path ----
MAIN_NF="/home/moyukh/storage_dshihlab/repos/nf-ffperase/main.nf"
MODEL="/home/moyukh/storage_dshihlab/repos/nf-ffperase/trained_models/model.snvs.pkl"

# ---- Defaults ----
VCF=""
BAM=""
REFERENCE=""
BED=""
COVERAGE=""
MEDIAN_INSERT=""
STEP="full"
MUTATION_TYPE="snvs"
MODEL_NAME="ARTIFACT"
SAMPLE_NAME=""        # auto-derived from VCF if empty
OUTDIR=""             # auto -> results/<sample>
RUN_DIR=""            # auto -> runs/<sample>
WORK_DIR=""           # auto -> nf-work/<sample>
PROFILE=""            # optional: cloud | hpc_slurm | test
RESUME=0
EXTRA_ARGS=()         # anything after `--` is forwarded to nextflow

usage() {
    cat <<EOF
Usage: $(basename "$0") [options] [-- <extra nextflow args>]

Required:
  --vcf PATH              Input VCF
  --bam PATH              Input BAM
  --reference PATH        Reference FASTA
  --bed PATH              BED file
  --coverage INT          Median coverage
  --median-insert INT     Median insert size

Optional (have defaults):
  --model PATH            Trained model pickle            [${MODEL}]
  --model-name NAME       Model name (column label)       [${MODEL_NAME}]
  --step STEP             preprocess|classify|train|full  [${STEP}]
  --mutation-type TYPE    snvs|indels                     [${MUTATION_TYPE}]
  --sample-name NAME      Override derived sample name    [auto from --vcf]
  --outdir DIR            Output directory                [results/<sample>]
  --run-dir DIR           Per-run logs/reports dir        [runs/<sample>]
  --work-dir DIR          Nextflow work dir               [nf-work/<sample>]
  --profile NAME          Nextflow -profile               [none]
  --resume                Pass -resume to nextflow
  -h, --help              Show this help and exit

Anything after a literal '--' is forwarded verbatim to nextflow, e.g.:
  $(basename "$0") --step classify -- -with-tower
EOF
}

# Helper: print error + usage + exit
die() {
    echo "ERROR: $*" >&2
    echo >&2
    usage >&2
    exit 2
}

# Helper: assert a variable is non-empty
require_param() {
    local flag="$1" value="$2"
    [[ -n "${value}" ]] || die "Missing required argument: ${flag}"
}

# Helper: format seconds as Hh Mm Ss
format_duration() {
    local total=$1
    local h=$(( total / 3600 ))
    local m=$(( (total % 3600) / 60 ))
    local s=$(( total % 60 ))
    printf '%dh%02dm%02ds' "${h}" "${m}" "${s}"
}

# ---- Parse args ----
while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf)            VCF="$2";            shift 2 ;;
        --bam)            BAM="$2";            shift 2 ;;
        --reference)      REFERENCE="$2";      shift 2 ;;
        --bed)            BED="$2";            shift 2 ;;
        --model)          MODEL="$2";          shift 2 ;;
        --model-name)     MODEL_NAME="$2";     shift 2 ;;
        --coverage)       COVERAGE="$2";       shift 2 ;;
        --median-insert)  MEDIAN_INSERT="$2";  shift 2 ;;
        --step)           STEP="$2";           shift 2 ;;
        --mutation-type)  MUTATION_TYPE="$2";  shift 2 ;;
        --sample-name)    SAMPLE_NAME="$2";    shift 2 ;;
        --outdir)         OUTDIR="$2";         shift 2 ;;
        --run-dir)        RUN_DIR="$2";        shift 2 ;;
        --work-dir)       WORK_DIR="$2";       shift 2 ;;
        --profile)        PROFILE="$2";        shift 2 ;;
        --resume)         RESUME=1;            shift   ;;
        -h|--help)        usage; exit 0 ;;
        --)               shift; EXTRA_ARGS=("$@"); break ;;
        *)
            die "Unknown argument: $1"
            ;;
    esac
done

# ---- Validate required arguments ----
require_param --vcf            "${VCF}"
require_param --bam            "${BAM}"
require_param --reference      "${REFERENCE}"
require_param --bed            "${BED}"
require_param --coverage       "${COVERAGE}"
require_param --median-insert  "${MEDIAN_INSERT}"

# ---- Validate that input files actually exist ----
for f in "${VCF}" "${BAM}" "${REFERENCE}" "${BED}" "${MODEL}"; do
    [[ -e "${f}" ]] || die "File not found: ${f}"
done

# ---- Derive sample name from VCF if not provided ----
# e.g. "testdata/FFPE-Colon-Tumoral_B00GXHK.vcf.gz" -> "FFPE-Colon-Tumoral_B00GXHK"
if [[ -z "${SAMPLE_NAME}" ]]; then
    SAMPLE_NAME="$(basename "${VCF}")"
    SAMPLE_NAME="${SAMPLE_NAME%%.*}"
fi

# ---- Defaulted output/run paths (derived after sample name is known) ----
RUN_ID="${SAMPLE_NAME}"
: "${OUTDIR:=results/${RUN_ID}}"
: "${RUN_DIR:=nf-runs/${RUN_ID}}"
: "${WORK_DIR:=nf-work/${RUN_ID}}"
mkdir -p "${RUN_DIR}"

# ---- Isolate Nextflow's .nextflow metadata (history/cache) per sample ----
export NXF_CACHE_DIR="$(readlink -f "${RUN_DIR}")/.nextflow"
mkdir -p "${NXF_CACHE_DIR}"

# ---- Runtime tracking output ----
RUNTIME_TSV="${OUTDIR}/${SAMPLE_NAME}.runtime.tsv"

# ---- Build optional nextflow args ----
NF_OPTS=()
[[ -n "${PROFILE}" ]] && NF_OPTS+=(-profile "${PROFILE}")
[[ "${RESUME}" -eq 1 ]] && NF_OPTS+=(-resume)

echo "Launching FFPERASE for sample: ${SAMPLE_NAME}"
echo "  Step:    ${STEP}"
echo "  VCF:     ${VCF}"
echo "  BAM:     ${BAM}"
echo "  Outdir:  ${OUTDIR}"
echo "  Run ID:  ${RUN_ID}"
[[ -n "${PROFILE}" ]] && echo "  Profile: ${PROFILE}"
echo "  Runtime log: ${RUNTIME_TSV}"

# ---- Time the nextflow execution ----
# Capture exit status without tripping `set -e`, so we always write the TSV.
START_EPOCH=$(date +%s)
START_ISO=$(date -u -d "@${START_EPOCH}" +"%Y-%m-%dT%H:%M:%SZ")

NF_EXIT=0
nextflow \
    -log "${RUN_DIR}/nextflow.log" \
    run "${MAIN_NF}" \
    -work-dir      "${WORK_DIR}" \
    -with-report   "${RUN_DIR}/report.html" \
    -with-timeline "${RUN_DIR}/timeline.html" \
    -with-dag      "${RUN_DIR}/dag.html" \
    -with-trace    "${RUN_DIR}/trace.txt" \
    "${NF_OPTS[@]}" \
    --step          "${STEP}" \
    --vcf           "${VCF}" \
    --bam           "${BAM}" \
    --reference     "${REFERENCE}" \
    --bed           "${BED}" \
    --outdir        "${OUTDIR}" \
    --coverage      "${COVERAGE}" \
    --medianInsert  "${MEDIAN_INSERT}" \
    --model         "${MODEL}" \
    --modelName     "${MODEL_NAME}" \
    --mutationType  "${MUTATION_TYPE}" \
    --sampleName    "${SAMPLE_NAME}" \
    "${EXTRA_ARGS[@]}" || NF_EXIT=$?

END_EPOCH=$(date +%s)
END_ISO=$(date -u -d "@${END_EPOCH}" +"%Y-%m-%dT%H:%M:%SZ")
DURATION_S=$(( END_EPOCH - START_EPOCH ))
DURATION_H=$(format_duration "${DURATION_S}")
STATUS=$([[ "${NF_EXIT}" -eq 0 ]] && echo "SUCCESS" || echo "FAILED")

# ---- Write runtime TSV ----
{
    printf 'sample_name\tstep\tmutation_type\tprofile\tstart_time_utc\tend_time_utc\tduration_seconds\tduration_human\texit_code\tstatus\n'
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${SAMPLE_NAME}" \
        "${STEP}" \
        "${MUTATION_TYPE}" \
        "${PROFILE:-none}" \
        "${START_ISO}" \
        "${END_ISO}" \
        "${DURATION_S}" \
        "${DURATION_H}" \
        "${NF_EXIT}" \
        "${STATUS}"
} > "${RUNTIME_TSV}"

echo
echo "Runtime: ${DURATION_H} (${DURATION_S}s) [${STATUS}]"
echo "Wrote: ${RUNTIME_TSV}"

exit "${NF_EXIT}"

