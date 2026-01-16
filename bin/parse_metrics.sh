#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# Step 18 — Compile Statistics
# Tool: parse_metrics.sh
#
# Usage:
#   bin/parse_metrics.sh <sample_id> [metrics_dir] [dedup_dir] [vcf_dir]
#
# Defaults assume your current pipeline layout:
#   metrics_dir = results/metrics
#   dedup_dir   = results/dedup
#   vcf_dir     = results/vcfs
#
# Output:
#   CSV to stdout (header + one row)
# ------------------------------------------------------------

SAMPLE="${1:-}"
METRICS_DIR="${2:-results/metrics}"
DEDUP_DIR="${3:-results/dedup}"
VCF_DIR="${4:-results/vcfs}"

if [[ -z "$SAMPLE" ]]; then
  echo "ERROR: missing sample id" >&2
  echo "Usage: $0 <sample_id> [metrics_dir] [dedup_dir] [vcf_dir]" >&2
  exit 2
fi

ALIGN_METRICS="${METRICS_DIR}/${SAMPLE}.alignment_metrics.txt"
INSERT_METRICS="${METRICS_DIR}/${SAMPLE}.insert_metrics.txt"
DEDUP_METRICS="${DEDUP_DIR}/${SAMPLE}.metrics.txt"
DEPTH_OUT="${METRICS_DIR}/${SAMPLE}.depth_out.txt"

RAW_SNPS="${VCF_DIR}/${SAMPLE}.raw_snps.vcf"
FILTERED_SNPS="${VCF_DIR}/${SAMPLE}.filtered_snps.vcf"
RAW_SNPS_RECAL="${VCF_DIR}/${SAMPLE}.raw_snps_recal.vcf"
FILTERED_SNPS_FINAL="${VCF_DIR}/${SAMPLE}.filtered_snps_final.vcf"

req() {
  local f="$1"
  if [[ ! -s "$f" ]]; then
    echo "ERROR: missing/empty required file: $f" >&2
    exit 3
  fi
}

req "$ALIGN_METRICS"
req "$INSERT_METRICS"
req "$DEDUP_METRICS"
req "$DEPTH_OUT"
req "$RAW_SNPS"
req "$FILTERED_SNPS"
req "$RAW_SNPS_RECAL"
req "$FILTERED_SNPS_FINAL"

# Extract a value from Picard-style metrics tables (tab-separated)
# selecting a given CATEGORY (e.g., PAIR or UNPAIRED) and a field name.
picard_metric() {
  local file="$1"
  local category="$2"
  local field="$3"

  awk -v cat="$category" -v key="$field" '
    BEGIN { FS="\t" }
    /^## METRICS CLASS/ { in_metrics=1; next }
    in_metrics && $0 ~ /^CATEGORY\t/ {
      for (i=1; i<=NF; i++) h[$i]=i
      next
    }
    in_metrics && $1==cat {
      if (key in h) print $(h[key])
      exit
    }
  ' "$file"
}

# Some Picard tables (e.g., InsertSizeMetrics, MarkDuplicates) don’t have CATEGORY.
picard_first_row_field() {
  local file="$1"
  local field="$2"

  awk -v key="$field" '
    BEGIN { FS="\t" }
    /^## METRICS CLASS/ { in_metrics=1; next }
    in_metrics && $1 ~ /^[A-Za-z_]+/ && $0 ~ /\t/ {
      for (i=1; i<=NF; i++) h[$i]=i
      getline
      if (key in h) print $(h[key])
      exit
    }
  ' "$file"
}

# Decide whether we have PAIR metrics; otherwise fall back to UNPAIRED.
HAS_PAIR="$(awk 'BEGIN{found=0} /^PAIR\t/ {found=1} END{print found}' "$ALIGN_METRICS")"
if [[ "$HAS_PAIR" == "1" ]]; then
  CAT="PAIR"
else
  CAT="UNPAIRED"
fi

TOTAL_READS="$(picard_metric "$ALIGN_METRICS" "$CAT" "TOTAL_READS")"
ALIGNED_READS="$(picard_metric "$ALIGN_METRICS" "$CAT" "PF_READS_ALIGNED")"
PCT_ALIGNED="$(picard_metric "$ALIGN_METRICS" "$CAT" "PCT_PF_READS_ALIGNED")"
READ_LENGTH="$(picard_metric "$ALIGN_METRICS" "$CAT" "READ_LENGTH")"

# Aligned bases: prefer PF_ALIGNED_BASES if present, else PF_HQ_ALIGNED_Q20_BASES if present.
ALIGNED_BASES="$(picard_metric "$ALIGN_METRICS" "$CAT" "PF_ALIGNED_BASES")"
if [[ -z "${ALIGNED_BASES:-}" ]]; then
  ALIGNED_BASES="$(picard_metric "$ALIGN_METRICS" "$CAT" "PF_HQ_ALIGNED_Q20_BASES")"
fi

# Convert PAIR counts to read counts (Picard reports PAIR as pairs, not reads).
if [[ "$CAT" == "PAIR" ]]; then
  READS="$(( TOTAL_READS * 2 ))"
  ALIGNED_READS_NUM="$(( ALIGNED_READS * 2 ))"
  PCT_PAIRED="1.0"
else
  READS="$TOTAL_READS"
  ALIGNED_READS_NUM="$ALIGNED_READS"
  PCT_PAIRED="0.0"
fi

# %Duplicate (MarkDuplicates)
PCT_DUPLICATION="$(picard_first_row_field "$DEDUP_METRICS" "PERCENT_DUPLICATION")"

# Mean insert size
MEAN_INSERT_SIZE="$(picard_first_row_field "$INSERT_METRICS" "MEAN_INSERT_SIZE")"

# VCF counts
vcf_total() {
  awk 'BEGIN{c=0} /^#/ {next} {c++} END{print c}' "$1"
}
vcf_filtered() {
  awk 'BEGIN{c=0} /^#/ {next} $7!="PASS" {c++} END{print c}' "$1"
}

SNPS_TOTAL="$(vcf_total "$RAW_SNPS")"
SNPS_FILTERED="$(vcf_filtered "$FILTERED_SNPS")"

SNPS_RECAL_TOTAL="$(vcf_total "$RAW_SNPS_RECAL")"
SNPS_RECAL_FILTERED="$(vcf_filtered "$FILTERED_SNPS_FINAL")"

# Average coverage from depth file (assumes 3 columns: chrom pos depth)
AVG_COV="$(awk 'BEGIN{sum=0; n=0} {sum+=$3; n++} END{ if(n>0) printf "%.6f", sum/n; else print "NA"}' "$DEPTH_OUT")"

# Output CSV
echo "sample,reads,aligned_reads,pct_aligned,aligned_bases,read_length,pct_paired,pct_duplicate,mean_insert_size,snps_total,snps_filtered,snps_recal_total,snps_recal_filtered,avg_coverage"
echo "${SAMPLE},${READS},${ALIGNED_READS_NUM},${PCT_ALIGNED},${ALIGNED_BASES},${READ_LENGTH},${PCT_PAIRED},${PCT_DUPLICATION},${MEAN_INSERT_SIZE},${SNPS_TOTAL},${SNPS_FILTERED},${SNPS_RECAL_TOTAL},${SNPS_RECAL_FILTERED},${AVG_COV}"