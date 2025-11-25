#!/usr/bin/env bash
set -euo pipefail

# Usage: ./download_fastq.sh <accessions_file> <outdir_fastq> <out_samples_tsv>
# Example: ./download_fastq.sh accessions.txt fastq/ samples.tsv

ACCESSIONS_FILE=${1:-accessions.txt}
OUTDIR=${2:-fastq}
OUT_SAMPLES=${3:-samples.tsv}

mkdir -p "$OUTDIR"
rm -f "$OUT_SAMPLES"
echo -e "sample\tr1\tr2" > "$OUT_SAMPLES"

while read -r ACC; do
  # skip empty or comment lines
  [[ -z "$ACC" || "$ACC" =~ ^# ]] && continue

  echo "Processing $ACC ..."

  # Use SRA Toolkit: prefetch + fasterq-dump (split paired files)
  # This stores .sra in ~/.ncbi/public/sra by default
  if command -v prefetch >/dev/null 2>&1 && command -v fasterq-dump >/dev/null 2>&1; then
    echo "Using SRA Toolkit for $ACC"
    prefetch --max-size 100G "$ACC" || { echo "prefetch failed for $ACC"; exit 1; }
    # convert to fastq, split into R1/R2 if paired
    fasterq-dump --split-files --threads 4 --outdir "$OUTDIR" "$ACC"
    # fasterq-dump writes files like ACC_1.fastq and ACC_2.fastq (or ACC.fastq for single-end)
    if [[ -f "${OUTDIR}/${ACC}_1.fastq" ]]; then
      gzip -f "${OUTDIR}/${ACC}_1.fastq"
      gzip -f "${OUTDIR}/${ACC}_2.fastq" || true
      R1="${OUTDIR}/${ACC}_1.fastq.gz"
      R2="${OUTDIR}/${ACC}_2.fastq.gz"
    elif [[ -f "${OUTDIR}/${ACC}.fastq" ]]; then
      gzip -f "${OUTDIR}/${ACC}.fastq"
      R1="${OUTDIR}/${ACC}.fastq.gz"
      R2=""
    else
      echo "Expected FASTQ output not found for $ACC"
      exit 1
    fi

  else
    # fallback: try ENA FTP (works for public ENA runs)
    # This will try to download paired files via FTP from ENA
    echo "SRA Toolkit not available; trying ENA for $ACC"
    # ENA run accession mapping: use simple pattern (may not work for all)
    # prefer to try ENA run direct URLs:
    ENA_BASE=https://www.ebi.ac.uk/ena/portal/api/fastq
    # request URLs for the run (comma-separated if multiple)
    curl -s "${ENA_BASE}?accession=${ACC}&download=true" -o "${OUTDIR}/${ACC}.fastq.gz" || true
    if [[ -f "${OUTDIR}/${ACC}.fastq.gz" ]]; then
      R1="${OUTDIR}/${ACC}.fastq.gz"
      R2=""
    else
      echo "ENA download failed for $ACC. Install sra-tools or check accession."
      exit 1
    fi
  fi

  # create a sample name — prefer SRR/ERR stripped to ID; you can change mapping later
  SAMPLE="${ACC}"
  # If R2 empty, put NA
  if [[ -n "${R2:-}" && -f "$R2" ]]; then
    echo -e "${SAMPLE}\t${R1}\t${R2}" >> "$OUT_SAMPLES"
  else
    echo -e "${SAMPLE}\t${R1}\tNA" >> "$OUT_SAMPLES"
  fi

done < "$ACCESSIONS_FILE"

echo "Done. Sample table: $OUT_SAMPLES"
