#!/usr/bin/env bash
set -euo pipefail

# 18_generate_rmats_event_set.sh
#
# Purpose:
#   Generate rMATS event-definition files for downstream MARVEL:
#     fromGTF.SE.txt
#     fromGTF.MXE.txt
#     fromGTF.A5SS.txt
#     fromGTF.A3SS.txt
#     fromGTF.RI.txt
#
# We are NOT using rMATS here for differential splicing statistics.
# We are only using it to define the event universe from:
#   - the GTF
#   - the STAR-aligned BAMs
#
# Dataset-specific assumptions:
#   - Smart-seq2
#   - single-end
#   - read length = 100
#   - likely unstranded
#
# Run from project root:
#   bash scripts/18_generate_rmats_event_set.sh

PROJECT_DIR="$HOME/Documents/mouse_oocyte_project"
cd "$PROJECT_DIR"

# -----------------------------
# paths / settings
# -----------------------------
METADATA="data/derived/metadata/sample_metadata_with_predicted_configuration.tsv"
GTF="reference/Mus_musculus.GRCm38.102.gtf"

OUTDIR="results/rmats_event_set"
TMPDIR="results/rmats_event_set_tmp"

B1="$OUTDIR/b1_young.txt"
B2="$OUTDIR/b2_old.txt"

READ_LENGTH=100
THREADS=8
LIBTYPE="fr-unstranded"

# Change this if needed.
# Examples:
#   RMATS_CMD="python /path/to/rmats.py"
#   RMATS_CMD="rmats.py"
#RMATS_CMD="rmats.py"
#RMATS_CMD="$HOME/rmats-turbo/run_rmats"

RMATS_PY="$HOME/rmats-turbo/rmats.py"
RMATS_PYTHON="$HOME/miniconda3/bin/python"

mkdir -p "$OUTDIR" "$TMPDIR"

# -----------------------------
# checks
# -----------------------------
if [[ ! -f "$METADATA" ]]; then
  echo "ERROR: metadata file not found: $METADATA" >&2
  exit 1
fi

if [[ ! -f "$GTF" ]]; then
  echo "ERROR: GTF not found: $GTF" >&2
  exit 1
fi

if [[ ! -f "$RMATS_PY" ]]; then
  echo "ERROR: rMATS script not found: $RMATS_PY" >&2
  exit 1
fi

if [[ ! -x "$RMATS_PYTHON" ]]; then
  echo "ERROR: Python executable not found: $RMATS_PYTHON" >&2
  exit 1
fi

# -----------------------------
# build BAM lists from metadata
# -----------------------------
python3 - <<'PY'
import pandas as pd
from pathlib import Path

metadata = Path("data/derived/metadata/sample_metadata_with_predicted_configuration.tsv")
outdir = Path("results/rmats_event_set")
b1 = outdir / "b1_young.txt"
b2 = outdir / "b2_old.txt"

df = pd.read_csv(metadata, sep="\t")

required = ["run_clean", "age_group_clean", "bam_path"]
missing = [c for c in required if c not in df.columns]
if missing:
    raise SystemExit(f"Missing required metadata columns: {missing}")

df = df[required].drop_duplicates()

young = df.loc[df["age_group_clean"] == "young", "bam_path"].astype(str).tolist()
old   = df.loc[df["age_group_clean"] == "old",   "bam_path"].astype(str).tolist()

if len(young) == 0:
    raise SystemExit("No young BAMs found in metadata.")
if len(old) == 0:
    raise SystemExit("No old BAMs found in metadata.")

for bam in young + old:
    if not Path(bam).exists():
        raise SystemExit(f"BAM file not found: {bam}")

# rMATS wants each file to contain one comma-separated list of BAMs
b1.write_text(",".join(young) + "\n")
b2.write_text(",".join(old) + "\n")

print(f"Wrote {len(young)} young BAMs to {b1}")
print(f"Wrote {len(old)} old BAMs to {b2}")
PY

echo
echo "Young BAM list file:"
cat "$B1"
echo
echo "Old BAM list file:"
cat "$B2"
echo

# -----------------------------
# run rMATS
# -----------------------------
"$RMATS_PYTHON" "$RMATS_PY" \
  --b1 "$B1" \
  --b2 "$B2" \
  --gtf "$GTF" \
  -t single \
  --readLength "$READ_LENGTH" \
  --libType "$LIBTYPE" \
  --nthread "$THREADS" \
  --od "$OUTDIR" \
  --tmp "$TMPDIR" \
  --statoff

# -----------------------------
# check outputs
# -----------------------------
echo
echo "Checking for expected rMATS event files..."

expected=(
  "fromGTF.SE.txt"
  "fromGTF.MXE.txt"
  "fromGTF.A5SS.txt"
  "fromGTF.A3SS.txt"
  "fromGTF.RI.txt"
)

missing_any=0
for f in "${expected[@]}"; do
  if [[ -f "$OUTDIR/$f" ]]; then
    echo "FOUND: $OUTDIR/$f"
  else
    echo "MISSING: $OUTDIR/$f"
    missing_any=1
  fi
done

echo
echo "Top-level output files:"
find "$OUTDIR" -maxdepth 1 -type f | sort

if [[ "$missing_any" -ne 0 ]]; then
  echo
  echo "WARNING: One or more expected fromGTF files are missing."
  echo "Inspect the rMATS output above and the files in:"
  echo "  $OUTDIR"
  exit 1
fi

echo
echo "Done."
echo "rMATS event-definition files are in:"
echo "  $OUTDIR"
