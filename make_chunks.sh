#!/usr/bin/env bash
# make_sbatch_for_chunks.sh
# Usage: ./make_sbatch_for_chunks.sh antismash.test.sbatch [out_dir] [glob]
# Defaults: out_dir=sbatch_chunks, glob="chunk_*"

set -euo pipefail

TEMPLATE="${1:?Usage: $0 TEMPLATE_SBATCH [out_dir] [glob]}"
OUTDIR="${2:-sbatch_chunks}"
GLOB="${3:-chunk_*}"

[[ -f "$TEMPLATE" ]] || { echo "Template not found: $TEMPLATE" >&2; exit 1; }
mkdir -p "$OUTDIR"

base="$(basename "$TEMPLATE")"
base_noext="${base%.sbatch}"

shopt -s nullglob
found=0
for chunkfile in $GLOB; do
  [[ -f "$chunkfile" ]] || continue
  found=1
  chunk="$(basename "$chunkfile")"
  out="$OUTDIR/${base_noext}_${chunk}.sbatch"

  awk -v chunk="$chunk" -v base="$base_noext" '
    BEGIN {
      in_header=1; injected=0
      jobname_done=0; redirect_done=0; runname_done=0
    }

    # Replace obvious placeholders early
    { gsub(/__CHUNK__|{{CHUNK}}|%%CHUNK%%/, chunk) }

    # Make antismash use scheduler threads, not a hardcoded number
    { gsub(/-c[[:space:]]*[0-9]+/, "-c \"$SLURM_CPUS_PER_TASK\"") }

    # Preserve shebang as-is
    NR==1 && /^#!/ { print; next }

    # While we are still in the SBATCH header, adjust job-name if present
    in_header==1 && /^#SBATCH/ {
      if ($0 ~ /--job-name=/) {
        line=$0
        sub(/.*--job-name[[:space:]]*=[[:space:]]*/, "", line)
        print "#SBATCH --job-name=" line "_" chunk
        jobname_done=1
      } else {
        print
      }
      next
    }

    # First non-#SBATCH line -> inject CHUNK/RUNNAME right before it
    in_header==1 {
      if (!injected) {
        print "CHUNK=\"" chunk "\""
        print "RUNNAME=\"" chunk "\""
        injected=1
      }
      in_header=0
    }

    # Normalize any later CHUNK/RUNNAME assignments
    /^[[:space:]]*CHUNK[[:space:]]*=/ { print "CHUNK=\"" chunk "\""; next }
    /^[[:space:]]*RUNNAME[[:space:]]*=/ { print "RUNNAME=\"" chunk "\""; runname_done=1; next }

    # while ... < file ; do  ->  while ... < "$CHUNK" ; do
    /^[[:space:]]*while[^#]*<[^>]*;[[:space:]]*do[[:space:]]*$/ {
      sub(/<.*/, "< \"$CHUNK\" ; do")
      print; redirect_done=1; next
    }

    # done < file  ->  done < "$CHUNK"
    /^[[:space:]]*done[[:space:]]*<.*$/ {
      sub(/<.*/, "< \"$CHUNK\"")
      print; redirect_done=1; next
    }

    # bare "done" -> add the redirection
    /^[[:space:]]*done[[:space:]]*$/ {
      print $0 " < \"$CHUNK\""
      redirect_done=1; next
    }

    { print }

    END {
      if (!jobname_done) {
        print "# NOTE: No #SBATCH --job-name= found; consider adding one (e.g., antismash_" chunk ")."
      }
      if (!redirect_done) {
        print "# NOTE: No while/done input redirection found to rewrite."
        print "# Add:   done < \"$CHUNK\"   at the end of your while loop."
      }
    }
  ' "$TEMPLATE" > "$out"

  chmod +x "$out"
  echo "Wrote: $out"
done
shopt -u nullglob

if [[ $found -eq 0 ]]; then
  echo "No chunk files matched: $GLOB" >&2
  exit 2
fi

echo "Done. Generated files in: $OUTDIR"



