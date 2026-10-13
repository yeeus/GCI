#!/usr/bin/env bash
# Benchmark GCI partitioned-hybrid at several worker counts on the same MF2 data.
# Run from GCI repository root: THREADS="4 8 16 32" bash scripts/benchmark_hybrid_threads.sh
set -euo pipefail

[[ -f GCI.py ]] || { echo "Run from the directory containing GCI.py" >&2; exit 2; }
: "${REFERENCE:=../../MF2_mat.v1.0.fa}"
: "${HIFI_BAM:=../../mat.winnowmap.hifi.bam}"
: "${HIFI_PAF:=../../mat.minimap2.hifi.paf}"
: "${ONT_BAM:=../../mat.winnowmap.ul.bam}"
: "${ONT_PAF:=../../mat.minimap2.ul.paf}"
: "${SCRATCH:=./scratch}"
: "${PARTITIONS:=256}"
: "${INPUT_SHARDS:=16}"
: "${THREADS:=4 8 16 32}"
: "${BASELINE_DIR:=test5}"
: "${RESULTS_PREFIX:=benchmark_hybrid_t}"

mkdir -p "$SCRATCH"
for item in "$REFERENCE" "$HIFI_BAM" "$HIFI_PAF" "$ONT_BAM" "$ONT_PAF"; do
    [[ -r "$item" ]] || { echo "Missing input: $item" >&2; exit 2; }
done

echo "Input shards: $INPUT_SHARDS; read-name partitions: $PARTITIONS"
echo "Scratch filesystem:"
df -Th "$SCRATCH"
echo "Threads to benchmark: $THREADS"

for t in $THREADS; do
    [[ "$t" =~ ^[1-9][0-9]*$ ]] || { echo "Invalid -t: $t" >&2; exit 2; }
    out="${RESULTS_PREFIX}${t}"
    log="${out}.log"
    echo "===== -t $t : $(date -Is) ====="
    /usr/bin/time -v env GCI_PROFILE=1 python -u GCI.py \
        -r "$REFERENCE" \
        --hifi "$HIFI_BAM" "$HIFI_PAF" \
        --nano "$ONT_BAM" "$ONT_PAF" \
        --filter-backend partitioned-hybrid \
        --scratch-dir "$SCRATCH" \
        --partitions "$PARTITIONS" \
        --input-shards "$INPUT_SHARDS" \
        -t "$t" \
        -d "$out" -o mat -f \
        2>&1 | tee "$log"

    if [[ -f "$BASELINE_DIR/mat.gci" ]]; then
        echo "Checking $out against $BASELINE_DIR ..."
        diff -q "$BASELINE_DIR/mat.gci" "$out/mat.gci"
        for type in hifi nano two_type; do
            diff -q "$BASELINE_DIR/mat_${type}.0.depth.bed" \
                    "$out/mat_${type}.0.depth.bed"
            gzip -t "$BASELINE_DIR/mat_${type}.depth.gz"
            gzip -t "$out/mat_${type}.depth.gz"
            cmp <(gzip -cd "$BASELINE_DIR/mat_${type}.depth.gz") \
                <(gzip -cd "$out/mat_${type}.depth.gz")
        done
        echo "Output comparison passed for -t $t."
    else
        echo "WARNING: $BASELINE_DIR/mat.gci absent; output equality NOT checked." >&2
    fi
    echo "Summary of -t $t:"
    grep -E 'Elapsed \(wall clock\) time|User time \(seconds\)|System time \(seconds\)|Percent of CPU this job got|Maximum resident set size \(kbytes\)' "$log" | tail -5
done

echo "Benchmarks completed. Retain logs and compare elapsed time, CPU seconds, peak RSS, and output equality."
