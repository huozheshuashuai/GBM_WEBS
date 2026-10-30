#!/bin/bash
# pget.sh <url> <output> [segments]
# Parallel segmented download using curl Range requests
set -e
URL="$1"
OUT="$2"
N="${3:-16}"

SIZE=$(curl -sI "$URL" | grep -i '^Content-Length:' | tail -1 | tr -dc '0-9')
if [ -z "$SIZE" ]; then echo "Cannot determine size for $URL"; exit 1; fi
echo "Downloading $OUT ($SIZE bytes, $N segments)"

CHUNK=$(( (SIZE + N - 1) / N ))
DIR="${OUT}.parts"
rm -rf "$DIR"; mkdir -p "$DIR"

pids=()
for ((i=0; i<N; i++)); do
  START=$((i * CHUNK))
  END=$(( (i+1) * CHUNK - 1 ))
  if [ $END -ge $SIZE ]; then END=$((SIZE - 1)); fi
  if [ $START -gt $END ]; then break; fi
  curl -s --retry 5 --retry-delay 3 -C - -H "Range: bytes=${START}-${END}" -o "$DIR/part_$(printf %03d $i)" "$URL" &
  pids+=($!)
done

fail=0
for p in "${pids[@]}"; do wait "$p" || fail=1; done

# verify sizes
TOTAL=0
for f in "$DIR"/part_*; do
  S=$(stat -c %s "$f"); TOTAL=$((TOTAL + S))
done
if [ "$TOTAL" -ne "$SIZE" ] || [ "$fail" -ne 0 ]; then
  echo "INCOMPLETE: got $TOTAL / $SIZE bytes (fail=$fail)"
  exit 2
fi

cat "$DIR"/part_* > "$OUT"
rm -rf "$DIR"
echo "DONE: $OUT ($SIZE bytes)"
