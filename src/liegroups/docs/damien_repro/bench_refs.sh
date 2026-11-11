#!/bin/bash
# usage: bench_refs.sh [ref ...]   (défaut : master de Damien, sa branche, la branche d'audit)
# Compile probe_bench.cpp en -O2 contre chaque ref. Lancer plusieurs fois : le temps varie de quelques %.
source "$(dirname "$0")/common.sh"
REFS=("$@"); [ ${#REFS[@]} -eq 0 ] && REFS=(fa71b1b defrost/damien-refactoring painless/liegroups-audit)
for ref in "${REFS[@]}"; do
  wt=$(checkout_ref "$ref") || continue
  echo "=== $ref ($(git -C "$REPO" rev-parse --short "$ref"))"
  if "$CXX" -std=gnu++20 -O2 -w -I"$EIGEN_INC" -I"$wt/src" "$HERE/probe_bench.cpp" -o "$WORK/b.bin" 2>"$WORK/b.log"; then
    for k in 1 2 3; do "$WORK/b.bin" | tr '\n' ' ' | sed 's/(sink=[^)]*)//'; echo; done
  else echo "NE COMPILE PAS"; grep error "$WORK/b.log" | head -2; fi
done
