#!/bin/bash
# usage: compare_refs.sh [ref ...]
# Compile probe_sizeof.cpp et probe_identity.cpp contre chaque ref et affiche les résultats.
# Refs par défaut : le merge-base (master), la branche de Damien, la branche d'audit.
source "$(dirname "$0")/common.sh"
REFS=("$@")
[ ${#REFS[@]} -eq 0 ] && REFS=(fa71b1b defrost/damien-refactoring painless/liegroups-audit)
for ref in "${REFS[@]}"; do
  echo "=================== $ref ($(git -C "$REPO" rev-parse --short "$ref"))"
  wt=$(checkout_ref "$ref") || continue
  for p in probe_sizeof probe_identity; do
    echo "--- $p"
    if build_probe "$wt" "$HERE/$p.cpp" "$WORK/$p.bin"; then "$WORK/$p.bin"; else echo "NE COMPILE PAS :"; grep error "$WORK/$p.bin.log" | head -3; fi
  done
done
