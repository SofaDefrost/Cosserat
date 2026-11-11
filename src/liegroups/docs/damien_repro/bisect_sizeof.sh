#!/bin/bash
# usage: bisect_sizeof.sh [base] [tip]
# Pour chaque commit de base..tip (défaut : master de Damien -> sa branche), compile les deux sondes
# et affiche sizeof(SE3<double>), polymorphisme et résultat de (g*g^-1).isApprox(Identity).
# Montre QUEL commit introduit le pointeur de table virtuelle.
source "$(dirname "$0")/common.sh"
BASE="${1:-fa71b1b}"; TIP="${2:-defrost/damien-refactoring}"
wt=$(checkout_ref "$BASE") || exit 1
printf "%-8s %-62s %-9s %-6s %s\n" commit sujet "sizeof SE3" "vptr" "g*g^-1==I"
for c in $BASE $(git -C "$REPO" rev-list --reverse "$BASE..$TIP"); do
  git -C "$wt" checkout -q --detach "$c" || continue
  subj=$(git -C "$REPO" log -1 --format=%s "$c" | cut -c1-60)
  sz="?"; vp="?"; id="?"
  if build_probe "$wt" "$HERE/probe_sizeof.cpp" "$WORK/s.bin"; then
    l=$("$WORK/s.bin" | grep '^SE3<double>'); sz=$(echo "$l" | sed 's/.*sizeof=\([0-9]*\).*/\1/'); vp=$(echo "$l" | sed 's/.*polymorphe=//')
  else sz="NO-COMPILE"; fi
  if build_probe "$wt" "$HERE/probe_identity.cpp" "$WORK/i.bin"; then
    id=$("$WORK/i.bin" | grep 'API' | sed 's/.*: //')
  fi
  printf "%-8s %-62s %-9s %-6s %s\n" "$(git -C "$REPO" rev-parse --short "$c")" "$subj" "$sz" "$vp" "$id"
done
