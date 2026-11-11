#!/bin/bash
# usage: contract_demo.sh   (ne nécessite que Eigen, pas le dépôt)
# Compare, sur une maquette, "virtual final + hooks purs virtuels" (Damien) et "CRTP + static_assert" :
#  - taille des objets ; - message d'erreur quand le dérivé masque inverse() ; - quand il oublie computeInverse().
source "$(dirname "$0")/common.sh"
try() {  # try <titre> <flags...>
  local title="$1"; shift
  if "$CXX" -std=gnu++20 -w -I"$EIGEN_INC" "$@" "$HERE/toy_contract.cpp" -o "$WORK/toy" 2>"$WORK/toy.log"; then
    printf "  %-34s OK      %s\n" "$title" "$("$WORK/toy")"
  else
    printf "  %-34s ERREUR  %s\n" "$title" "$(grep -m1 'error' "$WORK/toy.log" | sed -E "s/^.*error: //; s/^static assertion failed due to requirement '.*': /static_assert: /" | cut -c1-120)"
  fi
}
echo "== CRTP + static_assert (alternative)"
try "contrat respecté" 
try "dérivé masque inverse()"        -DCASE_HIDE
try "dérivé sans computeInverse()"   -DCASE_MISSING
echo "== virtual final + hooks purs virtuels (Damien)"
try "contrat respecté"               -DVIRTUAL_FINAL
try "dérivé masque inverse()"        -DVIRTUAL_FINAL -DCASE_HIDE
try "dérivé sans computeInverse()"   -DVIRTUAL_FINAL -DCASE_MISSING
