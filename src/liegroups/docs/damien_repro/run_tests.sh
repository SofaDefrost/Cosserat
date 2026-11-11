#!/bin/bash
# usage: [ALL=1] run_tests.sh <ref|dossier> [motif]
# Compile et exécute, hors SOFA, les tests d'intégration de liegroups (Tests/integration/*.cpp) et
# unit/test_ComponentOrder.cpp pour un ref git. Affiche par test : PASSED/FAILED, ou la cause d'échec
# de compilation / d'édition de liens. Un motif optionnel filtre les noms (ex. SGal3).
# <dossier> (ex. "." ) teste l'arbre local sans worktree. ALL=1 ajoute tous les Tests/unit/*.cpp
# (ceux qui exigent SOFA/Cosserat/config.h sont signalés "NE COMPILE PAS").
source "$(dirname "$0")/common.sh"
REF="${1:?usage: run_tests.sh <ref> [motif]}"; PAT="${2:-.}"
wt=$(checkout_ref "$REF") || exit 1
echo "=== $REF ($(git -C "$wt" rev-parse --short HEAD 2>/dev/null || git -C "$REPO" rev-parse --short "$REF"))"
if [ -n "$ALL" ]; then
  files=$( (ls "$wt"/src/liegroups/Tests/integration/*.cpp "$wt"/src/liegroups/Tests/unit/*.cpp) 2>/dev/null | grep -E "$PAT")
else
  files=$( (ls "$wt"/src/liegroups/Tests/integration/*.cpp; ls "$wt"/src/liegroups/Tests/unit/test_ComponentOrder.cpp) 2>/dev/null | grep -E "$PAT")
fi
one() {
  f="$1"; n=$(basename "$f" .cpp); out="$WORK/$n"
  if ! "$CXX" $CXXFLAGS -I"$wt/src" -I"$WORK/stub" "$f" -o "$out.bin" -L"$GTEST_LIB" -lgtest -lgtest_main >"$out.log" 2>&1; then
    if grep -q "Undefined symbols\|undefined reference" "$out.log"; then
      cause=$(grep -m1 -E '^  "' "$out.log" | sed 's/^  //' | cut -c1-90); echo "$n|LIEN IMPOSSIBLE|$cause"
    else
      cause=$(grep -m1 -E 'error' "$out.log" | sed -E 's/^.*(error:|erreur :) //' | cut -c1-90); echo "$n|NE COMPILE PAS|$cause"
    fi
  else
    r=$("$out.bin" 2>&1); p=$(echo "$r" | grep -E '^\[  PASSED' | grep -oE '[0-9]+' | head -1)
    fl=$(echo "$r" | grep -E '^\[  FAILED  \] [0-9]+ test' | grep -oE '[0-9]+' | head -1)
    t=$(echo "$r" | grep -E '^\[==========\] [0-9]+ tests? from .* ran' | grep -oE '[0-9]+' | head -1)
    echo "$n|${t:-0} tests|${p:-0} réussis, ${fl:-0} échoués"
  fi
}
for f in $files; do one "$f" & while [ "$(jobs -rp | wc -l)" -ge 4 ]; do sleep 0.2; done; done | sort | column -t -s'|'
wait
