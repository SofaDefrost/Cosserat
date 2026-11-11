#!/bin/bash
# Fonctions communes aux scripts de reproduction (à "source"r, pas à exécuter).
# Variables d'environnement (valeurs par défaut pour macOS/Homebrew) :
#   REPO      dépôt git contenant les branches (défaut : dépôt courant)
#   EIGEN_INC dossier d'include d'Eigen3      (défaut : /opt/homebrew/include/eigen3)
#   GTEST_INC / GTEST_LIB  include et lib de googletest (défaut : /opt/homebrew/{include,lib})
#   CXX       compilateur C++20               (défaut : c++)
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="${REPO:-$(git -C "$HERE" rev-parse --show-toplevel)}"
EIGEN_INC="${EIGEN_INC:-/opt/homebrew/include/eigen3}"
GTEST_INC="${GTEST_INC:-/opt/homebrew/include}"
GTEST_LIB="${GTEST_LIB:-/opt/homebrew/lib}"
CXX="${CXX:-c++}"
CXXFLAGS="-std=gnu++20 -O1 -w -I$EIGEN_INC -I$GTEST_INC"

# Dossier de travail jetable (supprimé en sortie du script appelant).
WORK="$(mktemp -d "${TMPDIR:-/tmp}/damien_repro.XXXXXX")"
mkdir -p "$WORK/stub/sofa/testing"
# Les tests d'intégration incluent <sofa/testing/BaseTest.h> : on le remplace par une fixture gtest.
cat > "$WORK/stub/sofa/testing/BaseTest.h" <<'EOT'
#pragma once
#include <gtest/gtest.h>
namespace sofa::testing { class BaseTest : public ::testing::Test {}; }
EOT
WT_LIST=""
cleanup() {
  for w in $WT_LIST; do git -C "$REPO" worktree remove --force "$w" >/dev/null 2>&1; done
  git -C "$REPO" worktree prune >/dev/null 2>&1
  rm -rf "$WORK"
}
trap cleanup EXIT

# checkout_ref <ref|dossier> -> affiche le chemin d'un arbre de travail à tester.
#   - un dossier existant (ex. "." ou "$REPO") est utilisé tel quel : on teste l'arbre LOCAL, rien n'est créé ;
#   - sinon, worktree détaché jetable (lecture seule) sur ce ref git.
checkout_ref() {
  if [ -d "$1" ]; then (cd "$1" && pwd -P); return 0; fi
  local ref="$1" dir="$WORK/wt_$(echo "$1" | tr '/:' '__')"
  git -C "$REPO" worktree add --detach "$dir" "$ref" >/dev/null 2>&1 || { echo "ref introuvable: $ref" >&2; return 1; }
  WT_LIST="$WT_LIST $dir"
  echo "$dir"
}

# build_probe <worktree> <source.cpp> <binaire>   (silencieux ; code retour = succès de la compilation)
build_probe() {
  "$CXX" $CXXFLAGS -I"$1/src" "$2" -o "$3" 2>"$3.log"
}
