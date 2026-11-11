#!/bin/bash
# Étape 2 du plan d'intégration (docs/INTEGRATION_DAMIEN_REFACTORING.md) : renommages mécaniques de Damien,
# rejoués par script sur NOTRE arbre au lieu d'une fusion textuelle (0 conflit).
#
#   RealSpace   -> Rn            (fichiers, classe, alias Rn3d, RnTest, RnTypes)
#   LieGroupBase-> LieGroup      (fichiers .h/.inl, classe, LieGroupTest)
#   _Scalar _Dim _AlgebraDim _ActionDim -> Scalar_ Dim_ AlgebraDim_ ActionDim_   (paramètres de template)
#   concept LieGroup (jamais utilisé) supprimé : il entrerait en collision avec le nom de la classe.
#
# usage : step2_renames.sh [--dry-run]       (à lancer depuis le dépôt)
# Exclus volontairement : docs/ et tous les *.md (repris de Damien à l'étape 1), ce dossier damien_repro.
# NON repris ici : le découpage Rn.h/Rn.inl (f32e897) et la base GroupProperties (887bb92) : étapes suivantes.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
DRY=0; [ "${1:-}" = "--dry-run" ] && DRY=1
run() { if [ $DRY = 1 ]; then echo "[dry-run] $*"; else "$@"; fi; }

# --- 1. renommage de fichiers (git mv conserve l'historique) ---
mv_() { if [ -e "$1" ]; then run git mv "$1" "$2"; else echo "déjà fait ou absent : $1"; fi; }
mv_ src/liegroups/RealSpace.h                               src/liegroups/Rn.h
mv_ src/liegroups/LieGroupBase.h                            src/liegroups/LieGroup.h
mv_ src/liegroups/LieGroupBase.inl                          src/liegroups/LieGroup.inl
mv_ src/liegroups/Tests/integration/RealSpaceTest.cpp       src/liegroups/Tests/integration/RnTest.cpp
mv_ src/liegroups/Tests/integration/LieGroupBaseTest.cpp    src/liegroups/Tests/integration/LieGroupTest.cpp
mv_ src/liegroups/Tests/unit/test_LieGroupBase_Generic.cpp  src/liegroups/Tests/unit/test_LieGroup_Generic.cpp

# --- 2. substitutions de texte dans les fichiers suivis (hors docs, *.md, ce dossier) ---
files=$(git grep -lIE 'RealSpace|LieGroupBase|(^|[^A-Za-z0-9_])_(Scalar|Dim|AlgebraDim|ActionDim)([^A-Za-z0-9_]|$)' \
        -- . ':!docs' ':!*.md' ':!src/liegroups/docs/damien_repro' || true)
echo "$files" | sed 's/^/  texte : /'
if [ $DRY = 0 ]; then
  for f in $files; do
    perl -pi -e '
      s/RealSpace/Rn/g;
      s/LieGroupBase/LieGroup/g;
      s/(?<![A-Za-z0-9_])_Scalar(?![A-Za-z0-9_])/Scalar_/g;
      s/(?<![A-Za-z0-9_])_AlgebraDim(?![A-Za-z0-9_])/AlgebraDim_/g;
      s/(?<![A-Za-z0-9_])_ActionDim(?![A-Za-z0-9_])/ActionDim_/g;
      s/(?<![A-Za-z0-9_])_Dim(?![A-Za-z0-9_])/Dim_/g;
    ' "$f"
  done
fi

# --- 3. le concept LieGroup (jamais utilisé) entre en collision avec la classe LieGroup ---
H=src/liegroups/LieGroup.h
if [ $DRY = 0 ] && [ -f "$H" ] && grep -q 'concept LieGroup = ' "$H"; then
  python3 - "$H" <<'PY'
import re, sys
p = sys.argv[1]; s = open(p).read()
pat = re.compile(r'\n[ \t]*/\*\*\s*\n[ \t]*\* @brief Concept to check if a type is a Lie group\s*\n[ \t]*\*/\s*\n[ \t]*template<typename T>\s*\n[ \t]*concept LieGroup = is_lie_group_v<T>;\s*\n')
s2, n = pat.subn('\n', s)
assert n == 1, "concept LieGroup introuvable (déjà supprimé ?)"
open(p, 'w').write(s2)
print("concept LieGroup supprimé de", p)
PY
fi
echo "terminé$([ $DRY = 1 ] && echo ' (dry-run : rien modifié)')"
