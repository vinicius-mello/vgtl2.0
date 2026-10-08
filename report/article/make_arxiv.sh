#!/bin/bash
# Builds arxiv.tar.gz: article.tex, article.bbl and the figures it includes,
# nothing else (no .bib, no aux files), then compiles that package in a clean
# directory the way arXiv does (pdflatex, no bibtex) and reports errors and
# undefined references. Run after a full local build, so article.bbl is current.
set -e
cd "$(dirname "$0")"
OUT=arxiv.tar.gz
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
FIGS=$(grep -o '\\includegraphics\[[^]]*\]{[^}]*}' article.tex | sed 's/.*{\(.*\)}/\1/')
mkdir -p "$TMP/pkg/figures"
cp article.tex article.bbl "$TMP/pkg/"
for f in $FIGS; do cp "$f" "$TMP/pkg/$f"; done
(cd "$TMP/pkg" && tar czf "$OLDPWD/$OUT" article.tex article.bbl figures)
mkdir "$TMP/test" && tar xzf "$OUT" -C "$TMP/test"
(cd "$TMP/test" && for i in 1 2 3; do pdflatex -interaction=nonstopmode article.tex > /dev/null || true; done
 echo "files: $(tar tzf "$OLDPWD/$OUT" | grep -v '/$' | tr '\n' ' ')"
 echo "size: $(du -h "$OLDPWD/$OUT" | cut -f1)"
 echo "errors: $(grep -c '^!' article.log)"
 echo "undefined: $(grep -c -i 'undefined' article.log)"
 echo "pages: $(pdfinfo article.pdf | sed -n 's/^Pages: *//p')")
