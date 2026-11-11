#!/usr/bin/env bash
# One page of a ROOT PDF -> a trimmed PNG for slides (2026-11-11).
#   pdfpage2png.sh <in.pdf> <page> <out.png> [dpi=200]
# Rasterised with Ghostscript, NOT pdftoppm: ROOT's PDFs reference the base-14
# fonts without embedding them, and pdftoppm substitutes a REGULAR face for
# Helvetica-BoldOblique, so the bold-italic "ATLAS" of ATLASLabel comes out
# plain. Ghostscript maps it to NimbusSans-BoldItalic (Helvetica metrics).
set -euo pipefail
in=$1; page=$2; out=$3; dpi=${4:-200}
gs -q -dNOPAUSE -dBATCH -dSAFER -sDEVICE=png16m -r"$dpi" -dTextAlphaBits=4 -dGraphicsAlphaBits=4 \
   -dFirstPage="$page" -dLastPage="$page" -sOutputFile="$out" "$in"
magick "$out" -trim +repage -bordercolor white -border 25 "$out"
