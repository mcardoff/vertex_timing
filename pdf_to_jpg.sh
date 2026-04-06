#!/usr/bin/env bash
# pdf_to_jpg.sh
#
# Convert PDF plot files to high-quality JPGs, one file per page.
# Whitespace/empty borders are automatically trimmed.
#
# Usage:
#   ./pdf_to_jpg.sh file1.pdf [file2.pdf ...]
#   ./pdf_to_jpg.sh figs/*.pdf
#
# Output: ~/project/vertex_timing/talk_figs/<basename>[_pN].jpg
#   Single-page PDFs  → <basename>.jpg
#   Multi-page PDFs   → <basename>_p1.jpg, <basename>_p2.jpg, ...

set -euo pipefail

# Ensure Homebrew tools (magick, gs) are on PATH
export PATH="/opt/homebrew/bin:$PATH"

MAGICK=/opt/homebrew/bin/magick
OUT_DIR="${HOME}/project/vertex_timing/talk_figs"
DENSITY=600    # DPI for rasterisation (600 = high-res, good for zoom)
QUALITY=95     # JPEG quality (0-100)
FUZZ=2         # trim fuzz % — catches near-white/near-black border pixels
BORDER=30      # pixels of white padding added back after trimming

# ---------------------------------------------------------------------------

if [[ $# -eq 0 ]]; then
    echo "Usage: $(basename "$0") file1.pdf [file2.pdf ...]"
    exit 1
fi

mkdir -p "$OUT_DIR"

n_converted=0
n_skipped=0

for pdf in "$@"; do
    if [[ ! -f "$pdf" ]]; then
        echo "WARNING: not found — $pdf"
        (( n_skipped++ )) || true
        continue
    fi

    base=$(basename "$pdf" .pdf)

    # Count pages using ImageMagick identify
    n_pages=$("$MAGICK" identify -format "%n\n" "$pdf" 2>/dev/null | head -1)
    n_pages=${n_pages:-1}

    echo "Converting: $pdf  ($n_pages page(s))"

    for (( i=0; i<n_pages; i++ )); do
        if [[ $n_pages -eq 1 ]]; then
            out="$OUT_DIR/${base}.jpg"
        else
            page_num=$(( i + 1 ))
            out="$OUT_DIR/${base}_p${page_num}.jpg"
        fi

        "$MAGICK" \
            -density "$DENSITY" \
            -background white \
            "${pdf}[${i}]" \
            -alpha remove -alpha off \
            -fuzz "${FUZZ}%" -trim +repage \
            -bordercolor white -border "${BORDER}" \
            -quality "$QUALITY" \
            "$out"

        echo "  -> $out"
        (( n_converted++ )) || true
    done
done

echo ""
echo "Done.  Converted: $n_converted file(s)   Skipped: $n_skipped"
echo "Output directory: $OUT_DIR"
