# gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 \
#    -dPDFSETTINGS=/prepress \
#    -dNOPAUSE -dQUIET -dBATCH \
#    -dDetectDuplicateImages=true \
#    -dCompressFonts=true \
#    -dSubsetFonts=true \
#    -sOutputFile=main_compressed.pdf LatexOut/main.pdf

git add --all
git commit -m "$1"
git push    