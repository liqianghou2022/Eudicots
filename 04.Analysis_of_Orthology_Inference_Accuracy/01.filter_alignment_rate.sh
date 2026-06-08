# Filter entries to keep alignment rate at or above 90%
for f in *.txt; do awk '$3 >= 90' "$f" > "${f%.txt}_filtered.txt"; done