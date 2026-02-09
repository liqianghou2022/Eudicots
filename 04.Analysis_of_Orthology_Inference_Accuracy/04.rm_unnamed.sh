# awk versio
awk 'NR>1{ gsub(/Unnamed: [0-9]+/, ""); print }' 37-ACEK-202_transposed.csv > 37-ACEK-202_transposed2.csv
mv 37-ACEK-202_transposed2.csv 37-ACEK-202_transposed.csv