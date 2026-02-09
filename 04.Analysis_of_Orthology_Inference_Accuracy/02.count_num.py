import csv, os, collections
from typing import Dict, List, Tuple

order_csv = "37-ACEK-202.id.csv"
out_file  = "all.202_summary.csv"
suffix    = ".final.csv"

# Read processing order from CSV
order: List[str] = []
with open(order_csv, newline="", encoding="utf-8") as f:
    r = csv.reader(f)
    for row in r:
        if not row:
            continue
        prefix = row[0].strip()
        if prefix and prefix.lower() != "gene_id":
            order.append(prefix + suffix)

with open(out_file, "w", newline="", encoding="utf-8") as fout:
    w = csv.writer(fout)
    w.writerow(["yes_count", "file", "files"])

    for f in order:
        if not os.path.exists(f):
            w.writerow(["0", f + " (MISSING)", ""])
            continue

        # species_code -> [(gene_id, yes/no)]
        by_sp: Dict[str, List[Tuple[str, str]]] = collections.defaultdict(list)
        file_value = ""

        with open(f, newline="", encoding="utf-8") as fin:
            r = csv.reader(fin)
            header = next(r, None)
            for row in r:
                if not row:
                    continue
                gid = (row[0] if len(row) >= 1 else "").strip()
                sp  = (row[1] if len(row) >= 2 else "").strip()
                yn  = "yes" if (len(row) >= 3 and row[2].strip().lower() == "yes") else "no"
                if sp:
                    by_sp[sp].append((gid or "", yn))

                # Reuse the last column as the files field
                fl = (row[-1] if len(row) >= 1 else "").strip()
                if fl and not file_value:
                    file_value = fl

        # Skip this block if files field is empty
        if not file_value:
            continue

        # Count: one vote per species_code
        yes_count = 0
        for sp, items in by_sp.items():
            # Whether this species_code has at least one "yes"
            any_yes = any(yn == "yes" for _, yn in items)
            if any_yes:
                yes_count += 1

        # Main summary row
        w.writerow([yes_count, f, file_value])

        # DUP expansion: show duplicated details per species_code (for inspection)
        for sp, items in by_sp.items():
            # Deduplicate gene_id for display
            seen = set()
            uniq = []
            for gid, yn in items:
                if gid not in seen:
                    seen.add(gid)
                    uniq.append((gid, yn))
            if len(uniq) >= 2:
                # For readability, show all entries under this species_code
                w.writerow([
                    "DUP",
                    f"{sp}=>" + "|".join(f"{gid or 'NA'}[{yn}]" for gid, yn in uniq),
                    ""
                ])

print(f"Wrote {out_file}")