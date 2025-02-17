import csv
import sys

with open(sys.argv[2], "r") as input, open(sys.argv[3], "w") as output:
    rows = list(csv.DictReader(input, delimiter="\t"))
    for i, row in enumerate(rows):
        copy_num = 0
        for chromatid in sys.argv[1]:
            copy_num += int(row["cn_" + chromatid])
        output.write(str(copy_num))
        if (
            i == len(rows) - 1 or
            rows[i + 1]["sample_id"] != rows[i]["sample_id"]
        ):
            output.write("\n")
        else:
            output.write(",")
