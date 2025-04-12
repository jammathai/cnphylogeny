import csv
import os
import subprocess
import sys
import urllib.request

# ensure base path exists
base_path = f"data/medicc2/{sys.argv[2]}"
if not os.path.exists(base_path): os.mkdir(base_path)

raw_cnps_filename = f"{base_path}/{sys.argv[2]}_final_cn_profiles.tsv"
raw_tree_filename = f"{base_path}/{sys.argv[2]}_final_tree.new"
input_basename = f"{base_path}/input"

# download MEDICC2 data
base_url = "https://api.bitbucket.org/2.0/repositories/schwarzlab/medicc2/src/master/examples/output_gundem_et_al_2015"
urllib.request.urlretrieve(f"{base_url}/{sys.argv[2]}_final_cn_profiles.tsv", raw_cnps_filename)
urllib.request.urlretrieve(f"{base_url}/{sys.argv[2]}_final_tree.new", raw_tree_filename)

# process CNP data
sample_names = []
with open(raw_cnps_filename, "r") as input, open(f"{input_basename}.csv", "w") as output:
    cnp_len = 0
    rows = list(csv.DictReader(input, delimiter="\t"))
    for i, row in enumerate(rows):
        if len(sample_names) == 0: cnp_len += 1
        copy_num = 0
        for chromatid in sys.argv[1]: copy_num += int(row["cn_" + chromatid])
        output.write(str(copy_num))
        if i == len(rows) - 1 or rows[i + 1]["sample_id"] != rows[i]["sample_id"]:
            sample_names.append(row["sample_id"])
            output.write("\n")
        else:
            output.write(",")

    # add root CNP
    for i in range(cnp_len):
        output.write("2")
        if (i == cnp_len - 1): output.write("\n")
        else: output.write(",")

# process tree data
with open(raw_tree_filename, "r") as input, open(f"{input_basename}.nwk", "w") as output:
    newick = input.read()

    # remove edge weights
    while (i := newick.find(":")) != -1:
        while newick[i] not in [",", ")", ";"]:
            newick = newick[:i] + newick[i + 1:]

    # replace sample names with indices
    for i, name in enumerate(sample_names):
        newick = newick.replace(name, str(i))

    # add root
    newick = newick.replace(";", str(len(sample_names)))

    output.write(newick)

# make and run cnphylogeny
subprocess.run(["make"])
subprocess.run(["build/cnphylogeny", "-o", f"{base_path}/output", f"{base_path}/input"])
