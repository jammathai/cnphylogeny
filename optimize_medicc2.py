import csv
import os
import subprocess
import sys
import urllib.request

ptx = sys.argv[1]

base_path = f"data/medicc2/{ptx}"
if not os.path.exists(base_path): os.mkdir(base_path)
if not os.path.exists(f"{base_path}/original"):
    os.mkdir(f"{base_path}/original")
if not os.path.exists(f"{base_path}/result"):
    os.mkdir(f"{base_path}/result")

raw_cnps = f"{base_path}/{ptx}_final_cn_profiles.tsv"
raw_tree = f"{base_path}/{ptx}_final_tree.new"

sample_names = []

def download_medicc2():
    base_url = "https://api.bitbucket.org/2.0/repositories/schwarzlab/medicc2/src/master/examples/output_gundem_et_al_2015"
    urllib.request.urlretrieve(f"{base_url}/{ptx}_final_cn_profiles.tsv", raw_cnps)
    urllib.request.urlretrieve(f"{base_url}/{ptx}_final_tree.new", raw_tree)

def process_cnps():
    with open(raw_cnps, "r") as input:
        rows = list(csv.DictReader(input, delimiter="\t"))
        for i in range(1, 23):
            sample_count = 0
            cnp_len = 0
            chr = f"chr{i}"
            with (
                open(f"{base_path}/original/{chr}_a.csv", "w") as output_a,
                open(f"{base_path}/original/{chr}_b.csv", "w") as output_b,
            ):
                chr_rows = [row for row in rows if row["chrom"] == chr]
                for j, row in enumerate(chr_rows):
                    cnp_len += 1
                    output_a.write(str(row["cn_a"]))
                    output_b.write(str(row["cn_b"]))
                    if (
                        j == len(chr_rows) - 1 or
                        chr_rows[j + 1]["sample_id"] != chr_rows[j]["sample_id"]
                    ):
                        if sample_count == len(sample_names):
                            sample_names.append(row["sample_id"])
                        sample_count += 1
                        output_a.write("\n")
                        output_b.write("\n")
                    else:
                        output_a.write(",")
                        output_b.write(",")

                # add root CNP
                cnp_len //= sample_count
                for i in range(cnp_len):
                    output_a.write("1")
                    output_b.write("1")
                    if (i == cnp_len - 1):
                        output_a.write("\n")
                        output_b.write("\n")
                    else:
                        output_a.write(",")
                        output_b.write(",")

def process_tree():
    with open(raw_tree, "r") as input, open(f"{base_path}/tree.nwk", "w") as output:
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

if __name__ == "__main__":
    download_medicc2()
    process_cnps()
    process_tree()
    subprocess.run(["make"])
    for i in range(1, 22):
        subprocess.run([
            "build/cnphylogeny",
            "-o", f"{base_path}/result/chr{i}_a.csv",
            "-s", "100000",
            f"{base_path}/tree.nwk", f"{base_path}/original/chr{i}_a.csv"
        ])
        subprocess.run([
            "build/cnphylogeny",
            "-o", f"{base_path}/result/chr{i}_b.csv",
            "-s", "100000",
            f"{base_path}/tree.nwk", f"{base_path}/original/chr{i}_b.csv"
        ])
