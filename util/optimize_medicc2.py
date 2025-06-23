import csv
import os
import subprocess
import sys
import urllib.request

def download_medicc2(ptx):
    base_path = f"data/medicc2/{ptx}"
    base_url = "https://api.bitbucket.org/2.0/repositories/schwarzlab/medicc2/src/master/examples/output_gundem_et_al_2015"
    raw_cnps = f"{ptx}_final_cn_profiles.tsv"
    raw_tree = f"{ptx}_final_tree.new"
    urllib.request.urlretrieve(f"{base_url}/{raw_cnps}", f"{base_path}/{raw_cnps}")
    urllib.request.urlretrieve(f"{base_url}/{raw_tree}", f"{base_path}/{raw_tree}")

def process_cnps(ptx):
    base_path = f"data/medicc2/{ptx}"
    raw_cnps = f"{ptx}_final_cn_profiles.tsv"

    sample_names = []

    with open(f"{base_path}/{raw_cnps}", "r") as input:
        rows = list(csv.DictReader(input, delimiter="\t"))
        for i in range(1, 23):
            sample_count = 0
            cnp_len = 0
            chr = f"chr{i}"
            with (
                open(f"{base_path}/{chr}a.csv", "w") as output_a,
                open(f"{base_path}/{chr}b.csv", "w") as output_b,
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

        return sample_names

def process_tree(ptx, sample_names):
    base_path = f"data/medicc2/{ptx}"
    raw_tree = f"{ptx}_final_tree.new"

    with (
        open(f"{base_path}/{raw_tree}", "r") as input,
        open(f"{base_path}/tree.nwk", "w") as output
    ):
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

def optimize_medicc2(id, iterations):
    pad = "0"
    ptx = f"PTX{id:{pad}>3}"

    print(f"{ptx}...")

    input_base_path = f"data/medicc2/{ptx}"
    if not os.path.exists(input_base_path): os.mkdir(input_base_path)
    output_base_path = f"data/results/{ptx}"
    if not os.path.exists(output_base_path): os.mkdir(output_base_path)

    download_medicc2(ptx)
    sample_names = process_cnps(ptx)
    process_tree(ptx, sample_names)
    for i in range(1, 23):
        subprocess.run([
            "build/cnphylogeny",
            "-o", f"{output_base_path}/chr{i}a.csv",
            "-s", str(iterations),
            f"{input_base_path}/tree.nwk", f"{input_base_path}/chr{i}a.csv",
        ])
        subprocess.run([
            "build/cnphylogeny",
            "-o", f"{output_base_path}/chr{i}b.csv",
            "-s", str(iterations),
            f"{input_base_path}/tree.nwk", f"{input_base_path}/chr{i}b.csv",
        ])

if __name__ == "__main__":
    subprocess.run(["make"])
    for i in range(4, 14):
        optimize_medicc2(i, int(sys.argv[1]))
