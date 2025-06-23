import csv
from DirectedCopyNumberDistance import DirectedCopyNumberDistanceLinear as med

class Node:
    def __init__(self, string):
        self.id = None
        self.left = None
        self.right = None

        if string[0] == "(":
            i = 1
            comma = None
            level = 1

            while True:
                if string[i] == "(":
                    level += 1
                elif string[i] == ")":
                    level -= 1
                    if level == 0: break
                elif string[i] == ',' and level == 1:
                    comma = i;
                i += 1

            if comma:
                self.left = Node(string[1:comma])
                self.right = Node(string[(comma + 1):i])
            else:
                self.left = Node(string[1:i])

            self.id = int(string[(i + 1):]);
        else:
            self.id = int(string);

    def score(self, cnps):
        result = 0
        if self.left:
            self_cnp = cnps[self.id][:]
            left_cnp = cnps[self.left.id][:]
            result += med(self_cnp, left_cnp) + self.left.score(cnps);
        if self.right:
            self_cnp = cnps[self.id][:]
            right_cnp = cnps[self.right.id][:]
            result += med(self_cnp, right_cnp) + self.right.score(cnps);
        return result

def load_cnps(path):
    with open(path, "r") as f:
        reader = csv.reader(f)
        return [
            [int(cn) for cn in row]
            for row in reader
        ]

def score_medicc2(id):
    pad = "0"
    ptx = f"PTX{id:{pad}>3}"
    input_base_path = f"data/medicc2/{ptx}"
    output_base_path = f"data/results/{ptx}"

    with open(f"{input_base_path}/tree.nwk", "r") as input:
        root = Node(input.read())
        original_total = 0
        result_total = 0
        for i in range(1, 23):
            original_cnps_a = load_cnps(f"{input_base_path}/chr{i}a.csv")
            original_cnps_b = load_cnps(f"{input_base_path}/chr{i}b.csv")
            result_cnps_a = load_cnps(f"{output_base_path}/chr{i}a.csv")
            result_cnps_b = load_cnps(f"{output_base_path}/chr{i}b.csv")
            original_score_a = root.score(original_cnps_a)
            original_score_b = root.score(original_cnps_b)
            result_score_a = root.score(result_cnps_a)
            result_score_b = root.score(result_cnps_b)
            original_total += original_score_a + original_score_b
            result_total += result_score_a + result_score_b

            print(f"{ptx},{i},A,{original_score_a},{result_score_a}")
            print(f"{ptx},{i},B,{original_score_b},{result_score_b}")

if __name__ == "__main__":
    print("Phylogeny,Chromosome,Chromatid,Original MED,Result MED")
    for i in range(4, 14):
        if i != 6: score_medicc2(i)
