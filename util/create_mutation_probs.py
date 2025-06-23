import csv
import math

def create_mutation_probs(max_copy_num, epsilon):
    mutation_probs = [[1]]
    for i in range(max_copy_num): mutation_probs[0].append(0)
    for i in range(1, max_copy_num + 1):
        mutation_probs.append([0] * (max_copy_num + 1))
        for j in range(max_copy_num + 1):
            if i != j: mutation_probs[i][j] = math.pow(epsilon, abs(i - j))
        mutation_probs[i][i] = 1 - sum(mutation_probs[i])

    with open("data/mutation-probs.csv", "w") as f:
        writer = csv.writer(f)
        writer.writerows(mutation_probs)

if __name__ == "__main__":
    create_mutation_probs(10, 0.000005)
