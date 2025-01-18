#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../include/cnphylogeny.h"


const size_t LIST_INIT_SIZE = 256;

size_t cnp_len;
copy_num max_copy_num;
double **mutation_probs;
double **neighbor_probs;

static int burn_in = 1000000;
static int sample_count = 1000000;
static char *mutation_probs_filename = "data/mutation-probs.csv";
static char *neighbor_probs_filename = "data/neighbor-probs.csv";
static char *output;


static void print_usage();
static double **prob_matrix_read(char *filename);


int main(int argc, char **argv)
{
    int opt;
    while ((opt = getopt (argc, argv, "b:c:hm:n:o:")) != -1) {
        switch (opt) {
            case 'b':
                burn_in = strtol(optarg, NULL, 10);
                break;
            case 'c':
                sample_count = strtol(optarg, NULL, 10);
                break;
            case 'h':
                print_usage();
                return EXIT_SUCCESS;
            case 'm':
                mutation_probs_filename = optarg;
                break;
            case 'n':
                neighbor_probs_filename = optarg;
                break;
            default:
                print_usage();
                return EXIT_FAILURE;
        }
    }

    mutation_probs = prob_matrix_read(mutation_probs_filename);
    neighbor_probs = prob_matrix_read(neighbor_probs_filename);

    return EXIT_SUCCESS;
}


static void print_usage()
{
    printf(
        "Usage: cnphylogeny [options] <phylogeny>\n"
        "\n"
        "Arguments:\n"
        "    <phylogeny>  The shared basename of the Newick file and CSV file that\n"
        "                 define a phylogeny\n"
        "\n"
        "Options:\n"
        "    -b <int>        Number of burn-in samples (default: %d)\n"
        "    -c <int>        Number of samples to record (default: %d)\n"
        "    -h              Print this message and exit\n"
        "    -m <csv>        Source mutation probabilities from the specified CSV file\n"
        "                    (default: %s)\n"
        "    -n <csv>        Source neighbor probabilities from the specified CSV file\n"
        "                    (default: %s)\n"
        "    -o <phylogeny>  Write the optimized phylogeny to <phylogeny>.nwk and\n"
        "                    <phylogeny>.csv (default: [YYYY]-[MM]-[DD]T[HH]:[MM]:[SS])\n",
        burn_in,
        sample_count,
        mutation_probs_filename,
        neighbor_probs_filename
    );
}


static double **prob_matrix_read(char *filename)
{
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error: Could not open %s\n", filename);
        exit(EXIT_FAILURE);
    }

    size_t probs_size = LIST_INIT_SIZE;
    int probs_len = 0;
    double *probs = malloc(probs_size * sizeof(double));
    double prob;

    while (fscanf(file, "%lf", &prob) != EOF) {
        getc(file);
        if (probs_len == probs_size) {
            probs_size *= 2;
            probs = realloc(probs, probs_size * sizeof(double));
        }
        probs[probs_len++] = prob;
    }

    fclose(file);

    int order = sqrt(probs_len);
    if (order * order != probs_len) {
        fprintf(
            stderr,
            "Error: Probability matrix %s is not square\n",
            filename
        );
        exit(EXIT_FAILURE);
    }

    if (max_copy_num) {
        if (order - 1 != max_copy_num) {
            fprintf(
                stderr,
                "Error: Probability matrix %s has order %d (expected %d)\n",
                filename,
                order,
                max_copy_num + 1
            );
            exit(EXIT_FAILURE);
        }
    }
    else {
        max_copy_num = order - 1;
    }

    return prob_matrix_new(probs);
}
