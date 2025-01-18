#include <stdio.h>

#include "../include/cnphylogeny.h"


size_t cnp_len;
copy_num max_copy_num;
double **mutation_probs;
double **neighbor_probs;

static int default_burn_in = 1000000;
static int default_sample_count = 1000000;
static char *default_mutation_probs = "mutation-probs.csv";
static char *default_neighbor_probs = "neighbor-probs.csv";


static void print_usage();


int main()
{
    print_usage();

    return 0;
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
        default_burn_in,
        default_sample_count,
        default_mutation_probs,
        default_neighbor_probs
    );
}
