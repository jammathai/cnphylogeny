#include <stdbool.h>
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../include/cnphylogeny.h"


const int LIST_INIT_SIZE = 256;

int cnp_len;
copy_num max_copy_num;
double **mutation_probs;
double **neighbor_probs;

static int burn_in = 1000000;
static int sample_count = 1000000;
static char *mutation_probs_filename = "data/mutation-probs.csv";
static char *neighbor_probs_filename = "data/neighbor-probs.csv";
static char *output;
static int cnps_len;
static copy_num *cnps;
static struct cnp_node *root;


static double **read_prob_matrix(char *filename);
static struct cnp_node *read_phylogeny(char *name);
static void write_phylogeny(char *name);
static copy_num *read_cnps(char *filename);
static struct cnp_node *parse_newick(char *start, char *end);
static void write_node(struct cnp_node *node, FILE *file);
static void print_phylogeny(
    struct cnp_node *node,
    struct cnp_node *parent,
    char *prefix
);
static void print_node(struct cnp_node *node, struct cnp_node *parent);
static void print_usage();
static FILE *file_open(char *filename, char *modes);


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
            case 'o':
                output = optarg;
                break;
            default:
                print_usage();
                return EXIT_FAILURE;
        }
    }

    mutation_probs = read_prob_matrix(mutation_probs_filename);
    neighbor_probs = read_prob_matrix(neighbor_probs_filename);

    if (optind >= argc) {
        print_usage();
        return EXIT_FAILURE;
    }

    root = read_phylogeny(argv[optind]);

    if (!output) {
        time_t now = time(NULL);
        struct tm *tm = localtime(&now);
        char timestamp[20];
        strftime(timestamp, 20, "%Y-%m-%dT%H:%M:%S", tm);
        output = timestamp;
    }

    printf("%s (before optimization):\n", argv[optind]);
    print_phylogeny(root, NULL, "");

    phylogeny_optimize(root, burn_in, sample_count);

    printf("\n%s (after optimization):\n", output);
    print_phylogeny(root, NULL, "");

    write_phylogeny(output);

    return EXIT_SUCCESS;
}


static double **read_prob_matrix(char *filename)
{
    FILE *file = file_open(filename, "r");

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


static struct cnp_node *read_phylogeny(char *name)
{
    char filename[strlen(name) + 5];
    strncpy(filename, name, strlen(name) + 1);

    strcat(filename, ".csv");
    read_cnps(filename);

    filename[strlen(name)] = '\0';
    strcat(filename, ".nwk");
    FILE *file = file_open(filename, "r");

    fseek(file, 0, SEEK_END);
    long file_size = ftell(file);
    fseek(file, 0, SEEK_SET);
    char newick[file_size + 1];
    fread(newick, 1, file_size, file);
    newick[file_size] = '\0';

    fclose(file);

    return parse_newick(newick, newick + file_size);
}


static void write_phylogeny(char *name)
{
    char filename[strlen(name) + 5];
    strncpy(filename, name, strlen(name) + 1);

    strcat(filename, ".nwk");
    FILE *file = file_open(filename, "w");
    write_node(root, file);
    fclose(file);

    filename[strlen(name)] = '\0';
    strcat(filename, ".csv");
    file = file_open(filename, "w");

    for (int i = 0; i < cnps_len / cnp_len; i++) {
        for (int j = 0; j < cnp_len; j++) {
            fprintf(file, "%hhu", cnps[i * cnp_len + j]);
            if (j < cnp_len - 1) putc(',', file);
            else putc('\n', file);
        }
    }

    fclose(file);
}


static copy_num *read_cnps(char *filename)
{
    FILE *file = file_open(filename, "r");

    int cnps_size = LIST_INIT_SIZE;
    cnps_len = 0;
    int row_count = 0;
    int col_count = 0;
    cnps = malloc(cnps_size * sizeof(copy_num));
    copy_num bin;
    char next_char;

    while (fscanf(file, "%hhu", &bin) != EOF) {
        col_count++;

        next_char = getc(file);
        if (next_char == '\n' || next_char == EOF) {
            row_count++;
            if (cnp_len) {
                if (col_count != cnp_len) {
                    fprintf(
                        stderr,
                        "Error: Row %d of %s has %d columns (expected %d)\n",
                        row_count,
                        filename,
                        col_count,
                        cnp_len
                    );
                    exit(EXIT_FAILURE);
                }
            }
            else {
                cnp_len = col_count;
            }
            col_count = 0;
        }

        if (cnps_len == cnps_size) {
            cnps_size *= 2;
            cnps = realloc(cnps, cnps_size * sizeof(copy_num));
        }

        cnps[cnps_len++] = bin;
    }

    fclose(file);
}


static struct cnp_node *parse_newick(char *start, char *end)
{
    int id;
    struct cnp_node *left = NULL;
    struct cnp_node *right = NULL;

    if (*start == '(') {
        char *i;
        char *comma = NULL;
        int level = 1;

        for (i = start + 1;; i++) {
            if (*i == '(') {
                level++;
            }
            else if (*i == ')') {
                level--;
                if (level == 0) break;
            }
            else if (*i == ',' && level == 1) {
                comma = i;
            }
        }

        if (comma) {
            left = parse_newick(start + 1, comma);
            right = parse_newick(comma + 1, i);
        }
        else {
            left = parse_newick(start + 1, i);
        }

        id = strtol(i + 1, NULL, 10);
    }
    else {
        id = strtol(start, NULL, 10);
    }

    return cnp_node_new(id, cnps + id * cnp_len, left, right);
}


static void write_node(struct cnp_node *node, FILE *file)
{
    if (!node) return;

    if (node->left) {
        putc('(', file);
        write_node(node->left, file);
        if (node->right) {
            putc(',', file);
            write_node(node->right, file);
        }
        putc(')', file);
    }

    fprintf(file, "%d", node->id);

    memcpy(cnps + node->id * cnp_len, node->bins, cnp_len);
}


static void print_phylogeny(
    struct cnp_node *node,
    struct cnp_node *parent,
    char *prefix
)
{
    if (!node) return;

    fputs(prefix, stdout);
    putc('-', stdout);

    print_node(node, parent);

    if (node->right) {
        char left_prefix[strlen(prefix) + 3];
        strcpy(left_prefix, prefix);
        strcat(left_prefix, " |");
        print_phylogeny(node->left, node, left_prefix);

        char right_prefix[strlen(prefix) + 3];
        strcpy(right_prefix, prefix);
        strcat(right_prefix, "  ");
        print_phylogeny(node->right, node, right_prefix);
    }
    else {
        char child_prefix[strlen(prefix) + 3];
        strcpy(child_prefix, prefix);
        strcat(child_prefix, "  ");
        print_phylogeny(node->left, node, child_prefix);
    }
}


static void print_node(struct cnp_node *node, struct cnp_node *parent)
{
    if (!parent) {
        puts(" (Root)");
        return;
    }

    int start = 0;
    int end;
    int i;
    int j;
    bool unchanged = true;

    while (start < cnp_len) {
        i = start;
        j = start;
        while (i < cnp_len - 1 && node->bins[i] == node->bins[i + 1]) i++;
        while (j < cnp_len - 1 && parent->bins[j] == parent->bins[j + 1]) j++;
        end = i < j ? i : j;

        if (node->bins[end] != parent->bins[end]) {
            printf(
                " [%d,%d]:%d->%d",
                start,
                end,
                parent->bins[end],
                node->bins[end]
            );
            unchanged = false;
        }

        start = end + 1;
    }

    if (unchanged) {
        fputs(" (Unchanged)", stdout);
    }

    putc('\n', stdout);
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


static FILE *file_open(char *filename, char *modes)
{
    FILE *file = fopen(filename, modes);
    if (!file) {
        fprintf(stderr, "Error: Could not open %s\n", filename);
        exit(EXIT_FAILURE);
    }
    return file;
}
