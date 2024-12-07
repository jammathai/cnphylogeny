#include "cnphylogeny.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


static double probs_lookup(struct prob_matrix *probs, copy_num i, copy_num j);
static double ulp(double x);


struct prob_matrix *prob_matrix_new(copy_num order, double *probs)
{
    int len = order * order;
    struct prob_matrix *matrix = malloc(
        sizeof(struct prob_matrix) + len * sizeof(double)
    );
    matrix->order = order;
    for (int i = 0; i < len; i++) matrix->probs[i] = log(probs[i]);
    return matrix;
}


bool prob_matrix_valid(struct prob_matrix *probs)
{
    for (copy_num i = 0; i < probs->order; i++) {
        double sum = 0;
        double epsilon = 0;
        for (copy_num j = 0; j < probs->order; j++) {
            double log_prob = probs_lookup(probs, i, j);
            if (log_prob != log_prob) {
                fprintf(
                    stderr,
                    "Invalid probability matrix: row %u, column %u is probably less than 1\n",
                    i, j
                );
                return false;
            }
            if (log_prob > 0) {
                fprintf(
                    stderr,
                    "Invalid probability matrix: row %u, column %u is greater than 1\n",
                    i, j
                );
                return false;
            }
            sum += exp(log_prob);
            epsilon += ulp(exp(log_prob));
        }
        if (fabs(sum - 1) >= epsilon) {
            fprintf(stderr, "Invalid probability matrix: sum of row %u is not 1\n", i);
            return false;
        }
    }
    return true;
}


static double probs_lookup(struct prob_matrix *probs, copy_num i, copy_num j)
{
    return probs->probs[i * probs->order + j];
}


static double ulp(double x)
{
    return fabs(nextafter(x, x > 0 ? INFINITY : -INFINITY) - x);
}
