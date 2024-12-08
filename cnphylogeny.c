#include "cnphylogeny.h"

#include <math.h>
#include <stdlib.h>

struct prob_matrix *prob_matrix_new(copy_num order, double *probs)
{
    int probs_len = order * order;
    struct prob_matrix *matrix = malloc(
        sizeof(struct prob_matrix) + probs_len * sizeof(double)
    );
    matrix->order = order;
    for (int i = 0; i < probs_len; i++) matrix->probs[i] = log(probs[i]);
    return matrix;
}
