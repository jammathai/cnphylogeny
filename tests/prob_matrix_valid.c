#include "../cnphylogeny.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>

int main()
{
    struct prob_matrix *probs = prob_matrix_new(
        3,
        (double []) {
            1, 0, 0,
            0.00000000005, 0.9999999999, 0.00000000005,
            0.00000000005, 0.00000000005, 0.9999999999,
        }
    );
    assert(prob_matrix_valid(probs));
    free(probs);

    probs = prob_matrix_new(
        3,
        (double []) {
            1 + DBL_EPSILON, 0, 0,
            0.05, 0.9, 0.05,
            0.05, 0.05, 0.9,
        }
    );
    assert(!prob_matrix_valid(probs));
    free(probs);

    probs = prob_matrix_new(
        3,
        (double []) {
            1, 0, 0,
            0.05, 0.9, 0.05,
            -0.05, 0.05, 0.9,
        }
    );
    assert(!prob_matrix_valid(probs));
    free(probs);

    probs = prob_matrix_new(
        3,
        (double []) {
            1, 0, 0,
            0.05 + DBL_EPSILON, 0.9, 0.05,
            0.05, 0.05, 0.9,
        }
    );
    assert(!prob_matrix_valid(probs));
    free(probs);

    probs = prob_matrix_new(
        3,
        (double []) {
            1, 0, 0,
            0.05, 0.9, 0.05,
            0.05, 0.05, 0.9 - DBL_EPSILON,
        }
    );
    assert(!prob_matrix_valid(probs));
    free(probs);

    return 0;
}
