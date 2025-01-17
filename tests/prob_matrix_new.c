#include "../include/cnphylogeny.h"

#include <assert.h>
#include <math.h>


size_t cnp_len;
copy_num max_copy_num = 1;
double **neighbor_probs;
double **mutation_probs;


int main()
{
    double **probs = prob_matrix_new((double []) {
        1, 0,
        0, 1,
    });

    assert(probs[0][0] == 0);
    assert(probs[0][1] == -INFINITY);
    assert(probs[1][0] == -INFINITY);
    assert(probs[1][1] == 0);

    free(probs);

    return 0;
}
