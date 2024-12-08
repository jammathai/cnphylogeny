#include "../cnphylogeny.h"

#include <assert.h>
#include <math.h>

int main()
{
    struct prob_matrix *probs = prob_matrix_new(
        2,
        (double []) {
            1, 0,
            0, 1,
        }
    );

    assert(probs->probs[0] == 0);
    assert(probs->probs[1] == -INFINITY);
    assert(probs->probs[2] == -INFINITY);
    assert(probs->probs[3] == 0);

    return 0;
}
