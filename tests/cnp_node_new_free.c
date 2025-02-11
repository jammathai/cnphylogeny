#include "../include/cnphylogeny.h"

#include <assert.h>


int cnp_len = 10;
copy_num max_copy_num;
double **neighbor_probs;
double **mutation_probs;


int main()
{
    struct cnp_node *node = cnp_node_new(
        0,
        NULL,
        cnp_node_new(
            1,
            (copy_num []) { 2, 2, 3, 3, 2, 1, 1, 1, 2, 2 },
            cnp_node_new(2, NULL, NULL, NULL),
            cnp_node_new(3, NULL, NULL, NULL)
        ),
        NULL
    );

    for (int i = 0; i < cnp_len; i++)
        assert(node->bins[i] == 0);

    assert(node->left->bins[0] == 2);
    assert(node->left->bins[1] == 2);
    assert(node->left->bins[2] == 3);
    assert(node->left->bins[3] == 3);
    assert(node->left->bins[4] == 2);
    assert(node->left->bins[5] == 1);
    assert(node->left->bins[6] == 1);
    assert(node->left->bins[7] == 1);
    assert(node->left->bins[8] == 2);
    assert(node->left->bins[9] == 2);

    assert(node->left);
    assert(node->left->left);
    assert(!node->left->left->left);
    assert(!node->left->left->right);
    assert(node->left->right);
    assert(!node->left->right->left);
    assert(!node->left->right->right);
    assert(!node->right);

    cnp_node_free(node);

    return 0;
}
