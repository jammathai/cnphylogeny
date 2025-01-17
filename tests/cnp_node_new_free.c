#include "../include/cnphylogeny.h"

#include <assert.h>


size_t cnp_len = 1000;
copy_num max_copy_num;
double **neighbor_probs;
double **mutation_probs;


int main()
{
    struct cnp_node *node = cnp_node_new(
        cnp_node_new(
            cnp_node_new(NULL, NULL),
            cnp_node_new(NULL, NULL)
        ),
        NULL
    );

    for (int i = 0; i < 1000; i++)
        assert(node->bins[i] == 0);

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
