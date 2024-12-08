#include "../cnphylogeny.h"

#include <assert.h>

int main()
{
    struct cnp_node *node = cnp_node_new(
        1000,
        cnp_node_new(
            1000,
            cnp_node_new(1000, NULL, NULL),
            cnp_node_new(1000, NULL, NULL)
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
