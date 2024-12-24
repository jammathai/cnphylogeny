#include "../cnphylogeny.h"

#include <assert.h>
#include <string.h>


size_t cnp_len = 5;
copy_num max_copy_num = 2;
double **neighbor_probs;
double **mutation_probs;

static struct cnp_node *root;
static struct cnp_node *interior1;
static struct cnp_node *leaf1;
static struct cnp_node *interior2;
static struct cnp_node *leaf2;
static struct cnp_node *leaf3;


void test_iteration()
{
    memcpy(root->bins, (copy_num []) { 2, 2, 2, 2, 2 }, 5);
    memcpy(interior1->bins, (copy_num []) { 0, 0, 0, 0, 0 }, 5);
    memcpy(leaf1->bins, (copy_num []) { 2, 2, 1, 1, 1 }, 5);
    memcpy(interior2->bins, (copy_num []) { 0, 0, 0, 0, 0 }, 5);
    memcpy(leaf2->bins, (copy_num []) { 2, 0, 2, 2, 2 }, 5);
    memcpy(leaf3->bins, (copy_num []) { 1, 1, 1, 2, 2 }, 5);

    phylogeny_optimize(root, 0, 1, 1);

    assert(!memcmp(interior1->bins, (copy_num []) { 2, 2, 1, 1, 1 }, 5));
    assert(!memcmp(interior2->bins, (copy_num []) { 0, 0, 0, 0, 0 }, 5));

    phylogeny_optimize(root, 0, 1, 1);

    assert(!memcmp(interior1->bins, (copy_num []) { 2, 2, 1, 1, 1 }, 5));
    assert(!memcmp(interior2->bins, (copy_num []) { 2, 1, 1, 2, 2 }, 5));

    phylogeny_optimize(root, 0, 1, 1);

    assert(!memcmp(root->bins, (copy_num []) { 2, 2, 2, 2, 2 }, 5));
    assert(!memcmp(interior1->bins, (copy_num []) { 2, 2, 1, 1, 2 }, 5));
    assert(!memcmp(leaf1->bins, (copy_num []) { 2, 2, 1, 1, 1 }, 5));
    assert(!memcmp(interior2->bins, (copy_num []) { 2, 1, 1, 2, 2 }, 5));
    assert(!memcmp(leaf2->bins, (copy_num []) { 2, 0, 2, 2, 2 }, 5));
    assert(!memcmp(leaf3->bins, (copy_num []) { 1, 1, 1, 2, 2 }, 5));
}

void test_burn_in()
{
    memcpy(root->bins, (copy_num []) { 2, 2, 2, 2, 2 }, 5);
    memcpy(interior1->bins, (copy_num []) { 0, 0, 0, 0, 0 }, 5);
    memcpy(leaf1->bins, (copy_num []) { 2, 2, 1, 1, 1 }, 5);
    memcpy(interior2->bins, (copy_num []) { 0, 0, 0, 0, 0 }, 5);
    memcpy(leaf2->bins, (copy_num []) { 2, 0, 2, 2, 2 }, 5);
    memcpy(leaf3->bins, (copy_num []) { 1, 1, 1, 2, 2 }, 5);

    phylogeny_optimize(root, 2, 1, 1);

    assert(!memcmp(root->bins, (copy_num []) { 2, 2, 2, 2, 2 }, 5));
    assert(!memcmp(interior1->bins, (copy_num []) { 2, 2, 1, 1, 2 }, 5));
    assert(!memcmp(leaf1->bins, (copy_num []) { 2, 2, 1, 1, 1 }, 5));
    assert(!memcmp(interior2->bins, (copy_num []) { 2, 1, 1, 2, 2 }, 5));
    assert(!memcmp(leaf2->bins, (copy_num []) { 2, 0, 2, 2, 2 }, 5));
    assert(!memcmp(leaf3->bins, (copy_num []) { 1, 1, 1, 2, 2 }, 5));
}

void test_mode() {
    memcpy(root->bins, (copy_num []) { 2, 2, 2, 2, 2 }, 5);
    memcpy(interior1->bins, (copy_num []) { 0, 0, 0, 0, 0 }, 5);
    memcpy(leaf1->bins, (copy_num []) { 2, 2, 1, 1, 1 }, 5);
    memcpy(interior2->bins, (copy_num []) { 0, 0, 0, 0, 0 }, 5);
    memcpy(leaf2->bins, (copy_num []) { 2, 0, 2, 2, 2 }, 5);
    memcpy(leaf3->bins, (copy_num []) { 1, 1, 1, 2, 2 }, 5);

    phylogeny_optimize(root, 0, 1, 3);

    assert(!memcmp(root->bins, (copy_num []) { 2, 2, 2, 2, 2 }, 5));
    assert(!memcmp(interior1->bins, (copy_num []) { 2, 2, 1, 1, 1 }, 5));
    assert(!memcmp(leaf1->bins, (copy_num []) { 2, 2, 1, 1, 1 }, 5));
    assert(!memcmp(interior2->bins, (copy_num []) { 2, 1, 1, 2, 2 }, 5));
    assert(!memcmp(leaf2->bins, (copy_num []) { 2, 0, 2, 2, 2 }, 5));
    assert(!memcmp(leaf3->bins, (copy_num []) { 1, 1, 1, 2, 2 }, 5));
}

void test_sample_rate() {
    memcpy(root->bins, (copy_num []) { 2, 2, 2, 2, 2 }, 5);
    memcpy(interior1->bins, (copy_num []) { 0, 0, 0, 0, 0 }, 5);
    memcpy(leaf1->bins, (copy_num []) { 2, 2, 1, 1, 1 }, 5);
    memcpy(interior2->bins, (copy_num []) { 0, 0, 0, 0, 0 }, 5);
    memcpy(leaf2->bins, (copy_num []) { 2, 0, 2, 2, 2 }, 5);
    memcpy(leaf3->bins, (copy_num []) { 1, 1, 1, 2, 2 }, 5);

    phylogeny_optimize(root, 0, 2, 2);

    assert(!memcmp(root->bins, (copy_num []) { 2, 2, 2, 2, 2 }, 5));
    assert(!memcmp(interior1->bins, (copy_num []) { 2, 2, 1, 1, 1 }, 5));
    assert(!memcmp(leaf1->bins, (copy_num []) { 2, 2, 1, 1, 1 }, 5));
    assert(!memcmp(interior2->bins, (copy_num []) { 0, 0, 0, 0, 0 }, 5));
    assert(!memcmp(leaf2->bins, (copy_num []) { 2, 0, 2, 2, 2 }, 5));
    assert(!memcmp(leaf3->bins, (copy_num []) { 1, 1, 1, 2, 2 }, 5));
}


int main()
{
    neighbor_probs = prob_matrix_new((double []) {
        0.9, 0.05, 0.05,
        0.05, 0.9, 0.05,
        0.05, 0.05, 0.9,
    });
    mutation_probs = prob_matrix_new((double []) {
        1, 0, 0,
        0.005, 0.99, 0.005,
        0.005, 0.005, 0.99,
    });

    root = cnp_node_new(
        cnp_node_new(
            cnp_node_new(NULL, NULL),
            cnp_node_new(
                cnp_node_new(NULL, NULL),
                cnp_node_new(NULL, NULL)
            )
        ),
        NULL
    );
    interior1 = root->left;
    leaf1 = root->left->left;
    interior2 = root->left->right;
    leaf2 = root->left->right->left;
    leaf3 = root->left->right->right;

    test_iteration();
    test_burn_in();
    test_mode();
    test_sample_rate();

    free(neighbor_probs);
    free(mutation_probs);
    cnp_node_free(root);

    return 0;
}
