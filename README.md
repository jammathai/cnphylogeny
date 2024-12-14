# cnphylogeny

**cnphylogeny** is a small C library for phylogeny inference using copy number
aberrations. Phylogenies are represented as Markov random fields and optimized
with Gibbs sampling.

## Getting Started

To use **cnphylogeny**, simply add `cnphylogeny.c` and `cnphylogeny.h` to your
project. Since **cnphylogeny** uses `<math.h>`, use the `-lm` flag during
compilation.

`cnphylogeny.h` contains detailed documentation in the form of docstrings. The
examples below provide the basic information to get started.

## Examples

### Creating a Node

`cnp_node_new()` creates a binary tree node that stores a copy number profile
(CNP). In the following example, we create a node with no children whose CNP
contains 1,000 bins (note: the bins of the new node's CNP are initialized to
zero):

```C
struct cnp_node *node = cnp_node_new(
    1000, // CNP length
    NULL, // Left child
    NULL // Right child
);
```

Often, we will want to set the new node's CNP contents. For example, we may be
reading CNPs from a file. In that case, we might do something like this:

```C
for (int i = 0; i < node->len; i++)
    node->bins[i] = /** Read a value */;
```

### Creating a Phylogeny

The simplest way to create a phylogeny is to use nested calls to
`cnp_node_new()` to define a tree structure, as shown below:

```C
struct cnp_node *root = cnp_node_new(
    1000,
    cnp_node_new(
        1000,
        cnp_node_new(1000, NULL, NULL),
        cnp_node_new(
            1000,
            cnp_node_new(1000, NULL, NULL),
            cnp_node_new(1000, NULL, NULL)
        )
    ),
    NULL
);
```

Note that, when a node has only one child, it should be the left child. This
allows us to assert that any node with no left child is a leaf node. The
phylogeny created in the example above looks like this:

```mermaid
graph TD;
    root-->A-->B;
    A-->C;
    C-->D;
    C-->E;
```

### Defining a Probability Matrix

`prob_matrix_new()` creates a new probability matrix. Note that probability
matrices must be square, right stochastic matrices. In other words, the number
of rows and columns must be equal and each row must sum to one. Create a
probability matrix like this:

```C
struct prob_matrix *mutation_probs = prob_matrix_new(5, (double []) {
    1, 0, 0, 0, 0,
    0.0025, 0.99, 0.0025, 0.0025, 0.0025,
    0.0025, 0.0025, 0.99, 0.0025, 0.0025,
    0.0025, 0.0025, 0.0025, 0.99, 0.0025,
    0.0025, 0.0025, 0.0025, 0.0025, 0.99,
});
```

The values of a probability matrix are stored as natural log probabilities.
Therefore, be careful when manually assigning probability values after a
probability matrix is created.

### Optimizing a Phylogeny

Optimize a phylogeny using `optimize_phylogeny()` as shown:

```C
phylogeny_optimize(
    root,
    neighbor_probs,
    mutation_probs,
    1000, // Ignore the first 1,000 iterations
    100, // Record every 100th iteration (after the first 1,000)
    50 // Stop after recording 50 iterations
);
```

`optimize_phylogeny()` uses Gibbs sampling to optimize the internal nodes of the
phylogeny (excluding the root). After sampling, the mode value for each bin is
selected. Be aware that `optimize_phylogeny()` modifies values in place.

## Thanks

Thanks to [Palash Sashittal](https://github.com/sashitt2) for guidance and
oversight.
