"""Cross two or more gene sets and return the commong ones."""
from functools import reduce


def cross_multiple(*gene_sets):
    """Return the common genes between two or more sets of genes."""
    return reduce(cross_two, gene_sets)


def cross_two(first, second):
    """Return the common genes between two sets of genes."""
    return list(set(first) & set(second))


def cross_multiple_files(*filenames):
    """Return a list with the common genes between multiple archives."""
    genes_lists = [load_list(f) for f in filenames]

    return reduce(cross_two, genes_lists)


def load_list(filename):
    """Return a list with the genes in a filename.

    Some gene names are multiple, like "ABCD, ABCE_1", and only the first one
      should be considered.

    """
    with open(filename) as f1:
        genes = []
        for line in f1:
            gene = line.rstrip().replace('"', "").split(";")[0]
            genes.append(gene)

    return genes
