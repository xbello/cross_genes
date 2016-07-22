"""Cross two or more gene sets and return the commong ones."""
from functools import reduce


def cross_multiple(*gene_sets):
    """Return the common genes between two or more sets of genes."""
    return reduce(cross_two, gene_sets)


def cross_two(first, second):
    """Return the common genes between two sets of genes."""
    return list(set(first) & set(second))


def cross_two_files(filename, filename2):
    """Return a list with the common genes between two archives.

    Some gene names are multiple, like "ABCD, ABCE_1", and only the first one
    should be considered.

    """
    with open(filename) as f1, open(filename2) as f2:
        genes_first = []
        genes_second = []
        for line in f1:
            gene = line.rstrip().replace('"', "").split(";")[0]
            genes_first.append(gene)
        for line in f2:
            gene = line.rstrip().replace('"', "").split(";")[0]
            genes_second.append(gene)

    print(genes_first)

    return cross_two(genes_first, genes_second)
