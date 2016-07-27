#!/usr/bin/env python3
"""Cross two or more gene sets and return the commong ones."""
from functools import reduce
from itertools import combinations
import os


def common_positions(*filenames):
    """Return a list with the common positions for two TSV files."""
    variants_lists = [load_variants(f) for f in filenames]

    return cross_multiple(*[_.keys() for _ in variants_lists])


def cross_combine(*filenames):
    """Return a dict with all combined two vs two crossings."""
    pair_dicts = {}
    for pair in combinations(filenames, 2):
        # Set the filenames to A-B.
        this_pair = "-".join(
            [os.path.splitext(os.path.basename(_))[0] for _ in pair])
        pair_dicts[this_pair] = cross_multiple_files(*pair)

    return pair_dicts


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

    Some gene names are multiple, like "ABCD, ABCE_1", everyone of them should
     be considered.

    """
    with open(filename) as f1:
        genes = []
        for line in f1:
            gene = line.rstrip().replace('"', "").split(";")
            genes.extend(gene)

    return genes


def load_variants(filename):
    """Return a dict with the variants in a filename."""
    with open(filename) as f1:
        variants = {}
        for line in f1:
            variant = line.rstrip().split("\t")
            variants[tuple(variant[:3])] = variant

    return variants


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Cross matches two or more files for common lines")
    parser.add_argument("filenames", nargs="+",
                        help="The name of the files to be crossmatched.")

    args = parser.parse_args()
    # Print the All-vs-All
    for gene in cross_multiple_files(*args.filenames):
        print(gene)
    # Print the One-vs-One
    for pair, genes in cross_combine(*args.filenames).items():
        print(pair)
        for gene in genes:
            print(gene)
