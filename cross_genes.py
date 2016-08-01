#!/usr/bin/env python3
"""Cross two or more gene sets and return the commong ones."""
from functools import reduce
from itertools import combinations
import os


def common_positions(*filenames):
    """Return a list with the common positions for two TSV files."""
    variants_lists = [load_variants(f) for f in filenames]

    return reduce(cross_two, variants_lists)


def cross_variants(*filenames, exclude=False):
    """Return a dict with the common variants for two TSV files.

    If exclude is True, return the variants unique to the first element instead.

    """
    variants_list = [load_variants(f) for f in filenames]

    if exclude:
        pos = reduce(difference_two, variants_list)
        if len(filenames) > 1:  # The header is not unique and was trimmed
            pos.append("header")
    else:
        pos = reduce(cross_two, variants_list)

    # The common positions HAS TO BE in the first, or any, variant_list.
    return {_: variants_list[0][_] for _ in pos}


def cross_combine(*filenames):
    """Return a dict with all combined two vs two crossings."""
    pair_dicts = {}
    for pair in combinations(filenames, 2):
        # Set the filenames to A-B.
        this_pair = "-".join(
            [os.path.splitext(os.path.basename(_))[0] for _ in pair])
        pair_dicts[this_pair] = cross_multiple_files(*pair)

    return pair_dicts


def cross_two(first, second):
    """Return the common genes between two sets of genes."""
    return list(set(first) & set(second))


def cross_multiple_files(*filenames):
    """Return a list with the common genes between multiple archives."""
    genes_lists = [load_list(f) for f in filenames]

    return reduce(cross_two, genes_lists)


def difference_two(first, second):
    """Return the variants unique to the first set of genes."""
    return list(set(first) - set(second))


def load_list(filename):
    """Return a list with the genes in a filename.

    Some gene names are multiple, like "ABCD, ABCE_1", everyone of them should
     be considered.

    """
    with open(filename, encoding="utf-8", errors="replace") as f1:
        genes = []
        for line in f1:
            gene = line.rstrip().replace('"', "").split(";")
            genes.extend(gene)

    return genes


def load_variants(filename):
    """Return a dict with the variants in a filename."""
    with open(filename, encoding="utf-8", errors="replace") as f1:
        variants = {"header": f1.readline().rstrip().split("\t")}
        for line in f1:
            variant = line.rstrip().split("\t")
            variants[tuple(variant[:3])] = variant

    return variants

def main(args):
    """Perform the CLI command."""
    if args.variants:
        var_bool = cross_variants(*args.filenames, exclude=args.exclusion)
        print("\t".join([_ for _ in var_bool.pop("header")]))
        for v in var_bool.values():
            print("\t".join(v))
    else:
        if args.exclusion:
            print("Not implemented for genes. Remove the --exclusion flag.")
        else:
            # Print the All-vs-All
            print("All-vs-All")
            for gene in cross_multiple_files(*args.filenames):
                print(gene)
            # Print the One-vs-One
            for pair, genes in cross_combine(*args.filenames).items():
                print()
                print(pair)
                for gene in genes:
                    print(gene)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Cross matches two or more files for common lines")
    parser.add_argument("filenames", nargs="+",
                        help="The name of the files to be crossmatched.")
    parser.add_argument("--variants", action="store_true",
                        help="""The files to be used are tsv/tab files with many
                        columns. The three first must be chromosome, start and
                        and end position.""")
    parser.add_argument("--exclusion", action="store_true",
                        help="""Perform an exclusion (items only in the first
                        file) instead a common position search.""")

    args = parser.parse_args()
    main(args)
