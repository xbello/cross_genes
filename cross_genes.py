#!/usr/bin/env python3
"""Cross two or more gene sets and return the commong ones."""
from functools import reduce
from itertools import combinations
import os
import warnings


def _deprecation(message):
    """Allow DeprecationWarnings to show."""
    warnings.simplefilter('always', DeprecationWarning)
    warnings.warn(message, category=DeprecationWarning)


def common_positions(*filenames):
    """Return a list with the common positions for two TSV files."""
    variants_lists = [load_variants(f) for f in filenames]

    return reduce(cross_two, variants_lists)


def cross_variants(*filenames, **kwargs):
    """Return a dict with the common variants for two TSV files.

    If exclude is True, return the variants unique to the first element.

    """
    variants_list = [load_variants(f) for f in filenames]

    if kwargs.get("exclude"):
        pos = reduce(difference_two, variants_list)
        if len(filenames) > 1:  # The header is not unique and was trimmed
            pos.append("header")
    else:
        pos = reduce(cross_two, variants_list)

    # The common positions HAS TO BE in the first, or any, variant_list.
    return {_: variants_list[0][_] for _ in pos}


def cross_combine(*filenames, **kwargs):
    """Return a dict with all combined two vs two crossings."""
    pair_dicts = {}
    for pair in combinations(filenames, 2):
        # Set the filenames to A-B.
        this_pair = "-".join(
            [os.path.splitext(os.path.basename(_))[0] for _ in pair])
        pair_dicts[this_pair] = cross_multiple_files(*pair)

    return pair_dicts


def cross_combine_variants(*filenames, **kwargs):
    """Return a dict with all combined two vs two crossings of variants."""
    pair_dicts = {}
    for pair in combinations(filenames, 2):
        # Set the filenames to A-B.
        this_pair = "-".join(
            [os.path.splitext(os.path.basename(_))[0] for _ in pair])
        pair_dicts[this_pair] = cross_variants(*pair, **kwargs)

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


def is_variants(filename):
    """Return true if filename has at leas five columns, else return false."""
    with open(filename) as f:
        cols = f.readline().split()
        if len(cols) > 5:
            return True
    return False


def load_list(filename):
    """Return a list with the genes in a filename.

    Some gene names are multiple, like "ABCD, ABCE_1", everyone of them should
     be considered.

    """
    try:
        with open(filename, encoding="utf-8", errors="replace") as f1:
            genes = []
            for line in f1:
                gene = line.rstrip().replace('"', "").split(";")
                genes.extend(gene)
    except TypeError:
        with open(filename) as f1:
            genes = []
            for line in f1:
                gene = line.rstrip().replace('"', "").split(";")
                genes.extend(gene)

    return genes


def load_variants(filename):
    """Return a dict with the variants in a filename."""
    try:
        with open(filename, encoding="utf-8", errors="replace") as f1:
            variants = {"header": f1.readline().rstrip().split("\t")}
            for line in f1:
                variant = line.rstrip().split("\t")
                variants[tuple(variant[:5])] = variant
    except TypeError:
        with open(filename) as f1:
            variants = {"header": f1.readline().rstrip().split("\t")}
            for line in f1:
                variant = line.rstrip().split("\t")
                variants[tuple(variant[:5])] = variant

    return variants


def print_variants(variants_dict):
    """Print the dictionary of variants."""
    for case, variants in variants_dict.items():
        print(case)
        print("\t".join([_ for _ in variants.pop("header")]))
        for variant in variants.values():
            print("\t".join(variant))

def main(args):
    """Perform the CLI command."""
    if all([is_variants(filename) for filename in args.filenames]):
        if len(args.filenames) > 2:
            if args.exclusion:
                print("Only can exclude 2 files. Remove the --exclusion flag.")
            else:
                var_bool = cross_combine_variants(*args.filenames)
                print_variants(var_bool)
        else:
            var_bool = cross_combine_variants(*args.filenames,
                                              exclude=args.exclusion)
            print_variants(var_bool)
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
    parser.add_argument(
        "--variants", action="store_true",
        help="""[Deprecated] The files to be used are tsv/tab files with many
        columns. The five first must be chromosome, start, end and two
        alleles.""")
    parser.add_argument("--exclusion", action="store_true",
                        help="""Perform an exclusion (items only in the first
                        file) instead a common position search.""")

    args = parser.parse_args()

    main(args)

    if args.variants:
        _deprecation("The --variants flag is no longer needed. Any file with" +
                     " five or more columns is assumed to be a variant file.")
