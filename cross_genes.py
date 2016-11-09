#!/usr/bin/env python3
"""Cross two or more gene sets and return the commong ones."""
from functools import reduce
from itertools import combinations
import os
import warnings

import utils


def _deprecation(message):
    """Allow DeprecationWarnings to show."""
    warnings.simplefilter('always', DeprecationWarning)
    warnings.warn(message, category=DeprecationWarning)


def open2(*args, **kwargs):
    """Return a handler to open files both in 2.7 and > 3."""
    try:
        return open(*args, **kwargs)
    except TypeError:
        return open(*args)


def cross_variants(*filenames, **kwargs):
    """Return a dict with the common variants for two TSV files.

    If exclude is True, return the variants unique to the first element.

    """
    variants_list = [load_variants(f, extra=kwargs.get("extra"))
                     for f in filenames]

    if kwargs.get("exclude"):
        pos = reduce(difference_two_genes, variants_list)
        if len(filenames) > 1:  # The header is not unique and was trimmed
            pos.append("header")
    else:
        pos = reduce(cross_two_genes, variants_list)

    # The common positions HAS TO BE in the first, or any, variant_list.
    return {_: variants_list[0][_] for _ in pos}


def cross_combine_genes(*filenames, **kwargs):
    """Return a dict with all combined two vs two crossings."""
    pair_dicts = {}
    for pair in combinations(filenames, 2):
        # Set the filenames to A-B.
        this_pair = "-".join(
            [os.path.splitext(os.path.basename(_))[0] for _ in pair])
        pair_dicts[this_pair] = cross_multiple_genes(*pair)

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


def cross_two_genes(first, second):
    """Return the common genes between two sets of genes."""
    return list(set(first) & set(second))


def cross_multiple_genes(*filenames):
    """Return a list with the common genes between multiple archives."""
    genes_lists = [load_genes(f) for f in filenames]

    return reduce(cross_two_genes, genes_lists)


def difference_two_genes(first, second):
    """Return the variants unique to the first set of genes."""
    return list(set(first) - set(second))


def is_variants(filename):
    """Return true if filename has at least five columns, else return false."""
    with open(filename) as f:
        cols = f.readline().split()
        if len(cols) >= 5:
            return True
    return False


def load_genes(filename):
    """Return a list with the genes in a filename.

    Some gene names are multiple, like "ABCD, ABCE_1", everyone of them should
     be considered.

    """
    with open2(filename, encoding="utf-8", errors="replace") as f1:
        genes = []
        for line in f1:
            gene = line.rstrip().replace('"', "").split(";")
            genes.extend(gene)

    return genes


def load_variants(filename, extra=""):
    """Return a dict with the variants in a filename."""
    extra_columns = []
    if extra:
        for column in extra.split(","):
            extra_columns.append(utils.get_column(filename, column))

    with open2(filename, encoding="utf-8", errors="replace") as f1:
        variants = {"header": f1.readline().rstrip().split("\t")}

        for line in f1:
            variant = line.split("\t")  # Some lines ends in lots of empty cols
                                        #  so we keep every one of them
            variant[-1] = variant[-1].rstrip()

            key = variant[:5]
            for column in extra_columns:
                if variant[column] in ["."]:
                    # Some columns, like "." should be forced to be empty
                    key.append("")
                else:
                    try:
                        key.append(round(float(variant[column]), 2))
                    except ValueError:
                        # If the column cannot be casted to a float, let it in
                        #  as a string
                        key.append(variant[column])

            variants[tuple(key)] = variant

    return variants


def print_genes(*filenames):
    """Print the genes filenames."""
    # Print the All-vs-All
    print("All-vs-All")
    for gene in cross_multiple_genes(*filenames):
        print(gene)
    # Print the One-vs-One
    for pair, genes in cross_combine_genes(*filenames).items():
        print()
        print(pair)
        for gene in genes:
            print(gene)


def print_variants(variants_dict):
    """Print the dictionary of variants."""
    for case, variants in variants_dict.items():
        with open(case + ".tsv", "w") as output:
            output.write("\t".join([_ for _ in variants.pop("header")]))
            output.write("\n")
            for variant in variants.values():
                output.write("\t".join(variant))
                output.write("\n")
        print("Writen {}".format(case + ".tsv"))


def argparser():
    """Return the parsed arguments."""
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
    parser.add_argument(
        "--extra",
        help="""When comparing two variant files only the first five columns are
        checked for equality. With --extra you can add column names separated by
        commas, e.g.: --extra ExAC_ALL,ExAC_NFE

        The columns are converted in floats with two decimals if possible.""")

    return parser


def main(n_args):
    """Perform the CLI command."""
    if all([is_variants(filename) for filename in n_args.filenames]):
        if len(n_args.filenames) > 2:
            if n_args.exclusion:
                raise Exception(
                    "Only can exclude 2 files. Remove the --exclusion flag.")
        var_bool = cross_combine_variants(*n_args.filenames,
                                          exclude=n_args.exclusion,
                                          extra=n_args.extra)
        print_variants(var_bool)
    else:
        if n_args.exclusion:
            raise Exception(
                "Not implemented for genes. Remove the --exclusion flag.")
        else:
            print_genes(*n_args.filenames)

    if n_args.variants:
        _deprecation("The --variants flag is no longer needed. Any file with" +
                     " five or more columns is assumed to be a variant file.")


if __name__ == "__main__":
    import sys
    main(argparser().parse_args(sys.argv[1:]))
