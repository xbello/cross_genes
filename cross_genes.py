"""Cross two or more gene sets and return the commong ones."""


def cross_two(first, second):
    """Return the common genes between two sets of genes."""
    return list(set(first) & set(second))
