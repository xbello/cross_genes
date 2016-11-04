def get_column(filename, column_name):
    """Return the index of the column name."""
    with open(filename) as f:
        for header in f:
            columns = header.rstrip().split("\t")
            return columns.index(column_name)
