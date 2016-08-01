# cross_genes
Check the common names in multiple files.

Install
=======

Just download the whole repo and execute `cross_genes.py` file with the python
interpreter. No libraries are needed.

Usage
=====


```
    usage: cross_genes.py [-h] [--variants] [--exclusion]
                          filenames [filenames ...]

    Cross matches two or more files for common lines

    positional arguments:
      filenames    The name of the files to be crossmatched.

    optional arguments:
      -h, --help   show this help message and exit
      --variants   The files to be used are tsv/tab files with many columns. The
                   three first must be chromosome, start and and end position.
      --exclusion  Perform an exclusion (items only in the first file) instead a
                   common position search.
```

Samples
=======

```

python3 cross_genes.py test_files/CASE1.variants.tsv test_files/CASE2.variants.tsv --variants --exclusion
```

Prints a .tsv output with the variants in CASE1 that are unique to that file.

```
python3 cross_genes.py test_files/CASE1.variants.tsv test_files/CASE2.variants.tsv --variants
```

Prints a .tsv output with the variants common to CASE1 and CASE2.

```
python3 cross_genes.py test_files/genes_list2.txt test_files/genes_list3.txt test_files/genes_list.txt
```

Prints a list with the genes common to all the gene_list files and three more
lists with the genes common to each posible pair (combination) of files.
