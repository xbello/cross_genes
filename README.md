# cross_genes

![python 2.7 3.4](https://img.shields.io/badge/python-2.7%203.4-blue.svg)

Check the common values in multiple files. 

Files can be a list of genes like [this](test_files/genes_list.txt), with one
gene name per line or multiple genes per line splitted with semicolons. Double
quotes are ignored.

Files can be a Tab-Separated-Values like [this](test_files/CASE1.variants.tsv),
with many columns. Only the first five columns are significant, i.e. if two
rows differ only in columns beyond the fifth column they are considered equal.

You can compare:

  - Multiple variant files showing the values *common* between the pairs.
  - Two variant files showing the values *only* in the first one (exclusion).
  - Multiple gene files showing the values *common* to all of them and in pairs.


Install
=======

Just download the whole repo and execute `cross_genes.py` file with the python
interpreter. No extra libraries are needed.

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
      --variants   [Deprecated] The files to be used are tsv/tab files with many
                   columns. The five first must be chromosome, start, end and two
                   alleles.
      --exclusion  Perform an exclusion (items only in the first file) instead a
                   common position search.
```

Samples
=======

Writes .tsv output with the variants in CASE1 that are unique to that file:

```
python3 cross_genes.py test_files/CASE1.variants.tsv test_files/CASE2.variants.tsv --exclusion
```

Writes a .tsv output with the variants common to CASE1 and CASE2:

```
python3 cross_genes.py test_files/CASE1.variants.tsv test_files/CASE2.variants.tsv
```

Prints a list with the genes common to all the gene_list files and three more
lists with the genes common to each posible pair (combination) of files:

```
python3 cross_genes.py test_files/genes_list2.txt test_files/genes_list3.txt test_files/genes_list.txt
```

The same as before, but using filematchers like * or ?:

```
python3 cross_genes.py test_files/*genes_list*.txt
```
