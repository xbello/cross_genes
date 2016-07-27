from unittest import TestCase

import cross_genes as cg

from cross_genes import cross_two
from cross_genes import cross_combine, cross_multiple_files
from cross_genes import load_list


class TestCrossGenes(TestCase):
    def setUp(self):
        self.genes_1 = ["GENE1", "GENE2", "GENE3"]
        self.genes_2 = ["GENE1", "GENE2", "GENE4"]
        self.genes_3 = ["GENE1", "GENE4", "GENE5"]

        self.filename1 = "test_files/genes_list.txt"
        self.filename2 = "test_files/genes_list2.txt"
        self.filename3 = "test_files/genes_list3.txt"

    def test_cross_two(self):
        self.assertCountEqual(cross_two(self.genes_1, self.genes_2),
                              ["GENE1", "GENE2"])
        self.assertCountEqual(cross_two(self.genes_1, self.genes_3),
                              ["GENE1"])

    def test_cross_multiple_files(self):
        self.assertCountEqual(cross_multiple_files(self.filename1,
                                                   self.filename2,
                                                   self.filename3),
                              ["PCDHA9", "PCDHA1", "PCDHA2", "PCDHA8", "NOP16"])

    def test_cross_one_file(self):
        self.assertCountEqual(
            cross_multiple_files(self.filename1),
            ["NOC2L", "SMYD2", "OR2T35", "PCDHA9", "PCDHA1", "PCDHA2", "PCDHA8",
             "NOP16"])

    def test_cross_two_files(self):
        self.assertCountEqual(
            cross_multiple_files(self.filename1, self.filename2),
            ["NOC2L", "PCDHA1", "PCDHA2", "PCDHA8", "PCDHA9", "NOP16"])
        self.assertCountEqual(
            cross_multiple_files(self.filename1, self.filename3),
            ["PCDHA9", "PCDHA1", "PCDHA2", "PCDHA8", "NOP16"])

    def test_cross_combine_multiple_files(self):
        # When comparing multiple files, the user want to know all-vs-all and
        #  one-vs-one individually.
        self.assertCountEqual(
            cross_combine(self.filename1, self.filename2, self.filename3),
            {"genes_list-genes_list2":
             ["NOC2L", "PCDHA1", "PCDHA2", "PCDHA8", "PCDHA9", "NOP16"],
             "genes_list2-genes_list3":
             ["CFAP74", "TRMT13", "SASS6", "PCDHA9", "PCDHA1", "PCDHA2",
              "PCDHA3", "PCDHA4", "PCDHA5", "PCDHA6", "PCDHA7", "PCDHA8",
              "NOP16", "SCARF2"],
             "genes_list-genes_list3":
             ["PCDHA9", "PCDHA1", "PCDHA2", "PCDHA8", "NOP16"]})

    def test_load_list(self):
        self.assertCountEqual(
            load_list(self.filename1),
            ["NOC2L", "SMYD2", "OR2T35", "PCDHA1", "PCDHA2", "PCDHA8", "PCDHA9",
             "NOP16"])


class TestCrossPositions(TestCase):
    def setUp(self):
        self.positions1 = [("chr1", "14930", "14930"),
                           ("chr1", "762592", "762592"),
                           ("chr1", "762601", "762601"),
                           ("chr1", "792263", "792263")]
        self.positions2 = [("chr1", "14930", "14930"),
                           ("chr1", "762273", "762273"),
                           ("chr1", "762601", "762601"),
                           ("chr1", "792263", "792263")]

        self.filename1 = "test_files/CASE1.variants.tsv"
        self.filename2 = "test_files/CASE2.variants.tsv"

    def test_common_positions_are_found(self):
        self.assertCountEqual(
            cg.common_positions(self.filename1, self.filename2),
            ["header",
             ("chr1", "14930", "14930"),
             ("chr1", "762273", "762273"),
             ("chr1", "762592", "762592"),
             ("chr1", "762601", "762601"),
             ("chr1", "792263", "792263")])

    def test_common_variants(self):
        self.assertCountEqual(
            cg.cross_variants(self.filename1, self.filename2).keys(),
            ["header",
             ("chr1", "14930", "14930"),
             ("chr1", "762273", "762273"),
             ("chr1", "762592", "762592"),
             ("chr1", "762601", "762601"),
             ("chr1", "792263", "792263")])

    def test_different_variants(self):
        self.assertCountEqual(
            cg.difference_two(self.positions1, self.positions2),
            [("chr1", "762592", "762592")])
        self.assertCountEqual(
            cg.difference_two(self.positions2, self.positions1),
            [("chr1", "762273", "762273")])

    def test_different_variants_from_files(self):
        self.assertCountEqual(
            cg.cross_variants(
                self.filename2, self.filename1, exclude=True).keys(),
            ["header",
             ("chr1", "792480", "792480"),
             ("chr1", "808922", "808922"),
             ("chr1", "808928", "808928"),
             ("chr1", "871334", "871334"),
             ("chr1", "876499", "876499"),
             ("chr1", "877831", "877831")])

    def test_load_variants(self):
        self.assertCountEqual(cg.load_variants(self.filename1).keys(),
                              ["header",
                               ("chr1", "14907", "14907"),
                               ("chr1", "14930", "14930"),
                               ("chr1", "69511", "69511"),
                               ("chr1", "762273", "762273"),
                               ("chr1", "762589", "762589"),
                               ("chr1", "762592", "762592"),
                               ("chr1", "762601", "762601"),
                               ("chr1", "762632", "762632"),
                               ("chr1", "792263", "792263")])
