from unittest import TestCase

from cross_genes import cross_multiple, cross_two
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

    def test_cross_multiple(self):
        self.assertCountEqual(cross_multiple(self.genes_1,
                                             self.genes_2,
                                             self.genes_3),
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
