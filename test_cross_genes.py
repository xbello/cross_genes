from unittest import TestCase

from cross_genes import cross_multiple, cross_two
from cross_genes import cross_multiple_files, cross_two_files
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

    def test_cross_two_files(self):
        self.assertCountEqual(cross_two_files(self.filename1, self.filename2),
                              ["NOC2L", "PCDHA9", "NOP16"])
        self.assertCountEqual(cross_two_files(self.filename1, self.filename3),
                              ["PCDHA9", "NOP16"])

    def test_cross_multiple_files(self):
        self.assertCountEqual(cross_multiple_files(self.filename1,
                                                   self.filename2,
                                                   self.filename3),
                              ["PCDHA9", "NOP16"])

    def test_load_list(self):
        self.assertCountEqual(
            load_list(self.filename1),
            ["NOC2L", "SMYD2", "OR2T35", "PCDHA9", "NOP16"])
