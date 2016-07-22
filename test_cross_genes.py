from unittest import TestCase

from cross_genes import cross_two


class TestCrossGenes(TestCase):
    def setUp(self):
        self.genes_1 = ["GENE1", "GENE2", "GENE3"]
        self.genes_2 = ["GENE1", "GENE2", "GENE4"]
        self.genes_3 = ["GENE1", "GENE4", "GENE5"]

    def test_cross_two(self):
        self.assertCountEqual(cross_two(self.genes_1, self.genes_2),
                              ["GENE1", "GENE2"])
        self.assertCountEqual(cross_two(self.genes_1, self.genes_3),
                              ["GENE1"])
