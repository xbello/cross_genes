from os.path import dirname, join
from unittest import TestCase
try:
    from unittest import mock
except ImportError:
    import mock

import cross_genes as cg


class TestWithCountItems(TestCase):
    def assertCountItemsEqual(self, *args, **kwargs):
        """Setup an intermediate for python 2.7 and 3 list testing."""
        try:
            self.assertCountEqual(*args, **kwargs)
        except AttributeError:
            self.assertItemsEqual(*args, **kwargs)


class TestCrossGenes(TestWithCountItems):
    def setUp(self):
        self.genes_1 = ["GENE1", "GENE2", "GENE3"]
        self.genes_2 = ["GENE1", "GENE2", "GENE4"]
        self.genes_3 = ["GENE1", "GENE4", "GENE5"]

        path = dirname(__file__)

        self.filename1 = join(path, "test_files/genes_list.txt")
        self.filename2 = join(path, "test_files/genes_list2.txt")
        self.filename3 = join(path, "test_files/genes_list3.txt")

    def test_cross_two(self):
        self.assertCountItemsEqual(
            cg.cross_two_genes(self.genes_1, self.genes_2), ["GENE1", "GENE2"])
        self.assertCountItemsEqual(
            cg.cross_two_genes(self.genes_1, self.genes_3), ["GENE1"])

    def test_cross_multiple_files(self):
        self.assertCountItemsEqual(
            cg.cross_multiple_genes(self.filename1,
                                    self.filename2,
                                    self.filename3),
            ["PCDHA9", "PCDHA1", "PCDHA2", "PCDHA8", "NOP16"])

    def test_cross_one_file(self):
        self.assertCountItemsEqual(
            cg.cross_multiple_genes(self.filename1),
            ["NOC2L", "SMYD2", "OR2T35", "PCDHA9", "PCDHA1", "PCDHA2",
             "PCDHA8", "NOP16"])

    def test_cross_two_files(self):
        self.assertCountItemsEqual(
            cg.cross_multiple_genes(self.filename1, self.filename2),
            ["NOC2L", "PCDHA1", "PCDHA2", "PCDHA8", "PCDHA9", "NOP16"])
        self.assertCountItemsEqual(
            cg.cross_multiple_genes(self.filename1, self.filename3),
            ["PCDHA9", "PCDHA1", "PCDHA2", "PCDHA8", "NOP16"])

    def test_cross_combine_multiple_files(self):
        # When comparing multiple files, the user want to know all-vs-all and
        #  one-vs-one individually.
        self.assertCountItemsEqual(
            cg.cross_combine_genes(
                self.filename1, self.filename2, self.filename3),
            {"genes_list-genes_list2":
             ["NOC2L", "PCDHA1", "PCDHA2", "PCDHA8", "PCDHA9", "NOP16"],
             "genes_list2-genes_list3":
             ["CFAP74", "TRMT13", "SASS6", "PCDHA9", "PCDHA1", "PCDHA2",
              "PCDHA3", "PCDHA4", "PCDHA5", "PCDHA6", "PCDHA7", "PCDHA8",
              "NOP16", "SCARF2"],
             "genes_list-genes_list3":
             ["PCDHA9", "PCDHA1", "PCDHA2", "PCDHA8", "NOP16"]})

    def test_load_list(self):
        self.assertCountItemsEqual(
            cg.load_genes(self.filename1),
            ["NOC2L", "SMYD2", "OR2T35", "PCDHA1", "PCDHA2", "PCDHA8",
             "PCDHA9", "NOP16"])


class TestCrossPositions(TestWithCountItems):
    def setUp(self):
        self.positions1 = [("chr1", "14930", "14930", "A", "G"),
                           ("chr1", "762592", "762592", "C", "G"),
                           ("chr1", "762601", "762601", "T", "C"),
                           ("chr1", "792263", "792263", "A", "G")]
        self.positions2 = [("chr1", "14930", "14930", "A", "T"),  # Different
                           ("chr1", "762273", "762273", "G", "A"),
                           ("chr1", "762601", "762601", "T", "C"),
                           ("chr1", "792263", "792263", "A", "G")]

        path = dirname(__file__)

        self.filename1 = join(path, "test_files/CASE1.variants.tsv")
        self.filename2 = join(path, "test_files/CASE2.variants.tsv")
        self.filename3 = join(path, "test_files/Extra1.tsv")
        self.filename4 = join(path, "test_files/Extra2.tsv")

    def test_common_variants(self):
        self.assertCountItemsEqual(
            cg.cross_variants(self.filename1, self.filename2).keys(),
            ["header",
             ("chr1", "14930", "14930", "A", "G"),
             ("chr1", "762273", "762273", "G", "A"),
             ("chr1", "762592", "762592", "C", "G"),
             ("chr1", "792263", "792263", "A", "G")])

    def test_different_variants(self):
        self.assertCountItemsEqual(
            cg.difference_two_genes(self.positions1, self.positions2),
            [("chr1", "14930", "14930", "A", "G"),
             ("chr1", "762592", "762592", "C", "G")])
        self.assertCountItemsEqual(
            cg.difference_two_genes(self.positions2, self.positions1),
            [("chr1", "14930", "14930", "A", "T"),
             ("chr1", "762273", "762273", "G", "A")])

    def test_different_variants_from_files(self):
        result = ["header",
                  ("chr1", "762632", "762632", "T", "A"),
                  ("chr1", "792480", "792480", "C", "T"),
                  ("chr1", "808922", "808922", "G", "A"),
                  ("chr1", "808928", "808928", "C", "T"),
                  ("chr1", "871334", "871334", "G", "T"),
                  ("chr1", "876499", "876499", "A", "G"),
                  ("chr1", "877831", "877831", "T", "C")]

        self.assertCountItemsEqual(
            cg.cross_variants(
                self.filename2, self.filename1, exclude=True).keys(),
            result)

    def test_load_variants(self):
        result = ["header",
                  ("chr1", "14907", "14907", "A", "G"),
                  ("chr1", "14930", "14930", "A", "G"),
                  ("chr1", "69511", "69511", "A", "G"),
                  ("chr1", "762273", "762273", "G", "A"),
                  ("chr1", "762589", "762589", "G", "C"),
                  ("chr1", "762592", "762592", "C", "G"),
                  ("chr1", "762601", "762601", "T", "C"),
                  ("chr1", "762632", "762632", "T", "C"),
                  ("chr1", "792263", "792263", "A", "G")]
        self.assertCountItemsEqual(cg.load_variants(self.filename1).keys(),
                                   result)

    def test_load_variants_with_extra_columns(self):
        result = ["header",
                  ("chr1", "14907", "14907", "A", "G", ""),
                  ("chr1", "14930", "14930", "A", "G", ""),
                  ("chr1", "69511", "69511", "A", "G", 0.94),
                  ("chr1", "324822", "324822", "A", "T", 0.07),
                  ("chr1", "664468", "664468", "G", "T", ""),
                  ("chr1", "762273", "762273", "G", "A", 0.81),
                  ("chr1", "762589", "762589", "G", "C", 0.98),
                  ("chr1", "762592", "762592", "C", "G", 0.78),
                  ("chr1", "762601", "762601", "T", "C", 0.78)]
        self.maxDiff = None
        self.assertCountItemsEqual(
            cg.load_variants(self.filename3, extra="ExAC_ALL").keys(), result)

    def test_different_variants_from_files_with_extra_columns(self):
        result = ["header",
                  ("chr1", "762589", "762589", "G", "C", 0.98)]
        self.maxDiff = None

        self.assertCountItemsEqual(
            cg.cross_variants(
                self.filename3,
                self.filename4,
                exclude=True,
                extra="ExAC_ALL").keys(),
            result)


class TestFileDetection(TestCase):
    def test_is_variants(self):
        path = dirname(__file__)
        genes = join(path, "test_files/genes_list.txt")
        variants = join(path, "test_files/CASE1.variants.tsv")

        self.assertTrue(cg.is_variants(variants))
        self.assertFalse(cg.is_variants(genes))


class TestCrossMultipleFiles(TestWithCountItems):
    def setUp(self):
        self.positions1 = [("chr1", "14930", "14930", "A", "G"),
                           ("chr1", "762592", "762592", "C", "G"),
                           ("chr1", "762601", "762601", "T", "C"),
                           ("chr1", "792263", "792263", "A", "G")]
        self.positions2 = [("chr1", "14930", "14930", "A", "T"),  # Different
                           ("chr1", "762273", "762273", "G", "A"),
                           ("chr1", "762601", "762601", "T", "C"),
                           ("chr1", "792263", "792263", "A", "G")]

        path = dirname(__file__)

        self.filename1 = join(path, "test_files/CASE1.variants.tsv")
        self.filename2 = join(path, "test_files/CASE2.variants.tsv")
        self.filename3 = join(path, "test_files/CASE3.variants.tsv")

    def test_common_variants(self):
        self.assertCountItemsEqual(
            cg.cross_variants(
                self.filename1, self.filename2, self.filename3).keys(),
            ["header",
             ("chr1", "14930", "14930", "A", "G"),
             ("chr1", "762273", "762273", "G", "A"),
             ("chr1", "762592", "762592", "C", "G")])

    def test_different_variants(self):
        variants = cg.cross_combine_variants(
            self.filename1, self.filename2, self.filename3)
        self.assertCountItemsEqual(variants.keys(),
                                   ["CASE1.variants-CASE2.variants",
                                    "CASE1.variants-CASE3.variants",
                                    "CASE2.variants-CASE3.variants"])
        self.assertEqual(len(variants["CASE1.variants-CASE2.variants"]), 5)
        self.assertEqual(len(variants["CASE1.variants-CASE3.variants"]), 9)
        self.assertEqual(len(variants["CASE2.variants-CASE3.variants"]), 4)


class TestMainEntry(TestCase):
    def setUp(self):
        path = dirname(__file__)

        self.genes1 = join(path, "test_files/genes_list.txt")
        self.genes2 = join(path, "test_files/genes_list2.txt")
        self.genes3 = join(path, "test_files/genes_list3.txt")

        self.variants1 = join(path, "test_files/CASE1.variants.tsv")
        self.variants2 = join(path, "test_files/CASE2.variants.tsv")

        self.argsparser = cg.argparser()

    @mock.patch("cross_genes.cross_combine_variants")
    @mock.patch("cross_genes.print_variants")
    def test_main_routes_properly_to_variants(self, patched_pv, patched_ccv):
        n_args = self.argsparser.parse_args(args=[
            self.variants1, self.variants2])
        cg.main(n_args)
        # Assert function was called correctly
        patched_ccv.assert_called_once_with(
            self.variants1, self.variants2, exclude=False)

        patched_pv.assert_called_once_with(patched_ccv.return_value)

        # Assert only one file raises an error.
        n_args = self.argsparser.parse_args(
            args=[self.variants1, self.variants1, self.variants1,
                  "--exclusion"])
        with self.assertRaises(Exception):
            cg.main(n_args)

    @mock.patch("cross_genes.print_genes")
    def test_main_routes_properly_to_genes(self, patched_pg):

        n_args = self.argsparser.parse_args(
            args=[self.genes1, self.genes2, self.genes3])

        cg.main(n_args)
        # Assert function was called correctly
        patched_pg.assert_called_once_with(
            self.genes1, self.genes2, self.genes3)

        # Assert exclusion raises an error
        n_args = self.argsparser.parse_args(
            args=[self.genes1, self.genes2, self.genes3,
                  "--exclusion"])
        with self.assertRaises(Exception):
            cg.main(n_args)
