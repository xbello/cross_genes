from os.path import dirname, join
from unittest import TestCase

import utils


class TestColumnOperations(TestCase):
    def test_column_position_is_get_correctly(self):
        path = dirname(__file__)
        annotation = join(path, "test_files/Extra1.tab")
        annotation2 = join(path, "test_files/Extra2.tab")

        self.assertEqual(33, utils.get_column(annotation, "ExAC_ALL"))
        self.assertEqual(34, utils.get_column(annotation2, "ExAC_ALL"))
