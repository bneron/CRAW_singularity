import os
import pandas as pd

from pandas.util.testing import assert_frame_equal

try:
    from tests import CRAWTest
except ImportError as err:
    msg = "Cannot import craw, check your installation or your CRAW_HOME variable : {0!s}".format(err)
    raise ImportError("Cannot import craw, check your installation or your CRAW_HOME variable : {0!s}".format(err))

import craw.heatmap as htmp


class TestHeatmap(CRAWTest):


    def test_get_data(self):
        data_expected = pd.DataFrame([
            ['S', 'name_d', 'RMD6', 'chrV', '+', 14415, 0, 1, 10, 100, 1000],
            ['AS', 'name_c', 'RMD6', 'chrV', '+', 14415, 1000, 100, 10, 1, 0],
            ['S', 'name_b', 'DLD3', 'chrV', '+', 17848, 10, 100, 10,  100, 10],
            ['AS', 'name_a', 'DLD3', 'chrV', '+', 17848, 1000, 100, 1000, 100, 1000]
        ],
            columns=['sense', 'name', 'gene', 'chromosome', 'strand', 'Position', '0', '1', '2', '3', '4'])
        data_received = htmp.get_data(os.path.join(self._data_dir, 'data.cov'))
        assert_frame_equal(data_expected, data_received)


    def test_split_data(self):
        expected_sense = pd.DataFrame([
            ['S', 'name_d', 'RMD6', 'chrV', '+', 14415, 0, 1, 10, 100, 1000],
            ['S', 'name_b', 'DLD3', 'chrV', '+', 17848, 10, 100, 10,  100, 10]
        ],
            columns=['sense', 'name', 'gene', 'chromosome', 'strand', 'Position', '0', '1', '2', '3', '4'],
            index=[0, 2])
        expected_antisense = pd.DataFrame([
            ['AS', 'name_c', 'RMD6', 'chrV', '+', 14415, 1000, 100, 10, 1, 0],
            ['AS', 'name_a', 'DLD3', 'chrV', '+', 17848, 1000, 100, 1000, 100, 1000]
        ],
            columns=['sense', 'name', 'gene', 'chromosome', 'strand', 'Position', '0', '1', '2', '3', '4'],
            index=[1, 3]
        )
        data = htmp.get_data(os.path.join(self._data_dir, 'data.cov'))
        received_sense, received_antisense = htmp.split_data(data)
        assert_frame_equal(expected_sense, received_sense)
        assert_frame_equal(expected_antisense, received_antisense)


    # def test_sort(self):
    #     pass
    #
    #
    # def test_sort_by_gene_size(self):
    #     pass
    #
    #
    # def test_sort_using_col(self):
    #     pass
    #
    #
    # def test_sort_using_file(self):
    #     pass
    #
    #

    def test_remove_metadata(self):
        expected_data = pd.DataFrame([
            [0, 1, 10, 100, 1000],
            [1000, 100, 10, 1, 0],
            [10, 100, 10,  100, 10],
            [1000, 100, 1000, 100, 1000]
            ],
            columns=['0', '1', '2', '3', '4'])
        data = htmp.get_data(os.path.join(self._data_dir, 'data.cov'))
        received_data = htmp.remove_metadata(data)
        assert_frame_equal(expected_data, received_data)


    def test_crop_matrix(self):
        expected_data = pd.DataFrame([
            [1, 10],
            [100, 10],
            [100, 10],
            [100, 1000]
            ],
            columns=['1', '2'])

        data = pd.DataFrame([
            [0, 1, 10, 100, 1000],
            [1000, 100, 10, 1, 0],
            [10, 100, 10,  100, 10],
            [1000, 100, 1000, 100, 1000]
            ],
            columns=['0', '1', '2', '3', '4'])
        received_data = htmp.crop_matrix(data, start_col='1', stop_col='2')
        assert_frame_equal(expected_data, received_data)
        self.assertIsNone(htmp.crop_matrix(None, start_col='1', stop_col='2'))
        empty_df = pd.DataFrame()
        assert_frame_equal(empty_df, htmp.crop_matrix(pd.DataFrame(), start_col='1', stop_col='2'))


    def test_normalize(self):
        data = pd.DataFrame([
            [0, 1, 10, 100, 1000],
            [1000, 100, 10, 1, 0],
            [10, 100, 10,  100, 10],
            [1000, 100, 1000, 100, 1000]
            ],
            columns=['0', '1', '2', '3', '4'])

        expected_data = pd.DataFrame([
            [0, 0.001, 0.01, 0.1, 1],
            [1, 0.1, 0.01, 0.001, 0],
            [0.01, 0.1, 0.01,  0.1, 0.01],
            [1, 0.1, 1, 0.1, 1]
            ],
            columns=['0', '1', '2', '3', '4'])

        received_data = htmp.normalize(data)
        assert_frame_equal(expected_data, received_data)
        self.assertIsNone(htmp.normalize(None))
        empty_df = pd.DataFrame()
        assert_frame_equal(empty_df, htmp.normalize(empty_df))


    def test_log_norm(self):
        data = pd.DataFrame([
            [0, 1, 10, 100, 1000],
            [1000, 100, 10, 1, 0],
            [10, 100, 10,  100, 10],
            [1000, 100, 1000, 100, 1000]
            ],
            columns=['0', '1', '2', '3', '4'])

        expected_data = pd.DataFrame([
            [0.000000,  0.000000,  0.333333,  0.666667,  1.000000],
            [1.000000,  0.666667,  0.333333,  0.000000,  0.000000],
            [0.333333,  0.666667,  0.333333,  0.666667,  0.333333],
            [1.000000,  0.666667,  1.000000,  0.666667,  1.000000],
            ],
            columns=['0', '1', '2', '3', '4'])

        received_data = htmp.log_norm(data)
        assert_frame_equal(expected_data, received_data)
        self.assertIsNone(htmp.log_norm(None))
        empty_df = pd.DataFrame()
        assert_frame_equal(empty_df, htmp.log_norm(empty_df))


    def test_normalize_row_by_row(self):
        data = pd.DataFrame([
            [0, 1, 10, 100, 1000],
            [1000, 100, 10, 1, 0],
            [10, 100, 10,  100, 10],
            [1000, 100, 1000, 100, 1000]
            ],
            columns=['0', '1', '2', '3', '4'])

        expected_data = pd.DataFrame([
            [0.000, 0.001, 0.01, 0.100, 1.0],
            [1.000, 0.100, 0.01, 0.001, 0.0],
            [0.000, 1.000, 0.00, 1.000, 0.0],
            [1.000, 0.000, 1.00, 0.000, 1.0]
            ],
            columns=['0', '1', '2', '3', '4'])

        received_data = htmp.normalize_row_by_row(data)
        assert_frame_equal(expected_data, received_data)
        self.assertIsNone(htmp.normalize_row_by_row(None))
        empty_df = pd.DataFrame()
        assert_frame_equal(empty_df, htmp.normalize_row_by_row(empty_df))


    def test_log_norm_row_by_row(self):
        data = pd.DataFrame([
            [0, 1, 10, 100, 1000],
            [1000, 100, 10, 1, 0],
            [10, 100, 10,  100, 10],
            [1000, 100, 1000, 100, 1000]
            ],
            columns=['0', '1', '2', '3', '4'])

        expected_data = pd.DataFrame([
            [0.000000, 0.000000, 0.333333, 0.666667, 1.000000],
            [1.000000, 0.666667, 0.333333, 0.000000, 0.000000],
            [0.000000, 1.000000, 0.000000, 1.000000, 0.000000],
            [1.000000, 0.000000, 1.000000, 0.000000, 1.000000],
        ],
            columns=['0', '1', '2', '3', '4'])

        received_data = htmp.log_norm_row_by_row(data)
        assert_frame_equal(expected_data, received_data)
        self.assertIsNone(htmp.log_norm_row_by_row(None))
        empty_df = pd.DataFrame()
        assert_frame_equal(empty_df, htmp.log_norm_row_by_row(empty_df))






