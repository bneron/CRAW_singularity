###########################################################################
#                                                                         #
# This file is part of Counter RNAseq Window (craw) package.              #
#                                                                         #
#    Authors: Bertrand Néron                                              #
#    Copyright © 2017  Institut Pasteur (Paris).                          #
#    see COPYRIGHT file for details.                                      #
#                                                                         #
#    craw is free software: you can redistribute it and/or modify         #
#    it under the terms of the GNU General Public License as published by #
#    the Free Software Foundation, either version 3 of the License, or    #
#    (at your option) any later version.                                  #
#                                                                         #
#    craw is distributed in the hope that it will be useful,              #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of       #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                 #
#    See the GNU General Public License for more details.                 #
#                                                                         #
#    You should have received a copy of the GNU General Public License    #
#    along with craw (see COPYING file).                                  #
#    If not, see <http://www.gnu.org/licenses/>.                          #
#                                                                         #
###########################################################################


import os
import logging
import pandas as pd
import matplotlib.pyplot as plt
from pandas.util.testing import assert_frame_equal

try:
    from tests import CRAWTest
except ImportError as err:
    msg = "Cannot import craw, check your installation or your CRAW_HOME variable : {0!s}".format(err)
    raise ImportError("Cannot import craw, check your installation or your CRAW_HOME variable : {0!s}".format(err))

import craw.heatmap as htmp


class TestHeatmap(CRAWTest):

    @classmethod
    def setUpClass(cls):
        htmp._log.setLevel(logging.ERROR)


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


    def test_sort(self):
        # all kind of test is test is independent unit test below
        # here I just test corner case
        self.assertIsNone(htmp.sort(None, 'by_gene_size'))
        empty_df = pd.DataFrame()
        assert_frame_equal(empty_df, htmp.sort(empty_df, 'by_gene_size'))
        with self.assertRaises(RuntimeError) as ctx:
            data = pd.DataFrame([
                ['S', 'name_d', 'RMD6', 'chrV', '+', 100, 500, 0, 1, 10, 100, 1000],
                ['AS', 'name_c', 'RMD6', 'chrV', '+', 100, 400, 1000, 100, 10, 1, 0],
                ['S', 'name_b', 'DLD3', 'chrV', '+', 100, 300, 10, 100, 10, 100, 10],
                ['AS', 'name_a', 'DLD3', 'chrV', '+', 100, 200, 1000, 100, 1000, 100, 1000]
                ])
            htmp.sort(data, 'foo_bar')
        self.assertEqual("The 'foo_bar' sorting method does not exists.", str(ctx.exception))


    def test_sort_by_gene_size(self):
        data = pd.DataFrame([
            ['S', 'name_d', 'RMD6', 'chrV', '+', 100, 500, 0, 1, 10, 100, 1000],
            ['AS', 'name_c', 'RMD6', 'chrV', '+', 100, 400, 1000, 100, 10, 1, 0],
            ['S', 'name_b', 'DLD3', 'chrV', '+', 100, 300, 10, 100, 10, 100, 10],
            ['AS', 'name_a', 'DLD3', 'chrV', '+', 100, 200, 1000, 100, 1000, 100, 1000]
        ],
            columns=['sense', 'name', 'gene', 'chromosome', 'strand', 'Position', 'gene_stop', '0', '1', '2', '3', '4'],
            index=[0, 1, 2, 3])

        expected_data = pd.DataFrame([
            ['AS', 'name_a', 'DLD3', 'chrV', '+', 100, 200, 1000, 100, 1000, 100, 1000],
            ['S', 'name_b', 'DLD3', 'chrV', '+', 100, 300, 10, 100, 10, 100, 10],
            ['AS', 'name_c', 'RMD6', 'chrV', '+', 100, 400, 1000, 100, 10, 1, 0],
            ['S', 'name_d', 'RMD6', 'chrV', '+', 100, 500, 0, 1, 10, 100, 1000]
        ],
            columns=['sense', 'name', 'gene', 'chromosome', 'strand', 'Position', 'gene_stop', '0', '1', '2', '3', '4'],
            index=[3, 2, 1, 0])

        received_data = htmp._sort_by_gene_size(data, start_col='Position', stop_col='gene_stop')
        assert_frame_equal(expected_data, received_data)
        received_data = htmp._sort_by_gene_size(expected_data, start_col='Position',
                                                stop_col='gene_stop', ascending=False)
        assert_frame_equal(data, received_data)


    def test_sort_using_col(self):
        data = pd.DataFrame([
            ['S', 'name_d', 'RMD6', 'chrV', '+', 14415, 0, 1, 10, 100, 1000],
            ['AS', 'name_c', 'RMD6', 'chrV', '+', 14415, 1000, 100, 10, 1, 0],
            ['S', 'name_b', 'DLD3', 'chrV', '+', 17848, 10, 100, 10,  100, 10],
            ['AS', 'name_a', 'DLD3', 'chrV', '+', 17848, 1000, 100, 1000, 100, 1000]
            ],
            columns=['sense', 'name', 'gene', 'chromosome', 'strand', 'Position', '0', '1', '2', '3', '4'],
            index=[0, 1, 2, 3])

        expected_data = pd.DataFrame([
            ['AS', 'name_a', 'DLD3', 'chrV', '+', 17848, 1000, 100, 1000, 100, 1000],
            ['S', 'name_b', 'DLD3', 'chrV', '+', 17848, 10, 100, 10, 100, 10],
            ['AS', 'name_c', 'RMD6', 'chrV', '+', 14415, 1000, 100, 10, 1, 0],
            ['S', 'name_d', 'RMD6', 'chrV', '+', 14415, 0, 1, 10, 100, 1000],
            ],
            columns=['sense', 'name', 'gene', 'chromosome', 'strand', 'Position', '0', '1', '2', '3', '4'],
            index=[3, 2, 1, 0])

        received_data = htmp._sort_using_col(data, col='name')
        assert_frame_equal(expected_data, received_data)


    def test_sort_using_file(self):
        data = pd.DataFrame([
            ['S', 'name_d', 'RMD6', 'chrV', '+', 14415, 0, 1, 10, 100, 1000],
            ['AS', 'name_c', 'RMD6', 'chrV', '+', 14415, 1000, 100, 10, 1, 0],
            ['S', 'name_b', 'DLD3', 'chrV', '+', 17848, 10, 100, 10,  100, 10],
            ['AS', 'name_a', 'DLD3', 'chrV', '+', 17848, 1000, 100, 1000, 100, 1000]
            ],
            columns=['sense', 'name', 'gene', 'chromosome', 'strand', 'Position', '0', '1', '2', '3', '4'],
            index=[0, 1, 2, 3])

        expected_data = pd.DataFrame([
            ['S', 'name_b', 'DLD3', 'chrV', '+', 17848, 10, 100, 10,  100, 10],
            ['AS', 'name_a', 'DLD3', 'chrV', '+', 17848, 1000, 100, 1000, 100, 1000],
            ['S', 'name_d', 'RMD6', 'chrV', '+', 14415, 0, 1, 10, 100, 1000],
            ['AS', 'name_c', 'RMD6', 'chrV', '+', 14415, 1000, 100, 10, 1, 0],
            ],
            columns=['sense', 'name', 'gene', 'chromosome', 'strand', 'Position', '0', '1', '2', '3', '4'],
            index=['name_b', 'name_a', 'name_d', 'name_c'])
        expected_data.index.name = 'name'
        received_data = htmp._sort_using_file(data, file=os.path.join(self._data_dir, 'sorting_file.txt'))
        assert_frame_equal(expected_data, received_data)


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
        self.assertIsNone(htmp.remove_metadata(None))

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


    def test_lin_norm(self):
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

        received_data = htmp.lin_norm(data)
        assert_frame_equal(expected_data, received_data)
        self.assertIsNone(htmp.lin_norm(None))
        empty_df = pd.DataFrame()
        assert_frame_equal(empty_df, htmp.lin_norm(empty_df))


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


    def test_lin_norm_row_by_row(self):
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

        received_data = htmp.lin_norm_row_by_row(data)
        assert_frame_equal(expected_data, received_data)
        self.assertIsNone(htmp.lin_norm_row_by_row(None))
        empty_df = pd.DataFrame()
        assert_frame_equal(empty_df, htmp.lin_norm_row_by_row(empty_df))


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


class TestMark(CRAWTest):

    @classmethod
    def setUpClass(cls):
        htmp._log.setLevel(logging.ERROR)

    def test_mark(self):
        data = pd.DataFrame([
            ['foo_1', '+', '214',  0, 1, 10, 100, 1000],
            ['foo_2', '-', '142', 1000, 100, 10, 1, 0],
            ['foo_3', '+', '241', 10, 100, 10,  100, 10],
            ['foo_4', '-', '421', 1000, 100, 1000, 100, 1000]
            ],
            columns=['name', 'strand', 'Position', '-1', '0', '1', '2', '3'])
        color_map = plt.cm.get_cmap('Blues')
        pos = 2
        m = htmp.Mark(pos, data, color_map)
        self.assertEqual(m.pos, pos)
        color = tuple([int(round(c * 255)) for c in color_map(1.0)][:-1])
        self.assertEqual(m._color, color)

        m = htmp.Mark(pos, data, color_map, color='red')
        self.assertEqual(m.pos, pos)
        self.assertEqual(m._color, (255, 0, 0))

        with self.assertRaises(ValueError) as ctx:
            htmp.Mark(-5, data, color_map, color='red')
        self.assertEqual("mark position must be -1 >= pos >= 3: provide -5", str(ctx.exception))

        with self.assertRaises(ValueError) as ctx:
            htmp.Mark(pos, data, color_map, color='foo')
        self.assertEqual("foo is not a valid color", str(ctx.exception))

    def test_rgb_int(self):
        data = pd.DataFrame([
            ['foo_1', '+', '214',  0, 1, 10, 100, 1000],
            ['foo_2', '-', '142', 1000, 100, 10, 1, 0],
            ['foo_3', '+', '241', 10, 100, 10,  100, 10],
            ['foo_4', '-', '421', 1000, 100, 1000, 100, 1000]
            ],
            columns=['name', 'strand', 'Position', '-1', '0', '1', '2', '3'])
        color_map = plt.cm.get_cmap('Blues')
        pos = 2
        m = htmp.Mark(pos, data, color_map, color='red')
        self.assertEqual(m.rgb_int(), (255, 0, 0))

    def test_rgb_float(self):
        data = pd.DataFrame([
            ['foo_1', '+', '214', 0, 1, 10, 100, 1000],
            ['foo_2', '-', '142', 1000, 100, 10, 1, 0],
            ['foo_3', '+', '241', 10, 100, 10, 100, 10],
            ['foo_4', '-', '421', 1000, 100, 1000, 100, 1000]
        ],
            columns=['name', 'strand', 'Position', '-1', '0', '1', '2', '3'])
        color_map = plt.cm.get_cmap('Blues')
        pos = 2
        m = htmp.Mark(pos, data, color_map, color='red')
        self.assertEqual(m.rgb_float(), (1.0, 0.0, 0.0))


    def test_get_matrix_bound(self):
        data = pd.DataFrame([
            ['foo_1', '+', '214',  0, 1, 10, 100, 1000],
            ['foo_2', '-', '142', 1000, 100, 10, 1, 0],
            ['foo_3', '+', '241', 10, 100, 10,  100, 10],
            ['foo_4', '-', '421', 1000, 100, 1000, 100, 1000]
            ],
            columns=['name', 'strand', 'Position', '-1', '0', '1', '2', '3'])

        color_map = plt.cm.get_cmap('Blues')
        pos = 2
        m = htmp.Mark(pos, data, color_map)

        mbb = m._get_matrix_bound(data)
        self.assertEqual(mbb, (-1, 3))


    def to_px(self):
        data = pd.DataFrame([
            ['foo_1', '+', '214',  0, 1, 10, 100, 1000],
            ['foo_2', '-', '142', 1000, 100, 10, 1, 0],
            ['foo_3', '+', '241', 10, 100, 10,  100, 10],
            ['foo_4', '-', '421', 1000, 100, 1000, 100, 1000]
            ],
            columns=['name', 'strand', 'Position', '-1', '0', '1', '2', '3'])

        color_map = plt.cm.get_cmap('Blues')
        pos = 2
        m = htmp.Mark(pos, data, color_map)
        self.assertEqual(m.to_px(), 3)
