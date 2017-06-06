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
import pysam
import logging
try:
    from tests import CRAWTest
except ImportError as err:
    msg = "Cannot import craw, check your installation or your CRAW_HOME variable : {0!s}".format(err)
    raise ImportError("Cannot import craw, check your installation or your CRAW_HOME variable : {0!s}".format(err))

from craw.wig import Genome, WigParser, _log
from craw.coverage import get_coverage_function, get_wig_coverage, get_bam_coverage
from craw.annotation import new_entry_type


class TestCoverage(CRAWTest):

    @classmethod
    def setUpClass(cls):
        _log.setLevel(logging.ERROR)

    def test_get_coverage_function(self):
        sam_path = os.path.join(self._data_dir, 'small.bam')
        bam_obj = pysam.AlignmentFile(sam_path, "rb")
        func = get_coverage_function(bam_obj)
        self.assertEqual(get_bam_coverage, func)
        genome = Genome()
        func = get_coverage_function(genome)
        self.assertEqual(func, get_wig_coverage)
        with self.assertRaises(RuntimeError) as ctx:
            get_coverage_function('foo')
        self.assertEqual(str(ctx.exception),
                         "get_coverage support only 'wig.Genome' or 'pysam.calignmentfile.AlignmentFile' "
                         "as Input, not str")


    def test_get_bam_coverage_fix_window(self):
        sam_path = os.path.join(self._data_dir, 'small.bam')
        sam_file = pysam.AlignmentFile(sam_path, "rb")
        annot_fields = ['name', 'gene', 'chromosome', 'strand', 'Position']
        entry_cls_name = 'foo'
        ref_col = 'Position'
        ne_class = new_entry_type(entry_cls_name, annot_fields, ref_col)
        value_lines = [['YEL072W', 'RMD6', 'chrV', '+', 14415],
                       ['YEL071W', 'DLD3', 'chrV', '+', 17848],
                       ['YEL071W', 'DLD3', 'chrV', '+', 4]
                       ]

        expected = [{'for': [0, 0, 0, 0, 0, 0, 0, 0],
                    'rev': [0, 0, 0, 0, 0, 0, 0, 0]
                     },
                    {'for': [227, 227, 227, 227, 227, 226, 225, 224],
                     'rev': [0, 0, 0, 0, 0, 0, 0, 0]
                     },
                    {'for': [None, None, 0, 0, 0, 0, 0, 0],
                     'rev': [0, 0, 0, 0, 0, 0, None, None]
                     }
                    ]
        # get_bam_coverage work with 0-based positions
        # whereas annot_entry with 1-based positions
        for values, exp_val in zip(value_lines, expected):
            annot_entry = ne_class([str(v) for v in values])
            forward_cov, reverse_cov = get_bam_coverage(sam_file,
                                                        annot_entry,
                                                        start=values[-1] - 5 - 1,
                                                        stop=values[-1] + 3 - 1,
                                                        qual_thr=0,
                                                        max_left=0,
                                                        max_right=0)
            self.assertListEqual(forward_cov, exp_val['for'])
            self.assertListEqual(reverse_cov, exp_val['rev'])


    def test_get_bam_coverage_strand_rev(self):
        sam_path = os.path.join(self._data_dir, 'small.bam')
        sam_file = pysam.AlignmentFile(sam_path, "rb")
        annot_fields = ['name', 'gene', 'chromosome', 'strand', 'Position']
        entry_cls_name = 'foo'
        ref_col = 'Position'
        ne_class = new_entry_type(entry_cls_name, annot_fields, ref_col)
        value_lines = [['YEL071W', 'DLD3', 'chrV', '+', 17848],
                       ['YEL077C', 'YEL077C', 'chrV', '-', 264]]

        expected = [
                    {'for': [227, 227, 227, 227, 227, 226, 225, 224, 224],
                     'rev': [0, 0, 0, 0, 0, 0, 0, 0, 0]
                     },
                    {'for': [0, 0, 0, 0, 0, 0, 0, 0, 0],
                     'rev': [12, 12, 12, 12, 12, 12, 12, 12, 8]
                    }
                    ]
        # get_bam_coverage work with 0-based positions
        # whereas annot_entry with 1-based positions
        for values, exp_val in zip(value_lines, expected):
            annot_entry = ne_class([str(v) for v in values])
            forward_cov, reverse_cov = get_bam_coverage(sam_file,
                                                        annot_entry,
                                                        start=values[-1] - 5 - 1,
                                                        stop=values[-1] + 3,
                                                        qual_thr=0,
                                                        max_left=0,
                                                        max_right=0)
            self.assertListEqual(forward_cov, exp_val['for'])
            self.assertListEqual(reverse_cov, exp_val['rev'])


    def test_get_bam_coverage_var_window(self):
        sam_path = os.path.join(self._data_dir, 'small.bam')
        sam_file = pysam.AlignmentFile(sam_path, "rb")
        annot_fields = ['name', 'gene', 'chromosome', 'strand', 'Position','beg', 'end']
        entry_cls_name = 'foo'
        ref_col = 'Position'
        ne_class = new_entry_type(entry_cls_name, annot_fields, ref_col, start_col='beg', stop_col='end')
        value_lines = [['YEL072W', 'RMD6', 'chrV', '+', 14415, 14412, 14419],
                       ['YEL071W', 'DLD3', 'chrV', '+', 17848, 17840, 17850]]

        expected = [{'for': [None, None, None, None, 0, 0, 0, 0, 0, 0, 0, None, None],
                     'rev': [None, None, None, None, None, 0, 0, 0, 0, 0, 0, 0, None]
                     },
                    {'for': [227, 227, 227, 227, 227, 227, 227, 227, 226, 225, None, None, None, None],
                     'rev': [None, None, None, None, None, None, None, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                     }
                    ]
        # get_bam_coverage work with 0-based positions
        # whereas annot_entry with 1-based positions
        for values, exp_val in zip(value_lines, expected):
            annot_entry = ne_class([str(v) for v in values])
            forward_cov, reverse_cov = get_bam_coverage(sam_file,
                                                        annot_entry,
                                                        start=annot_entry.start - 1,
                                                        stop=annot_entry.stop - 1,
                                                        qual_thr=0,
                                                        max_left=8,
                                                        max_right=5)
            self.assertListEqual(forward_cov, exp_val['for'])
            self.assertListEqual(reverse_cov, exp_val['rev'])


    def test_get_wig_coverage_fix_window(self):
        wig_parser = WigParser(os.path.join(self._data_dir, 'small_fixed.wig'))
        genome = wig_parser.parse()
        annot_fields = ['name', 'gene', 'chromosome', 'strand', 'Position']
        entry_cls_name = 'foo'
        ref_col = 'Position'
        ne_class = new_entry_type(entry_cls_name, annot_fields, ref_col)
        value_lines = [['YEL072W', 'RMD6', 'chrV', '+', 15],
                       ['YEL071W', 'DLD3', 'chrV', '+', 20],
                       ['YEL071W', 'DLD3', 'chrV', '+', 4]
                       ]
        exp_values = [
            {'for': [0., 150., 150., 150., 150., 150., 100., 100.],
             'rev': [0., 0., 0., 0., 0., 0., 0., 0.]
             },
            {'for': [150., 100., 100., 100., 100., 100., 5., 5.],
             'rev': [0., 0., 0., 0., 0., 0., 0., 0.]
             },
            {'for': [None, None, 0., 0., 0., 0., 0., 0.],
             'rev': [None, None, 0., 0., 0., 0., 0., 0.]
             },
        ]

        # get_wig_coverage work with 0-based positions
        # whereas annot_entry with 1-based positions
        for values, exp_val in zip(value_lines, exp_values):
            annot_entry = ne_class([str(v) for v in values])
            start = values[-1] - 5
            stop = values[-1] + 3
            forward_cov, reverse_cov = get_wig_coverage(genome,
                                                        annot_entry,
                                                        start=start - 1,
                                                        stop=stop - 1,
                                                        qual_thr=0,
                                                        max_left=0,
                                                        max_right=0)

            self.assertListEqual(forward_cov, exp_val['for'])
            self.assertListEqual(reverse_cov, exp_val['rev'])


    def test_get_wig_coverage_var_window(self):
        wig_parser = WigParser(os.path.join(self._data_dir, 'small_variable.wig'))
        genome = wig_parser.parse()
        annot_fields = ['name', 'gene', 'chromosome', 'strand', 'Position', 'beg', 'end']
        entry_cls_name = 'foo'
        ref_col = 'Position'
        ne_class = new_entry_type(entry_cls_name, annot_fields, ref_col, start_col='beg', stop_col='end')
        value_lines = [['YEL072W', 'RMD6', 'chrV', '+', 15, 12, 19],
                       ['YEL071W', 'DLD3', 'chrV', '+', 8, 5, 17]]

        exp_values = [
            {'for': [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, None, None, None, None, None, None, None, None, None],
             'rev': [12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, None, None, None, None, None, None, None, None, None]
             }
            ]

        # get_wig_coverage work with 0-based positions
        # whereas annot_entry with 1-based positions
        for values, exp_val in zip(value_lines, exp_values):
            annot_entry = ne_class([str(v) for v in values])
            forward_cov, reverse_cov = get_wig_coverage(genome,
                                                        annot_entry,
                                                        start=annot_entry.start - 1,
                                                        stop=annot_entry.stop - 1,
                                                        qual_thr=0,
                                                        max_left=3,
                                                        max_right=12)
            self.assertListEqual(forward_cov, exp_val['for'])
            self.assertListEqual(reverse_cov, exp_val['rev'])



