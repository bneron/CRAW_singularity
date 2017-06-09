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
try:
    from tests import CRAWTest
except ImportError as err:
    msg = "Cannot import craw, check your installation or your CRAW_HOME variable : {0!s}".format(err)
    raise ImportError("Cannot import craw, check your installation or your CRAW_HOME variable : {0!s}".format(err))

from craw.annotation import Entry, Idx, new_entry_type, AnnotationParser


class TestEntry(CRAWTest):

    def test_new_entry_type(self):
        name = 'toto'
        fields = ['name', 'gene', 'chromosome', 'strand', 'Position']
        fields_idx = {'ref': Idx('Position', 4),
                      'strand': Idx('strand', 3),
                      'chr': Idx('chromosome', 2)
                      }
        ref_col = 'Position'
        ne = new_entry_type(name, fields, ref_col)
        self.assertTrue(issubclass(ne, Entry))
        self.assertListEqual(fields, ne._fields)
        self.assertDictEqual(fields_idx, ne._fields_idx)

        fields = ['beg', 'end', 'name', 'gene', 'chromosome', 'strand', 'Position']
        fields_idx = {'ref': Idx('Position', 6),
                      'strand': Idx('strand', 5),
                      'start': Idx('beg', 0),
                      'stop': Idx('end', 1),
                      'chr': Idx('chromosome', 4)
                      }
        ref_col = 'Position'
        ne_class = new_entry_type(name, fields, ref_col, start_col='beg', stop_col='end')
        self.assertTrue(issubclass(ne_class, Entry))
        self.assertListEqual(fields, ne_class._fields)
        self.assertDictEqual(fields_idx, ne_class._fields_idx)

        name = 'toto'
        fields = ['beg', 'end', 'name', 'gene', 'chromosome', 'strand', 'Position']
        fields_idx = {'ref': Idx('Position', 6),
                      'strand': Idx('strand', 5),
                      'start': Idx('beg', 0),
                      'stop': Idx('end', 1),
                      'chr': Idx('chromosome', 4)
                      }
        ref_col = 'foo'
        with self.assertRaises(RuntimeError) as ctx:
            ne_class = new_entry_type(name, fields, ref_col)
        self.assertEqual(str(ctx.exception),
                         "The ref_col 'foo' does not match any fields: 'beg, end, name, gene, chromosome, strand, Position'\n"
                         "You must specify the '--ref-col' option")

        ref_col = 'Position'
        strand_col = 'foo'
        with self.assertRaises(RuntimeError) as ctx:
            ne_class = new_entry_type(name, fields, ref_col, strand_col=strand_col)
        self.assertEqual(str(ctx.exception),
                         "The strand_col 'foo' does not match any fields: 'beg, end, name, gene, chromosome, strand, Position'")

        ref_col = 'Position'
        chr_col = 'foo'
        with self.assertRaises(RuntimeError) as ctx:
            ne_class = new_entry_type(name, fields, ref_col, chr_col=chr_col)
        self.assertEqual(str(ctx.exception),
                         "The chr_col 'foo' does not match any fields: 'beg, end, name, gene, chromosome, strand, Position'\n"
                         "You must specify the '--chr-col' option")

        with self.assertRaises(RuntimeError) as ctx:
            ne_class = new_entry_type(name, fields, ref_col, start_col='foo')
        self.assertEqual(str(ctx.exception),
                         "if start_col is specified stop_col must be specified too and vice versa")

        with self.assertRaises(RuntimeError) as ctx:
            ne_class = new_entry_type(name, fields, ref_col, stop_col='foo')
        self.assertEqual(str(ctx.exception),
                         "if start_col is specified stop_col must be specified too and vice versa")

        with self.assertRaises(RuntimeError) as ctx:
            ne_class = new_entry_type(name, fields, ref_col, start_col='foo', stop_col='end')
        self.assertEqual(str(ctx.exception),
                         "The start_col 'foo' does not match any fields: 'beg, end, name, gene, chromosome, strand, Position'")

        with self.assertRaises(RuntimeError) as ctx:
            ne_class = new_entry_type(name, fields, ref_col, start_col='beg', stop_col='foo')
        self.assertEqual(str(ctx.exception),
                         "The stop_col 'foo' does not match any fields: 'beg, end, name, gene, chromosome, strand, Position'")


    def test_entry(self):
        name = 'toto'
        ref_col = 'Position'
        fields = ['name', 'gene', 'chromosome', 'strand', 'Position']
        ne_class = new_entry_type(name, fields, ref_col)
        values = ['YEL072W', 'RMD6', 'chrV', '+', 14415]
        ne = ne_class([str(v) for v in values])
        self.assertListEqual(values, ne._values)

        with self.assertRaises(RuntimeError) as ctx:
            extra_values = values[:]
            extra_values.append('extra')
            ne_class(extra_values)
        self.assertEqual(str(ctx.exception),
                         "the number of values (6) does not match with number of fields (5): ['YEL072W', 'RMD6', 'chrV', '+', 14415, 'extra']")

        fields = ['beg', 'end', 'name', 'gene', 'chromosome', 'strand', 'Position']
        ne_class = new_entry_type(name, fields, ref_col, start_col='beg', stop_col='end')
        values = [14000, 15000, 'YEL072W', 'RMD6', 'chrV', '+', 14415]
        ne = ne_class([str(v) for v in values])
        self.assertListEqual(values, ne._values)

        # test stop > start on forward strand
        values = [15000, 14000, 'YEL072W', 'RMD6', 'chrV', '+', 14415]
        with self.assertRaises(RuntimeError) as ctx:
            ne = ne_class([str(v) for v in values])
        self.assertEqual(str(ctx.exception),
                         "error in line '15000\t14000\tYEL072W\tRMD6\tchrV\t+\t14415': beg:15000 > end: 14000on forward strand")

        # test switch of start, stop for reverse strand
        values = [15000, 14000, 'YEL072W', 'RMD6', 'chrV', '-', 14415]
        ne = ne_class([str(v) for v in values])
        self.assertListEqual([14000, 15000, 'YEL072W', 'RMD6', 'chrV', '-', 14415], ne._values)

        with self.assertRaises(RuntimeError) as ctx:
            values = [14000, 15000, 'YEL072W', 'RMD6', 'chrV', '=', 14415]
            ne = ne_class([str(v) for v in values])
        self.assertEqual(str(ctx.exception),
                         "strand must be '+/-', '1/-1' or 'for/rev' got '='")

        with self.assertRaises(RuntimeError) as ctx:
            values = [14450, 15000, 'YEL072W', 'RMD6', 'chrV', '+', 14415]
            ne = ne_class([str(v) for v in values])
        self.assertEqual(str(ctx.exception),
                        "error in line '14450\t15000\tYEL072W\tRMD6\tchrV\t+\t14415': Position 14415 is not between beg: 14450 and end: 15000"
        )

        with self.assertRaises(RuntimeError) as ctx:
            values = [14405, 14000, 'YEL072W', 'RMD6', 'chrV', '+', 14415]
            ne = ne_class([str(v) for v in values])
        self.assertEqual(str(ctx.exception),
                        "error in line '14405\t14000\tYEL072W\tRMD6\tchrV\t+\t14415': Position 14415 is not between beg: 14405 and end: 14000"
        )


    def test_chromosome(self):
        name = 'toto'
        ref_col = 'pos_ref'
        fields = ['name', 'gene', 'chr', 'strand', 'pos_ref']
        ne_class = new_entry_type(name, fields, ref_col, chr_col='chr')
        values = ['YEL072W', 'RMD6', 'chrV', '+', 14415]
        ne = ne_class([str(v) for v in values])
        self.assertEqual(ne.chromosome, 'chrV')


    def test_ref(self):
        name = 'toto'
        ref_col = 'pos_ref'
        fields = ['name', 'gene', 'chromosome', 'strand', 'pos_ref']
        ne_class = new_entry_type(name, fields, ref_col)
        values = ['YEL072W', 'RMD6', 'chrV', '+', 14415]
        ne = ne_class([str(v) for v in values])
        self.assertEqual(ne.ref, 14415)


    def test_start(self):
        name = 'toto'
        ref_col = 'pos_ref'
        fields = ['beg', 'end', 'name', 'gene', 'chromosome', 'strand', 'pos_ref']
        ne_class = new_entry_type(name, fields, ref_col, start_col='beg', stop_col='end')
        values = [14000, 15000, 'YEL072W', 'RMD6', 'chrV', '+', 14415]
        ne = ne_class([str(v) for v in values])
        self.assertEqual(ne.start, 14000)
        ne_class = new_entry_type(name, fields, ref_col)
        values = [14000, 15000, 'YEL072W', 'RMD6', 'chrV', '+', 14415]
        ne = ne_class([str(v) for v in values])
        self.assertIsNone(ne.start)


    def test_stop(self):
        name = 'toto'
        ref_col = 'pos_ref'
        fields = ['beg', 'end', 'name', 'gene', 'chromosome', 'strand', 'pos_ref']
        ne_class = new_entry_type(name, fields, ref_col, start_col='beg', stop_col='end')
        values = [14000, 15000, 'YEL072W', 'RMD6', 'chrV', '+', 14415]
        ne = ne_class([str(v) for v in values])
        self.assertEqual(ne.stop, 15000)
        ne_class = new_entry_type(name, fields, ref_col)
        values = [14000, 15000, 'YEL072W', 'RMD6', 'chrV', '+', 14415]
        ne = ne_class([str(v) for v in values])
        self.assertIsNone(ne.stop)


    def test_strand(self):
        name = 'toto'
        ref_col = 'pos_ref'
        fields = ['beg', 'end', 'name', 'gene', 'chromosome', 'strand', 'pos_ref']
        ne_class = new_entry_type(name, fields, ref_col, start_col='beg', stop_col='end')
        values = [14000, 15000, 'YEL072W', 'RMD6', 'chrV', '+', 14415]
        ne = ne_class([str(v) for v in values])
        self.assertEqual(ne.strand, '+')

        name = 'toto'
        ref_col = 'pos_ref'
        fields = ['beg', 'end', 'name', 'gene', 'chromosome', 'strand', 'pos_ref']
        ne_class = new_entry_type(name, fields, ref_col, start_col='beg', stop_col='end')
        values = [14000, 15000, 'YEL072W', 'RMD6', 'chrV', '1', 14415]
        ne = ne_class([str(v) for v in values])
        self.assertEqual(ne.strand, '+')

        name = 'toto'
        ref_col = 'pos_ref'
        fields = ['beg', 'end', 'name', 'gene', 'chromosome', 'strand', 'pos_ref']
        ne_class = new_entry_type(name, fields, ref_col, start_col='beg', stop_col='end')
        values = [14000, 15000, 'YEL072W', 'RMD6', 'chrV', 'for', 14415]
        ne = ne_class([str(v) for v in values])
        self.assertEqual(ne.strand, '+')

        fields = ['beg', 'end', 'name', 'gene', 'chromosome', 'strand', 'pos_ref']
        ne_class = new_entry_type(name, fields, ref_col, start_col='beg', stop_col='end')
        values = [14000, 15000, 'YEL072W', 'RMD6', 'chrV', '-', 14415]
        ne = ne_class([str(v) for v in values])
        self.assertEqual(ne.strand, '-')

        fields = ['beg', 'end', 'name', 'gene', 'chromosome', 'strand', 'pos_ref']
        ne_class = new_entry_type(name, fields, ref_col, start_col='beg', stop_col='end')
        values = [14000, 15000, 'YEL072W', 'RMD6', 'chrV', 'rev', 14415]
        ne = ne_class([str(v) for v in values])
        self.assertEqual(ne.strand, '-')

        fields = ['beg', 'end', 'name', 'gene', 'chromosome', 'strand', 'pos_ref']
        ne_class = new_entry_type(name, fields, ref_col, start_col='beg', stop_col='end')
        values = [14000, 15000, 'YEL072W', 'RMD6', 'chrV', '-1', 14415]
        ne = ne_class([str(v) for v in values])
        self.assertEqual(ne.strand, '-')

        fields = ['beg', 'end', 'name', 'gene', 'chromosome', 'brin', 'pos_ref']
        ne_class = new_entry_type(name, fields, ref_col, start_col='beg', stop_col='end', strand_col='brin')
        values = [14000, 15000, 'YEL072W', 'RMD6', 'chrV', '+', 14415]
        ne = ne_class([str(v) for v in values])
        self.assertEqual(ne.strand, '+')


    def test_str(self):
        name = 'toto'
        ref_col = 'pos_ref'
        fields = ['beg', 'end', 'name', 'gene', 'chromosome', 'strand', 'pos_ref']
        ne_class = new_entry_type(name, fields, ref_col, start_col='beg', stop_col='end')
        values = [14000, 15000, 'YEL072W', 'RMD6', 'chrV', '+', 14415]
        ne = ne_class([str(v) for v in values])
        self.assertEqual(str(ne), '\t'.join([str(v) for v in values]))


    def test_header(self):
        name = 'toto'
        ref_col = 'pos_ref'
        fields = ['transcription start', 'end of gene', 'name', 'gene', 'chromosome', 'strand', 'pos_ref']
        ne_class = new_entry_type(name, fields, ref_col)
        values = [14000, 15000, 'YEL072W', 'RMD6', 'chrV', '+', 14415]
        ne = ne_class([str(v) for v in values])
        self.assertEqual(ne.header, '\t'.join(fields))


    def test_eq(self):
        name = 'toto'
        ref_col = 'pos_ref'
        fields = ['transcription start', 'end of gene', 'name', 'gene', 'chromosome', 'strand', 'pos_ref']
        ne_class = new_entry_type(name, fields, ref_col)
        values = [14000, 15000, 'YEL072W', 'RMD6', 'chrV', '+', 14415]
        ne_1 = ne_class([str(v) for v in values])
        ne_2 = ne_class([str(v) for v in values])
        self.assertEqual(ne_1, ne_2)
        values = [15000, 15000, 'YEL072W', 'RMD6', 'chrV', '+', 14415]
        ne_3 = ne_class([str(v) for v in values])
        self.assertNotEqual(ne_1, ne_3)


class TestAnnotationParser(CRAWTest):

    def test_max(self):
        ap = AnnotationParser(os.path.join(self._data_dir, 'annotation_w_start.txt'),
                              'Position',
                              start_col='beg',
                              stop_col='end')
        self.assertEqual(ap.max(), (10, 10))

        ap = AnnotationParser(os.path.join(self._data_dir, 'annotation_w_start.txt'), 'Position')
        self.assertEqual(ap.max(), (0, 0))


    def test_get_annotations_header_with_start(self):
        ap = AnnotationParser(os.path.join(self._data_dir, 'annotation_w_start.txt'),
                              'Position',
                              start_col='beg',
                              stop_col='end')
        ne_class = new_entry_type('toto',
                                  ['name', 'gene', 'chromosome', 'strand', 'Position', 'beg', 'end'],
                                  'Position',
                                  start_col='beg',
                                  stop_col='end')
        entries = [ne_class(['YEL072W', 'RMD6', 'chrV', '+', '14415', '14405', '14425']),
                   ne_class(['YEL071W', 'DLD3', 'chrV', '+', '17848', '17839', '17853']),
                   ne_class(['snR67', 'SNR67', 'chrV', '+', '61433', '61425', '61439']),
                   ne_class(['YEL043W', 'YEL043W', 'chrV', '+', '73348', '73345', '73350']),
                   ne_class(['YPR036W', 'VMA13', 'chrXVI', '+', '645272', '645270', '645272']),
                   ne_class(['YEL077C', 'YEL077C', 'chrV', '-', '264', '259', '268']),
                   ]
        it = ap.get_annotations()
        for i, e in enumerate(it):
            self.assertEqual(entries[i], e, "\n{}\n{}".format(entries[i], e))


    def test_get_annotations_header_w_spaces(self):
        ap = AnnotationParser(os.path.join(self._data_dir, 'annotation_fields_with_spaces.txt'),
                              'Position',
                              start_col='gene start',
                              stop_col='gene end')
        ne_class = new_entry_type('toto',
                                  ['name', 'gene', 'chromosome', 'strand', 'Position', 'gene start', 'gene end'],
                                  'Position',
                                  start_col='gene start',
                                  stop_col='gene end')
        entries = [ne_class(['YEL072W', 'RMD6', 'chrV', '+', '14415', '14405', '14425']),
                   ne_class(['YEL071W', 'DLD3', 'chrV', '+', '17848', '17839', '17853']),
                   ne_class(['snR67', 'SNR67', 'chrV', '+', '61433', '61425', '61439']),
                   ne_class(['YEL043W', 'YEL043W', 'chrV', '+', '73348', '73345', '73350']),
                   ne_class(['YPR036W', 'VMA13', 'chrXVI', '+', '645272', '645270', '645272']),
                   ]
        it = ap.get_annotations()
        for i, e in enumerate(it):
            self.assertEqual(entries[i], e, "\n{}\n{}".format(entries[i], e))


    def test_annotation_parser_wo_start(self):
        ap = AnnotationParser(os.path.join(self._data_dir, 'annotation_wo_start.txt'), 'Position')

        ne_class = new_entry_type('toto',
                                  ['name', 'gene', 'chromosome', 'strand', 'Position'],
                                  'Position')
        entries = [ne_class(['YEL072W', 'RMD6', 'chrV', '+', '14415']),
                   ne_class(['YEL071W', 'DLD3', 'chrV', '+', '17848']),
                   ne_class(['snR67', 'SNR67', 'chrV', '+', '61433']),
                   ne_class(['YEL043W', 'YEL043W', 'chrV', '+', '73348']),
                   ne_class(['YPR036W', 'VMA13', 'chrXVI', '+', '645272']),
                   ne_class(['YEL077C', 'YEL077C', 'chrV', '-', '262'])
                   ]
        it = ap.get_annotations()
        for i, e in enumerate(it):
            self.assertEqual(entries[i], e, "\n{}\n{}".format(entries[i], e))

        ap = AnnotationParser(os.path.join(self._data_dir, 'annotation_bad_header.txt'),
                              'Position',
                              start_col='beg',
                              stop_col='end')

        it = ap.get_annotations()
        with self.assertRaises(RuntimeError) as ctx:
            it.__next__()
        self.assertEqual(str(ctx.exception),
                         "the number of values (5) does not match with number of fields (7): ['YEL072W', 'RMD6', 'chrV', '+', '14415']")
