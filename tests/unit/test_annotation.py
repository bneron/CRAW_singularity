import os
try:
    from tests import CRAWTest
except ImportError as err:
    # TODO
    raise

from craw.annotation import Entry, Idx, new_entry_type, AnnotationParser

class TestEntry(CRAWTest):

    def test_new_entry_type(self):
        name = 'toto'
        fields = ['name', 'gene', 'chromosome', 'strand', 'Position']
        fields_idx = {'ref': Idx('Position', 4),
                      'strand': Idx('strand', 3),
                      }
        ref_col = 'Position'
        ne = new_entry_type(name, '\t'.join(fields), ref_col)
        self.assertTrue(issubclass(ne, Entry))
        self.assertListEqual(fields, ne._fields)
        self.assertDictEqual(fields_idx, ne._fields_idx)

        fields = ['beg', 'end', 'name', 'gene', 'chromosome', 'strand', 'Position']
        fields_idx = {'ref': Idx('Position', 6),
                      'strand': Idx('strand', 5),
                      'start': Idx('beg', 0),
                      'stop': Idx('end', 1)
                      }
        ref_col = 'Position'
        ne_class = new_entry_type(name, '\t'.join(fields), ref_col, start_col='beg', stop_col='end')
        self.assertTrue(issubclass(ne_class, Entry))
        self.assertListEqual(fields, ne_class._fields)
        self.assertDictEqual(fields_idx, ne_class._fields_idx)

    def test_entry(self):
        name = 'toto'
        ref_col = 'Position'
        fields = ['name', 'gene', 'chromosome', 'strand', 'Position']
        ne_class = new_entry_type(name, '\t'.join(fields), ref_col)
        values = ['YEL072W', 'RMD6', 'chrV', '+', 14415]
        ne = ne_class([str(v) for v in values])
        self.assertListEqual(values, ne._values)

        with self.assertRaises(RuntimeError) as ctx:
            extra_values = values[:]
            extra_values.append('extra')
            ne_class(extra_values)
        self.assertEqual(str(ctx.exception),
                         'the number of values does not match with number of fields')

        fields = ['beg', 'end', 'name', 'gene', 'chromosome', 'strand', 'Position']
        ne_class = new_entry_type(name, '\t'.join(fields), ref_col, start_col='beg', stop_col='end')
        values = [14000, 15000, 'YEL072W', 'RMD6', 'chrV', '+', 14415]
        ne = ne_class([str(v) for v in values])
        self.assertListEqual(values, ne._values)

        with self.assertRaises(RuntimeError) as ctx:
            values = [14000, 15000, 'YEL072W', 'RMD6', 'chrV', '=', 14415]
            ne = ne_class([str(v) for v in values])
        self.assertEqual(str(ctx.exception),
                         "strand must be '+ or '-' got '='")


    def test_ref(self):
        name = 'toto'
        ref_col = 'pos_ref'
        fields = ['name', 'gene', 'chromosome', 'strand', 'pos_ref']
        ne_class = new_entry_type(name, '\t'.join(fields), ref_col)
        values = ['YEL072W', 'RMD6', 'chrV', '+', 14415]
        ne = ne_class([str(v) for v in values])
        self.assertEqual(ne.ref, 14415)

    def test_start(self):
        name = 'toto'
        ref_col = 'pos_ref'
        fields = ['beg', 'end', 'name', 'gene', 'chromosome', 'strand', 'pos_ref']
        ne_class = new_entry_type(name, '\t'.join(fields), ref_col, start_col='beg', stop_col='end')
        values = [14000, 15000, 'YEL072W', 'RMD6', 'chrV', '+', 14415]
        ne = ne_class([str(v) for v in values])
        self.assertEqual(ne.start, 14000)

    def test_stop(self):
        name = 'toto'
        ref_col = 'pos_ref'
        fields = ['beg', 'end', 'name', 'gene', 'chromosome', 'strand', 'pos_ref']
        ne_class = new_entry_type(name, '\t'.join(fields), ref_col, start_col='beg', stop_col='end')
        values = [14000, 15000, 'YEL072W', 'RMD6', 'chrV', '+', 14415]
        ne = ne_class([str(v) for v in values])
        self.assertEqual(ne.stop, 15000)

    def test_strand(self):
        name = 'toto'
        ref_col = 'pos_ref'
        fields = ['beg', 'end', 'name', 'gene', 'chromosome', 'strand', 'pos_ref']
        ne_class = new_entry_type(name, '\t'.join(fields), ref_col, start_col='beg', stop_col='end')
        values = [14000, 15000, 'YEL072W', 'RMD6', 'chrV', '+', 14415]
        ne = ne_class([str(v) for v in values])
        self.assertEqual(ne.strand, '+')

        fields = ['beg', 'end', 'name', 'gene', 'chromosome', 'brin', 'pos_ref']
        ne_class = new_entry_type(name, '\t'.join(fields), ref_col, start_col='beg', stop_col='end', strand_col='brin')
        values = [14000, 15000, 'YEL072W', 'RMD6', 'chrV', '+', 14415]
        ne = ne_class([str(v) for v in values])
        self.assertEqual(ne.strand, '+')

    def test_str(self):
        name = 'toto'
        ref_col = 'pos_ref'
        fields = ['beg', 'end', 'name', 'gene', 'chromosome', 'strand', 'pos_ref']
        ne_class = new_entry_type(name, '\t'.join(fields), ref_col, start_col='beg', stop_col='end')
        values = [14000, 15000, 'YEL072W', 'RMD6', 'chrV', '+', 14415]
        ne = ne_class([str(v) for v in values])
        self.assertEqual(str(ne), '\t'.join([str(v) for v in values]))


class TestAnnotationParser(CRAWTest):

    def test_max(self):
        ap = AnnotationParser(os.path.join(self._data_dir, 'annotation_w_start.txt'),
                              'Position',
                              start_col='beg',
                              stop_col='end')
        self.assertEqual(ap.max(), (845, 903))