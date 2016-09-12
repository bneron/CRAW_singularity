from test import CRAWTest
from craw.annotation import Entry, Idx, new_entry_type

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
                      'stop':Idx('end', 1)
                      }
        ref_col = 'Position'
        ne = new_entry_type(name, '\t'.join(fields), ref_col, start_col='beg', stop_col='end')
        self.assertTrue(issubclass(ne, Entry))
        self.assertListEqual(fields, ne._fields)
        self.assertDictEqual(fields_idx, ne._fields_idx)

    def test_entry(self):
        pass
    def test_ref(self):
        pass
    def test_start(self):
        pass
    def test_stop(self):
        pass
    def test_str(self):
        pass

class TestAnnotationParser(CRAWTest):
    pass