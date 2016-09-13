from collections import namedtuple


class Entry:

    def __init__(self, values):
        if len(values) != len(self._fields):
            raise RuntimeError("the number of values does not match with number of fields")
        self._values = [self._convert(f, v) for f, v in zip(self._fields, values)]

    def _convert(self, field, value):
        if field == self._fields_idx['ref'].col_name:
            value = int(value)
        elif field == self._fields_idx['strand'].col_name:
            if value not in ('+', '-'):
                raise RuntimeError("strand must be '+ or '-' got '{}'".format(value))
        elif 'start' in self._fields_idx and field == self._fields_idx['start'].col_name:
            value = int(value)
        elif 'stop' in self._fields_idx and field == self._fields_idx['stop'].col_name:
            value = int(value)
        return value

    @property
    def ref(self):
        return self._values[self._fields_idx['ref'].idx]

    @property
    def strand(self):
        return self._values[self._fields_idx['strand'].idx]

    @property
    def start(self):
        if 'start' in self._fields_idx:
            return self._values[self._fields_idx['start'].idx]
        else:
            return None
    @property
    def stop(self):
        if 'stop' in self._fields_idx:
            return self._values[self._fields_idx['stop'].idx]
        else:
            return None

    def __str__(self):
        return '\t'.join([str(v) for v in self._values])


Idx = namedtuple('Idx', ('col_name', 'idx'))


def new_entry_type(name, header, ref_col, strand_col='strand', start_col=None, stop_col=None):
    fields = header.split()
    fields_idx = {}
    try:
        fields_idx['ref'] = Idx(ref_col, fields.index(ref_col))
    except ValueError:
        raise RuntimeError("The ref_col {} does not match any fields: {}".format(ref_col, header))
    try:
        fields_idx['strand'] = Idx(strand_col, fields.index(strand_col))
    except ValueError:
        raise RuntimeError("The strand_col {} does not match any fields: {}".format(strand_col, header))
    if start_col:
        try:
            fields_idx['start'] = Idx(start_col, fields.index(start_col))
        except ValueError:
            raise RuntimeError("The start_col {} does not match any fields: {}".format(start_col, header))
        try:
            fields_idx['stop'] = Idx(stop_col, fields.index(stop_col))
        except ValueError:
            raise RuntimeError("The stop_col {} does not match any fields: {}".format(stop_col, header))
    return type(name, (Entry,), {'_fields_idx': fields_idx, '_fields': fields})


class AnnotationParser:

    def __init__(self, path, ref_col, strand_col='strand', start_col=None, stop_col=None):
        self.path = path
        self.ref_col = ref_col
        self.strand_col = strand_col
        self.start_col = start_col
        self.stop_col = stop_col

    def annot_iterator(self):
        """
        parse an annotation file and yield a :class:`Entry` for each line of the file.

        :param path: the path of the annotation file to parse.
        :type path: string
        :return: a generator on a annotation file.
        """
        with open(self.path, 'r') as annot_file:
            header = annot_file.readline().strip()
            MyEntryClass = new_entry_type('MyEntry', header, self.ref_col,
                                          strand_col=self.strand_col,
                                          start_col=self.start_col,
                                          stop_col=self.stop_col)
            for line in annot_file:
                yield MyEntryClass(line.split())

    def max(self):
        if self.start_col is not None:
            max_left = max_right = 0
            for entry in self.annot_iterator():
                left = entry.ref - entry.start
                right = entry.stop - entry.ref
                if left > max_left:
                    max_left = left
                if right > max_right:
                    max_right = max
            return max_left, max_right

if __name__ == '__main__':

    ap = AnnotationParser('../../data/annotations_cerevisiae.txt', 'Position')
    for e in ap.annot_iterator():
        print(e)