
class Entry:

    def __init__(self, values):
        values = [self._convert(f, v) for f, v in zip(self._fields, values)]
        self._values = values

    def _convert(self, field, value):
        if field == self._idx['ref'][0]:
            value = int(value)
        elif field == self._idx['strand'][0]:
            if value not in ('+', '-'):
                raise RuntimeError("strand must be '+ or '-' got {}".format(value))
        elif field == self._idx['start'][0]:
            value = int(value)
        elif field == self._idx['stop'][0]:
            value = int(value)
        return value

    def __str__(self):
        attrs = [(attr, val) for attr, val in self.__dict__.items()]
        attrs.sort(key=lambda x: x[0])
        l = []
        for attr, value in attrs:
            l.append('{}={}'.format(attr, value))
        return '{}({})'.format(self.__class__.__name__, ', '.join(l))


def generate_entry(header, ref_col, strand_col='strand', start_col=None, stop_col=None):
    entry_class = type(Entry)
    entry_class._fields = header.split()
    entry_class._idx = {}
    entry_class._idx['ref'] = (ref_col, entry_class._fields.index(ref_col))
    entry_class._idx['strand'] = (strand_col, entry_class._fields.index(strand_col))
    if start_col:
        entry_class._idx['start'] = (start_col, entry_class._fields.index(start_col))
        entry_class._idx['stop'] = (stop_col, entry_class._fields.index(stop_col))

    return entry_class



def annot_iter(path):
    """
    parse an annotation file and yield a :class:`Entry` for each line of the file.

    :param path: the path of the annotation file to parse.
    :type path: string
    :return: a generator on a annotation file.
    """
    with open(path, 'r') as annot_file:
        header = annot_file.readline().strip()
        fields = header.split()
        line = annot_file.readline()
        while line:
            values = line.split()
            new_entry = Entry(**{k: v for k, v in zip(fields, values)})
            yield new_entry
            line = annot_file.readline()



def max(annot_file_path, ref, start, stop):
    with open(annot_file_path, 'r') as annot_file:
        headers = annot_file.readline().split()
        start_idx = headers.index(start)
        stop_idx = headers.index(stop)
        ref_idx = headers.index(ref)
        max_left = max_right = 0
        with annot_file as line:
            values = line.split()
            ref = int(values[ref_idx])
            start = int(values[start_idx])
            stop = int(values[stop_idx])
            left = ref - start
            right = stop - ref
            if left > max_left:
                max_left = left
            if right > max_right:
                max_right = max
    return (max_left, max_right)