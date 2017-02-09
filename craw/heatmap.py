from inspect import isfunction
import sys
import logging

import matplotlib.pyplot as plt
import matplotlib as mtp
import pandas as pd

_log = logging.getLogger(__name__)
_log.setLevel(logging.NOTSET)


def get_data(coverage_file):
    """

    :param coverage_file: the path of the coverage file to parse.
    :type coverage_file: str
    :return: the data as 2 dimension dataframe
    :rtype: a :class:`pandas.DataFrame` object
    """
    data = pd.read_table(coverage_file, comment='#', na_values='None')
    return data


def split_data(data):
    """
    Split the matrix in 2 matrices one for sens the other for antisense.

    :param data: the coverage data to split
    :type data: a 2 dimension :class:`pandas.DataFrame` object
    :return: two matrix
    :rtype: tuple of two :class:`pandas.DataFrame` object
    """
    sense = data.loc[data['sense'] == 'S']
    antisense = data.loc[data['sense'] == 'AS']
    return sense, antisense


def sort(data, criteria, **kwargs):
    """
    Sort the matrix in function of criteria.
    This function act as proxy for several specific sorting functions

    :param data: the data to sort.
    :type data: :class:`pandas.DataFrame`.
    :param criteria: which criteria to use to sort the data (by_gene_size, using_col, using_file).
    :type criteria: string.
    :param kwargs: depending of the criteria
     - start_col, stop_col for sort_by_gene_size
     - col for using_col
     - file for using file
    :return: sorted data.
    :rtype: a :class:`pandas.DataFrame` object.
    """
    if data is None or data.empty:
        return data
    func_name = '_sort_' + criteria
    all_func = globals()
    if func_name in all_func and isfunction(all_func[func_name]):
        s_d = globals()[func_name](data, **kwargs)
        return s_d
    else:
        raise RuntimeError('BLABLA')


def _sort_by_gene_size(data, start_col=None, stop_col=None):
    """
    Sort the matrix in function of the gene size.

    :param data: the data to sort.
    :type data: :class:`pandas.DataFrame`.
    :param start_col: the name of the column representing the beginning of the gene.
    :type start_col: string.
    :param stop_col: the name of the column representing the end of the gene.
    :type stop_col: string
    :return: sorted data.
    :rtype: a :class:`pandas.DataFrame` object.
    """
    _log.info("Sorting data by gene size using cols {}:{}".format(start_col, stop_col))
    data['gene_len'] = abs(data[stop_col] - data[start_col])
    data = data.sort_values('gene_len', axis='index')
    del data['gene_len']
    return data


def _sort_using_col(data, col=None):
    """
    Sort the matrix in function of the column col

    :param data: the data to sort.
    :type data: :class:`pandas.DataFrame`.
    :param col: the name of the column to use for sorting the data.
    :type col: string.
    :return: sorted data.
    :rtype: a :class:`pandas.DataFrame` object.
    """
    _log.info("Sorting data using col {}".format(col))
    data = data.sort_values(col, axis='index')
    return data


def _sort_using_file(data, file=None):
    """
    Sort the matrix in function of file.
    The file must have the following structure
    the first line must be the name of the column
    the following lines must be the values, one per line
    each line starting by '#' will be ignore.

    :param data: the data to sort.
    :type data: :class:`pandas.DataFrame`.
    :param file: The file to use as guide to sort the data.
    :type file: a file like object.
    :return: sorted data.
    :rtype: a :class:`pandas.DataFrame` object.
    """
    _log.info("Sorting data using file {}".format(file))
    ref = pd.read_table(file, comment="#")
    col_name = ref.columns[0]

    # change the index of the data using the col_name
    data.set_index(data[col_name], inplace=True)

    # reindex the data according the ref dataframe
    reindexed_data = data.reindex(ref[col_name])
    return reindexed_data


def crop_matrix(data, start_col, stop_col):
    """
    Crop matrix (remove columns). The resulting matrix will be [start_col, stop_col]

    :param data: the data to sort.
    :type data: :class:`pandas.DataFrame`.
    :param start_col: The name of the first column to keep.
    :type start_col: string.
    :param stop_col: The name of the last column to keep.
    :type stop_col: string.
    :return: sorted data.
    :rtype: a :class:`pandas.DataFrame` object.
    """
    if data is None or data.emtpy:
        return data
    return data.loc[:, start_col:stop_col]


def remove_metadata(data):
    """
    Remove all information which is not coverage value (as chromosome, strand, name, ...)

    :param data: the data to sort.
    :type data: :class:`pandas.DataFrame`.
    :return: sorted data.
    :rtype: a :class:`pandas.DataFrame` object.
    """
    if data is None or data.empty:
        return data

    def find_col_2_split(data):
        prev_col = None
        for col in data.columns:
            try:
                int(col)
            except ValueError:
                prev_col = col
                continue
            else:
                break
        return prev_col, col

    last_metadata_col, first_cov_col = find_col_2_split(data)
    coverage_data = data.loc[:, first_cov_col:]
    return coverage_data


def draw_one_matrix(mat, ax, cmap=plt.cm.Blues, y_label=None):
    """
    Draw a matrix sing matplotlib imshow object

    :param mat: the data to represent graphically.
    :type mat: a :class:`pandas.DataFrame` object.
    :param ax: the axis where to represent the data
    :type ax: a :class:`matplotlib.axis` object
    :param cmap: the color map to use to represent the data.
    :type cmap: a :class:`matplotlib.pyplot.cm` object.
    :param y_label: the label for the data draw on y-axis.
    :type y_label: string
    :return: the mtp image corresponding to data
    :rtype: a :class:`matplotlib.image` object.
    """

    row_num, col_num = mat.shape
    mat_img = ax.imshow(mat,
                        cmap=cmap,
                        origin='upper',
                        interpolation='none',
                        aspect=col_num / row_num,
                        extent=[int(mat.columns[0]), int(mat.columns[-2]) + 1, row_num, 0],
                        )
    for ylabel_i in ax.axes.get_yticklabels():
        ylabel_i.set_visible(False)
    for tick in ax.axes.get_yticklines():
        tick.set_visible(False)
    if y_label is not None:
        ax.set_ylabel(y_label, size='large')
    return mat_img


def draw_heatmap(sense, antisense, color_map=plt.cm.Blues, title='', sense_on='top', norm=None, size=None):
    """
    Create a figure with subplot to represent the data as heat map.

    :param sense: the data representing coverage on sense.
    :type sense: a :class:`pandas.DataFrame` object.
    :param antisense: the data representing coverage on anti sense.
    :type sense: a :class:`pandas.DataFrame` object.
    :param color_map: the color map to use to represent the data.
    :type color_map: a :class:`matplotlib.pyplot.cm` object.
    :param title: the figure title (by default the same as the coverage file).
    :type title: string.
    :param sense_on: specify the lay out. Where to place the heat map representing the sense data.
     the available values are: 'left', 'right', 'top', 'bottom' (default = 'top').
    :type sense_on: string.
    :param norm: a color normalisation.
    :type norm: :class:`mtp.colors.Normalize` object.
    :param size: the size of the figure in inches (wide, height).
    :type size: tuple of 2 float.
    :return: The figure.
    :rtype: a :class:`matplotlib.pyplot.Figure` object.
    """

    draw_sense = True
    draw_antisense = True
    if sense is None or sense.empty:
        draw_sense = False
    if antisense is None or antisense.empty:
        draw_antisense = False

    if all((draw_sense, draw_antisense)):
        if sense_on in ('bottom', 'top'):
            if size is None:
                size = (7, 10)
            fig, axes_array = plt.subplots(nrows=2, ncols=1, figsize=size)
            if sense_on == 'top':
                sense_subplot, antisense_subplot = axes_array
            else:
                antisense_subplot, sense_subplot = axes_array
        elif sense_on in ('left', 'right'):
            if size is None:
                size = (10, 7)
            fig, axes_array = plt.subplots(nrows=1, ncols=2, figsize=size)
            if sense_on == 'left':
                sense_subplot, antisense_subplot = axes_array
            else:
                antisense_subplot, sense_subplot = axes_array

    elif any((draw_sense, draw_antisense)):
        if size is None:
            size = (6, 6)
        fig, axes_array = plt.subplots(nrows=1, ncols=1, figsize=size)
        if draw_sense:
            sense_subplot = axes_array
        else:
            antisense_subplot = axes_array
    else:
        _log.warning("No matrix to draw")
        return

    fig.suptitle(title, fontsize='large')


    if draw_sense:
        _log.info("Drawing sense matrix")
        sense_img = draw_one_matrix(sense, sense_subplot, cmap=color_map, y_label="Sense")
        if norm:
            sense_img.set_norm(norm)
    if draw_antisense:
        _log.info("Drawing antisense matrix")
        antisense_img = draw_one_matrix(antisense, antisense_subplot, cmap=color_map, y_label="Anti sense")
        if norm:
            antisense_img.set_norm(norm)

    fig.suptitle(title,
                 fontsize=12,
                 horizontalalignment='center',
                 verticalalignment='top')

    fig.tight_layout()
    if sense_on in ('top', 'bottom'):
        fig.subplots_adjust(top=0.95)
    fig.canvas.set_window_title(title)
    return fig


def normalize(data):
    pass

def log_norm(data):
    pass

def normalize_row_by_row(data):
    pass

def draw_raw_image(data, color_map=lt.cm.Blues, format='png'):
    pass