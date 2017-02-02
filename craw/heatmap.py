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

    :param coverage_file:
    :return:
    """
    data = pd.read_table(coverage_file, comment='#', na_values='None')
    return data


def split_data(data):
    """

    :param data:
    :return:
    """
    sense = data.loc[data['sense'] == 'S']
    antisense = data.loc[data['sense'] == 'AS']
    return sense, antisense


def sort(data, criteria, **kwargs):
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
    _log.info("Sorting data by gene size using cols {}:{}".format(start_col, stop_col))
    data['gene_len'] = abs(data[stop_col] - data[start_col])
    data = data.sort_values('gene_len', axis='index')
    del data['gene_len']
    return data


def _sort_using_col(data, col=None):
    _log.info("Sorting data using col {}".format(col))
    data = data.sort_values(col, axis='index')
    return data


def _sort_using_file(data, file=None):
    _log.info("Sorting data using file {}".format(file))
    ref = pd.read_table(file, comment="#")
    col_name = ref.columns[0]

    # change the index of the data using the col_name
    data.set_index(data[col_name], inplace=True)

    # reindex the data according the ref dataframe
    reindexed_data = data.reindex(ref[col_name])
    return reindexed_data


def crop_matrix(data, start_col, stop_col):
    if data is None or data.emtpy:
        return data
    return data.loc[:, start_col:stop_col]


def remove_metadata(data):
    """

    :param data:
    :return:
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


def draw_one_matrix(mat, ax, cmap=plt.cm.seismic, y_label=None):
    """

    :param mat:
    :param cmap:
    :return:
    """
    row_num, col_num = mat.shape
    mat_img = ax.imshow(mat,
                        cmap=cmap,
                        origin='upper',
                        interpolation='nearest',
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


def draw_heatmap(sense, antisense, color_map=plt.cm.Blues, title='', sense_on='top', norm=None):
    """

    :param sense:
    :param antisense:
    :param color_map:
    :return:
    """
    draw_sense = True
    draw_antisense = True
    if sense is None or sense.empty:
        draw_sense = False
    if antisense is None or antisense.empty:
        draw_antisense = False

    if all((draw_sense, draw_antisense)):
        if sense_on in ('bottom', 'top'):
            fig, axes_array = plt.subplots(nrows=2, ncols=1, figsize=(7, 10))
            if sense_on == 'top':
                sense_subplot, antisense_subplot = axes_array
            else:
                antisense_subplot, sense_subplot = axes_array
        elif sense_on in ('left', 'right'):
            fig, axes_array = plt.subplots(nrows=1, ncols=2, figsize=(10, 7))
            if sense_on == 'left':
                sense_subplot, antisense_subplot = axes_array
            else:
                antisense_subplot, sense_subplot = axes_array

    elif any((draw_sense, draw_antisense)):
        fig, axes_array = plt.subplots(nrows=1, ncols=1, figsize=(7, 7))
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
