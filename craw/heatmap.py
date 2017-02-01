from inspect import isfunction
import sys

import matplotlib.pyplot as plt
import matplotlib as mtp
import pandas as pd


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
    func_name = '_sort_' + criteria
    all_func = globals()
    if func_name in all_func and isfunction(all_func[func_name]):
        s_d = globals()[func_name](data, **kwargs)
        return s_d
    else:
        raise RuntimeError('BLABLA')

def _sort_by_gene_size(data, start_col=None, stop_col=None):
    data['gene_len'] = abs(data[stop_col] - data[start_col])
    data = data.sort_values('gene_len', axis='index')
    del data['gene_len']
    return data


def _sort_using_col(data, col=None):
    data = data.sort_values(col, axis='index')
    return data


def _sort_using_file(data, file=None):
    ref = pd.read_table(file, comment="#")
    col_name = ref.columns[0]

    # change the index of the data using the col_name
    data.set_index(data[col_name], inplace=True)

    # reindex the data according the ref dataframe
    reindexed_data = data.reindex(ref[col_name])
    return reindexed_data


def crop_matrix(data, start_col, stop_col):
    return data.loc[:, start_col:stop_col]


def remove_metadata(data):
    """

    :param data:
    :return:
    """
    def find_col_2_split(data):
        prev_col = None
        for col in data.columns:
            try:
                int(col)
            except:
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


def draw_heatmap(sense, antisense, color_map=plt.cm.Blues):
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
        print("there is n")
    if antisense is None or antisense.empty:
        draw_antisense = False

    if all((draw_sense, draw_antisense)):
        print("all")
        fig, axes_array = plt.subplots(nrows=2, ncols=1, figsize=(7, 10))
    elif any((draw_sense, draw_antisense)):
        print("any")
        fig, axes_array = plt.subplots(nrows=1, ncols=1, figsize=(7, 10))
    else:
        print("WARNING: no matrix to draw", file=sys.stderr)
        return
    fig.suptitle("figure name", fontsize='large')

    log_norm = mtp.colors.LogNorm()

    if draw_sense:
        print("INFO: drawing sense matrix", file=sys.stderr)
        sense_subplot = axes_array[0]
        sense_img = draw_one_matrix(sense, sense_subplot, cmap=color_map, y_label="Sense")
        sense_img.set_norm(log_norm)
    if draw_antisense:
        print("INFO: drawing sense matrix", file=sys.stderr)
        if draw_sense:
            antisense_subplot = axes_array[1]
        else:
            antisense_subplot = axes_array[0]
        antisense_img = draw_one_matrix(antisense, antisense_subplot, cmap=color_map, y_label="Anti sense")
        antisense_img.set_norm(log_norm)

    fig.suptitle('this is the figure title',
                 fontsize=12,
                 horizontalalignment='center',
                 verticalalignment='bottom')

    fig.tight_layout()
    fig.subplots_adjust(top=0.95)
    fig.canvas.set_window_title('window name')
    return fig
