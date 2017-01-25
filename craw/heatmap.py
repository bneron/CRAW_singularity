import matplotlib.pyplot as plt
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


def sort(data, start_col, stop_col):
    data['gene_len'] = abs(data[stop_col] - data[start_col])
    data = data.sort_values('gene_len', axis='index')
    del data['gene_len']
    return data


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


def heatmap(sense, antisense, color_map=plt.cm.seismic):
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

    fig = plt.figure()
    if draw_sense:
        plt.subplot(111)
        draw_one_matrix(sense, y_label='sense')
    if draw_antisense:
        if draw_sense:
            plt.subplot(121)
        else:
            plt.subplot(111)
        draw_one_matrix(antisense, ylabel='anti sense')

    fig.suptitle('this is the figure title',
                 fontsize=12,
                 horizontalalignment='center',
                 verticalalignment='bottom')

    plt.show()

