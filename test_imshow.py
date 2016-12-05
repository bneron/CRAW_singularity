import os
import matplotlib.pyplot as plt
import pandas as pd

def get_data(coverage_file):
    data = pd.read_table(coverage_file, comment='#', na_values='None')
    return data

def split_data(data):
    sense = data.loc[data['sense'] == 'S']
    antisense = data.loc[data['sense'] == 'AS']
    return sense, antisense

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


def draw_one_matrix(mat, ax, cmap=plt.cm.seismic):
    row_num, col_num = mat.shape
    mat_img = ax.imshow(mat,
                        cmap=cmap,
                        origin='upper',
                        interpolation='nearest',
                        aspect=col_num / row_num,
                        extent=[int(mat.columns[0]), int(mat.columns[-1]), row_num, 0],
                        )

    for ylabel_i in ax.axes.get_yticklabels():
        ylabel_i.set_visible(False)

    for tick in ax.axes.get_yticklines():
        tick.set_visible(False)


def heatmap(sense, antisense, color_map=plt.cm.seismic, title=None):

    draw_sense = True
    draw_antisense = True
    if sense is None or sense.empty:
        draw_sense = False
    if antisense is None or antisense.empty:
        draw_antisense = False

    fig = plt.figure()
    if title is not None:
        fig.suptitle(title)

    if draw_sense and not draw_antisense:
        sense_plot = plt.subplot(111)
        sense_plot.set_title("sense")
        draw_one_matrix(sense, sense_plot, cmap=color_map)
    elif draw_antisense and not draw_sense:
        antisense_plot = plt.subplot(111)
        antisense_plot.set_title("anti sense")
        draw_one_matrix(antisense, antisense_plot, cmap=color_map)
    elif draw_sense and draw_antisense:
        sense_plot = plt.subplot(121)
        sense_plot.set_title("sense")
        draw_one_matrix(sense, sense_plot, cmap=color_map)

        antisense_plot = plt.subplot(122)
        antisense_plot.set_title("anti-sense")
        draw_one_matrix(antisense, antisense_plot, cmap=color_map)
    title = "unknow" if title is None else title.replace(' ', '_')
    #fig.savefig("{}.png".format(title))

cov_file = '/home/bneron/Projects/gwenael/src/crac_tac4_window_200.cov'
cov_file = 'WTE1_0+2000.cov'
d = get_data(cov_file)
sense, antisense = split_data(d)
sense = remove_metadata(sense)
antisense = remove_metadata(antisense)
heatmap(sense, antisense, title=os.path.splitext(os.path.basename(cov_file))[0])
plt.show()
