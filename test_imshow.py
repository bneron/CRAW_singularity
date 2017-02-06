import os
import matplotlib.pyplot as plt
import matplotlib as mtp
import pandas as pd


def get_data(coverage_file):
    data = pd.read_table(coverage_file, comment='#', na_values='None')
    return data


def split_data(data):
    sense = data.loc[data['sense'] == 'S']
    antisense = data.loc[data['sense'] == 'AS']
    return sense, antisense


def sort(data, start_col, stop_col):
    data['gene_len'] = abs(data[stop_col] - data[start_col])
    data = data.sort_values('gene_len', axis='index')
    del data['gene_len']
    return data


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


def draw_one_matrix(mat, ax, cmap=plt.cm.seismic, vmin=None, vmax=None):
    row_num, col_num = mat.shape
    mat_img = ax.imshow(mat,
                        cmap=cmap,
                        origin='upper',
                        interpolation='nearest',
                        aspect=col_num / row_num,
                        extent=[int(mat.columns[0]), int(mat.columns[-2])+ 1, row_num, 0],
                        )
    for ylabel_i in ax.axes.get_yticklabels():
        ylabel_i.set_visible(False)
    for tick in ax.axes.get_yticklines():
        tick.set_visible(False)
    return mat_img


def heatmap(sense, antisense, color_map=plt.cm.seismic, title=None, sense_title=None, vmin=None, vmax=None):

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
        sense_plot.set_title(sense_title if sense_title is not None else "sense")
        sense_im = draw_one_matrix(sense, sense_plot, cmap=color_map, vmin=vmin, vmax=vmax)
        fig.colorbar(sense_im)
    elif draw_antisense and not draw_sense:
        antisense_plot = plt.subplot(111)
        antisense_plot.set_title("anti sense")
        antisense_im = draw_one_matrix(antisense, antisense_plot, cmap=color_map, vmin=vmin, vmax=vmax)
        fig.colorbar(antisense_im)
    elif draw_sense and draw_antisense:
        sense_plot = plt.subplot(121)
        sense_plot.set_title(sense_title if sense_title is not None else "sense")
        sense_im = draw_one_matrix(sense, sense_plot, cmap=color_map, vmin=vmin, vmax=vmax)
        fig.colorbar(sense_im)

        antisense_plot = plt.subplot(122)
        antisense_plot.set_title("anti-sense")
        antisense_im = draw_one_matrix(antisense, antisense_plot, cmap=color_map, vmin=vmin, vmax=vmax)
        fig.colorbar(antisense_im)
    title = "unknow" if title is None else title.replace(' ', '_')
    w = 16
    h = 12
    #fig.set_size_inches(16,12)
    #fig.savefig("{}_{}x{}.png".format(title, w * 2.54, h * 2.54))
    return fig, sense_im


cov_file = '/home/bneron/Projects/gwenael/src/crac_tac4_window_200.cov'
cov_file = 'WTE1_0+2000.test'
d = get_data(cov_file)

sense, antisense = split_data(d)
sense = sort(sense, 'annotation_start', 'annotation_end')
sense = remove_metadata(sense)

antisense = sort(antisense, 'annotation_start', 'annotation_end')
antisense = remove_metadata(antisense)


#sense['total'] = sense.sum(axis=1)


color_map = plt.cm.Blues
#color_map = plt.cm.seismic

##################
# sense
##################

#fig = plt.figure("data normalization")
cols = ('max = 1000', 'logNorm', 'normalized data')
rows = ('Sense', 'AntiSense')

fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(12, 8))

print("axes = ", axes)

for ax, col in zip(axes[0], cols):
    ax.set_title(col)
for ax, row in zip(axes[:, 0], rows):
    ax.set_ylabel(row, rotation=0, size='large')
# max_plot = plt.subplot(2, 3, 1)
# max_plot.set_title("max = 100")
# s_im = draw_one_matrix(sense, max_plot, cmap=color_map)
# norm_100 = mtp.colors.Normalize(vmin=0, vmax=100)
# s_im.set_norm(norm_100)

#max_plot = plt.subplot(2, 3, 1)
#max_plot.set_title("max = 1000")
max_plot = axes[0][0]
s_im = draw_one_matrix(sense, max_plot, cmap=color_map)
norm_1000 = mtp.colors.Normalize(vmin=0, vmax=1000)
s_im.set_norm(norm_1000)

# max_plot = plt.subplot(2, 3, 3)
# max_plot.set_title("max = 10000")
# s_im = draw_one_matrix(sense, max_plot, cmap=color_map)
# norm_10000 = mtp.colors.Normalize(vmin=0, vmax=10000)
# s_im.set_norm(norm_10000)

#max_plot = plt.subplot(2, 3, 2)
#max_plot.set_title("logNorm")
max_plot = axes[0][1]
s_im = draw_one_matrix(sense, max_plot, cmap=color_map)
log_norm = mtp.colors.LogNorm()
s_im.set_norm(log_norm)

import math
sense_log = sense.applymap(lambda x: math.log(x, 10) if x > 0 else 0)
#max_plot = fig.add_subplot(2, 3, 3)
#max_plot.set_title("normalized data (log 10)")
max_plot = axes[0][2]
s_im = draw_one_matrix(sense_log, max_plot, cmap=color_map)


#####################
# antisense
#####################

#fig = plt.figure("data normalization antisense")


# max_plot = plt.subplot(2, 3, 1)
# max_plot.set_title("max = 100")
# s_im = draw_one_matrix(antisense, max_plot, cmap=color_map)
# norm_100 = mtp.colors.Normalize(vmin=0, vmax=100)
# s_im.set_norm(norm_100)

#max_plot = plt.subplot(2, 3, 4)
#max_plot.set_title("max = 1000")
max_plot = axes[1][0]
s_im = draw_one_matrix(antisense, max_plot, cmap=color_map)
norm_1000 = mtp.colors.Normalize(vmin=0, vmax=1000)
s_im.set_norm(norm_1000)

# max_plot = plt.subplot(2, 3, 5)
# max_plot.set_title("max = 10000")
# s_im = draw_one_matrix(antisense, max_plot, cmap=color_map)
# norm_10000 = mtp.colors.Normalize(vmin=0, vmax=10000)
# s_im.set_norm(norm_10000)

#max_plot = plt.subplot(2, 3, 5)
#max_plot.set_title("logNorm")
max_plot = axes[1][1]
s_im = draw_one_matrix(antisense, max_plot, cmap=color_map)
log_norm = mtp.colors.LogNorm()
s_im.set_norm(log_norm)

import math
antisense_log = antisense.applymap(lambda x: math.log(x, 10) if x > 0 else 0)
#max_plot = fig.add_subplot(2, 3, 6)
#max_plot.set_title("normalized data")
max_plot = axes[1][2]
s_im = draw_one_matrix(antisense_log, max_plot, cmap=color_map)

print("axes = ", axes)

fig.tight_layout()
plt.show()
