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

def crop(data, start_col, stop_col):
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

if __name__ == '__main__':
    import sys
    cov_file = sys.argv[1]
    #cov_file = 'WTE1_0+2000.new.cov'
    #cov_file = 'WTE1_var_window.cov'

    d = get_data(cov_file)

    sense, antisense = split_data(d)
    sense = sort(sense, 'annotation_start', 'annotation_end')
    sense = remove_metadata(sense)

    antisense = sort(antisense, 'annotation_start', 'annotation_end')
    antisense = remove_metadata(antisense)

    sense = crop(sense, '0', '2000')
    antisense = crop(antisense, '0', '2000')


    import math
    sense_log = sense.applymap(lambda x: math.log(x, 10) if x > 0 else 0)
    antisense_log = antisense.applymap(lambda x: math.log(x, 10) if x > 0 else 0)

    #color_map = plt.cm.Blues
    #color_map = plt.cm.seismic
    #map_2_test = [plt.cm.Blues, plt.cm.YlOrRd, plt.cm.cool, plt.cm.BrBG, plt.cm.seismic]
    map_2_test = [plt.cm.Blues]
    for color_map in map_2_test:
        ##################
        # sense
        ##################

        #fig = plt.figure("data normalization")
        cols = ('max = 1000', 'logNorm', 'normalized data')
        rows = ('Sense', 'AntiSense')

        fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(14, 10))
        fig.suptitle(color_map.name, fontsize='large')
        #print("axes = ", axes)

        #max_plot = plt.subplot(2, 3, 1)
        #max_plot.set_title("max = 1000")
        max_plot = axes[0][0]
        s_im = draw_one_matrix(sense, max_plot, cmap=color_map)
        norm_1000 = mtp.colors.Normalize(vmin=0, vmax=1000)
        s_im.set_norm(norm_1000)

        #max_plot = plt.subplot(2, 3, 2)
        #max_plot.set_title("logNorm")
        max_plot = axes[0][1]
        s_im = draw_one_matrix(sense, max_plot, cmap=color_map)
        log_norm = mtp.colors.LogNorm()
        s_im.set_norm(log_norm)

        #max_plot = fig.add_subplot(2, 3, 3)
        #max_plot.set_title("normalized data (log 10)")
        max_plot = axes[0][2]
        s_im = draw_one_matrix(sense_log, max_plot, cmap=color_map)

        #####################
        # antisense
        #####################

        #max_plot = plt.subplot(2, 3, 4)
        #max_plot.set_title("max = 1000")
        max_plot = axes[1][0]
        s_im = draw_one_matrix(antisense, max_plot, cmap=color_map)
        norm_1000 = mtp.colors.Normalize(vmin=0, vmax=1000)
        s_im.set_norm(norm_1000)

        #max_plot = plt.subplot(2, 3, 5)
        #max_plot.set_title("logNorm")
        max_plot = axes[1][1]
        s_im = draw_one_matrix(antisense, max_plot, cmap=color_map)
        log_norm = mtp.colors.LogNorm()
        s_im.set_norm(log_norm)

        #max_plot = fig.add_subplot(2, 3, 6)
        #max_plot.set_title("normalized data")
        max_plot = axes[1][2]
        s_im = draw_one_matrix(antisense_log, max_plot, cmap=color_map)

        for ax, col in zip(axes[0], cols):
            ax.set_title(col)
        for ax, row in zip(axes[:, 0], rows):
            ax.set_ylabel(row, size='large')

        fig.tight_layout()

    plt.show()
