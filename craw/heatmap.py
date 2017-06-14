###########################################################################
#                                                                         #
# This file is part of Counter RNAseq Window (craw) package.              #
#                                                                         #
#    Authors: Bertrand Néron                                              #
#    Copyright © 2017  Institut Pasteur (Paris).                          #
#    see COPYRIGHT file for details.                                      #
#                                                                         #
#    craw is free software: you can redistribute it and/or modify         #
#    it under the terms of the GNU General Public License as published by #
#    the Free Software Foundation, either version 3 of the License, or    #
#    (at your option) any later version.                                  #
#                                                                         #
#    craw is distributed in the hope that it will be useful,              #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of       #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                 #
#    See the GNU General Public License for more details.                 #
#                                                                         #
#    You should have received a copy of the GNU General Public License    #
#    along with craw (see COPYING file).                                  #
#    If not, see <http://www.gnu.org/licenses/>.                          #
#                                                                         #
###########################################################################

from inspect import isfunction
import logging

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageColor

_log = logging.getLogger(__name__)


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
    Split the matrix in 2 matrices one for sense the other for antisense.

    :param data: the coverage data to split
    :type data: a 2 dimension :class:`pandas.DataFrame` object
    :return: two matrix
    :rtype: tuple of two :class:`pandas.DataFrame` object (sense pandas.DataFrame, antisense pandas.DataFrame)
    """
    grp = data.groupby(by=['sense'])
    return grp.get_group('S').copy(), grp.get_group('AS').copy()


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
        raise RuntimeError("The '{}' sorting method does not exists.".format(criteria))


def _sort_by_gene_size(data, start_col=None, stop_col=None, ascending=True):
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
    sorted_data = data.sort_values('gene_len', axis='index', ascending=ascending)
    del data['gene_len']
    del sorted_data['gene_len']
    return sorted_data


def _sort_using_col(data, col=None, ascending=True):
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
    if col is None:
        raise RuntimeError("You must specify the column used to sort.")
    data = data.sort_values(col, axis='index', ascending=ascending)
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
    :type data: a 2D :class:`pandas.DataFrame` object.
    :param start_col: The name of the first column to keep.
    :type start_col: string.
    :param stop_col: The name of the last column to keep.
    :type stop_col: string.
    :return: sorted data.
    :rtype: a 2D :class:`pandas.DataFrame` object or None if data is None.
    """
    if data is None or data.empty:
        return data
    return data.loc[:, start_col:stop_col]


def lin_norm(data):
    """
    Normalize data with linear algorithm.
    The formula applied to obtain the results is:

        zi = xi - min(x) / max(x) - min(x)

    where x=(x1,...,xn) and zi is now your with normalized data.
    Ensure that the resulting values are comprise between 0 and 1.
    return None if data is None, return empty :class:`pd.DataFrame` object if data is empty.

    :param data: the data to normalize, this 2D matrix must contains only coverage scores (no more metadata).
    :type data: a 2D :class:`pandas.DataFrame` object.
    :return: a normalize matrix, where  0 <= zi <=1 where z=(z1, ..., zn)
    :rtype: a 2D :class:`pandas.DataFrame` object or None if data is None.
    """
    if data is None or data.empty:
        return data
    # data is a 2D DataFrame
    # so min() return a Series of min by columns
    # and min.min is the min of the matrix
    # idem for max
    vmin = data.min().min()
    vmax = data.max().max()
    data -= vmin
    data /= vmax - vmin
    return data


def log_norm(data):
    """
    The base 10 logarithm is compute for all values before
    a normalization (see :func:`normalize` ) to ensure
    that all values are comprise between 0 and 1 .

    .. note::
        coverage scores are integers >= 0.
        log10(0) = -inf or warning in macos
        prior to normalize data the 0 values are replace by 1.

    :param data: the data to normalize, this 2D matrix must contains only coverage scores (no more metadata).
    :type data: a 2D :class:`pandas.DataFrame` object.
    :return: a normalize matrix, where  0 <= zi <=1 where z=(z1, ..., zn)
    :rtype: a 2D :class:`pandas.DataFrame` object or None if data is None.
    """
    if data is None or data.empty:
        return data
    # log10(0) = -inf or warning on macos
    # so transform all value below 1 to 1
    # normally only 0 because coverage scores are integers >=0
    data = np.maximum(data, 1)
    data = np.log10(data)
    data = lin_norm(data)
    return data


def lin_norm_row_by_row(data):
    """
    Normalize data with linear algorithm but
    instead to normalize all the matrix, the normalization formula
    (see :func:`normalize`) is applied row by row.
    It ensure that all values are between 0 and 1.

    :param data: the data to normalize, this 2D matrix must contains only coverage scores (no more metadata).
    :type data: a 2D :class:`pandas.DataFrame` object.
    :return: a normalize matrix, where  0 <= zi <=1 where z=(z1, ..., zn)
    :rtype: a 2D :class:`pandas.DataFrame` object or None if data is None.
    """
    if data is None or data.empty:
        return data
    # axis=0 operation on columns
    # axis=1 operation on rows

    # subtract for each value the min of the corresponding row
    norm_data = data.sub(data.min(axis=1), axis=0)
    # compute a series for each row value = max - min
    amplitude = data.max(axis=1) - data.min(axis=1)
    # for each value divide this value by the corresponding amplitude of the row
    norm_data = norm_data.div(amplitude, axis=0)
    return norm_data


def log_norm_row_by_row(data):
    """
    as :func:`normalize_row_by_row` but prior normalisation
    a 10 base logarithm is applied.

    .. note::
        coverage scores are integers >= 0.
        log10(0) = -inf
        to normalize data the -inf value are change in 0.

    :param data: the data to normalize, this 2D matrix must contains only coverage scores (no more metadata).
    :type data: a 2D :class:`pandas.DataFrame` object.
    :return: a normalize matrix, where  0 <= zi <=1 where z=(z1, ..., zn)
    :rtype: a 2D :class:`pandas.DataFrame` object or None if data is None.
    """
    if data is None or data.empty:
        return data
    # log10(0) = -inf or warning on macos
    # so transform all negative values below 1
    # normally only 0 because coverage scores are integers >=0
    data = np.maximum(data, 1)
    data = np.log10(data)
    data = lin_norm_row_by_row(data)
    return data


def remove_metadata(data):
    """
    Remove all information which is not coverage value (as chromosome, strand, name, ...)

    :param data: the data coming from a coverage file parsing containing coverage information and metadata chromosome,
                 gene name , ...
    :type data: :class:`pandas.DataFrame`.
    :return: sorted data.
    :rtype: a 2D :class:`pandas.DataFrame` object or None if data is None.
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


def draw_one_matrix(mat, ax, cmap=plt.cm.Blues, y_label=None, marks=None):
    """
    Draw a matrix using matplotlib imshow object

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
    # Classes that are ‘array-like’ such as pandas data objects and np.matrix may or may not work as intended.
    # It is best to convert these to np.array objects prior to plotting.
    # http://matplotlib.org/faq/usage_faq.html#types-of-inputs-to-plotting-functions
    mat_img = ax.imshow(mat.values,
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
    if marks:
        for mark in marks:
            ax.axvline(x=mark.to_px(), linewidth=0.5, color=mark.rgb_float())
    return mat_img


def draw_heatmap(sense, antisense, color_map=plt.cm.Blues, title='', sense_on='top', size=None, marks=None):
    """
    Create a figure with subplot to represent the data as heat map.

    :param sense: the data normalized (xi in [0,1]) representing coverage on sense.
    :type sense: a :class:`pandas.DataFrame` object.
    :param antisense: the data normalized (xi in [0,1]) representing coverage on anti sense.
    :type sense: a :class:`pandas.DataFrame` object.
    :param color_map: the color map to use to represent the data.
    :type color_map: a :class:`matplotlib.pyplot.cm` object.
    :param title: the figure title (by default the same as the coverage file).
    :type title: string.
    :param sense_on: specify the lay out. Where to place the heat map representing the sense data.
     the available values are: 'left', 'right', 'top', 'bottom' (default = 'top').
    :type sense_on: string.
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
        sense_img = draw_one_matrix(sense, sense_subplot, cmap=color_map, y_label="Sense", marks=marks)
    if draw_antisense:
        _log.info("Drawing antisense matrix")
        antisense_img = draw_one_matrix(antisense, antisense_subplot, cmap=color_map, y_label="Anti sense", marks=marks)

    fig.suptitle(title,
                 fontsize=12,
                 horizontalalignment='center',
                 verticalalignment='top')

    fig.tight_layout()
    if sense_on in ('top', 'bottom'):
        fig.subplots_adjust(top=0.95)
    fig.canvas.set_window_title(title)
    return fig


def draw_raw_image(data, out_name, color_map=plt.cm.Blues, format='PNG', marks=None):
    """
    Generate an image file with one pixel for each values of the data matrix.
    the data can be either the coverage on sense or on antisense.

    :param data: a **Normalized** (where all values are between 0 and 1) matrix.
    :type data: 2D :class:`pandas.DataFrame` or :class:`numpy.array` object
    :param out_name: The name of the generated graphic file.
    :type out_name: string
    :param color_map:
    :type color_map:
    :param format: the format of the result png, jpeg, ... (see pillow supported formats)
    :type format: string
    :raise: RuntimeError if data are not normalized.
    """
    # the data was normalized to ensure that values are in [0,1]
    # I apply directly a matplotlib.colormap on this values
    # then rescale the values to [0, 255]
    # and convert them in integer 8 bits
    # save it with pillow
    # see the whole recipe at
    # http://stackoverflow.com/questions/10965417/how-to-convert-numpy-array-to-pil-image-applying-matplotlib-colormap
    if data.max().max() > 1 or data.min().min() < 0:
        raise RuntimeError("data must be normalized (between [0,1])")
    im = Image.fromarray(np.uint8(color_map(data)*255))
    if marks:
        width, height = im.size
        draw = ImageDraw.Draw(im)
        for mark in marks:
            draw.line([(mark.to_px(), 0), (mark.to_px(), height)], fill=mark.rgb_int)
    im.save(out_name, format=format)


class Mark:
    """A mark is a position and a color tight together. 
    It is used to draw a colored vertical line at the given position on the heatmap"""

    def __init__(self, pos, data, color_map, color=None):
        """
        
        :param pos: The position where to draw a mark, the position is relative to the reference position (0)
        :type pos: int
        :param data: the coverage matrix 
        :type data: :class:`pandas.DataFrame` object
        :param color_map: the color map used to draw the heatmap
        :type color_map: class`:matplotlib.pyplot.ColorMap` object
        :param color: the color of the line, the supported formats are
                        - hexadecimal values as #rgb or #rrggbb, for instance #ff0000 is pure red.
                        - common html color names
                        
        """
        self._min, self._max = self._get_matrix_bound(data)
        if self._min <= pos <= self._max:
            self.pos = pos
        else:
            raise ValueError("mark position must be {} >= pos >= {}: provide {}".format(self._min, self._max, pos))
        self.color_map = color_map
        self._color = self._color_converter(color, data)

    def rgb_int(self):
        return self._color

    def rgb_float(self):
        return tuple([c / 255 for c in self._color])


    def _color_converter(self, color, data):
        """
        
        :param color: the color of the line, the supported formats are
                        - hexadecimal values as #rgb or #rrggbb, for instance #ff0000 is pure red.
                        - common html color names
        :type color: string
        :param data: the matrix coverage
        :type data: :class:`pandas.DataFrame` object
        :return: rgb color
        :rtype: tuple with 3 int between 0 and 255  
        """
        if color:
            try:
                color = ImageColor.getrgb(color)
            except ValueError:
                raise ValueError("{} is not a valid color".format(color))
        else:
            color = tuple([int(round(c * 255)) for c in self.color_map(1.0)][:-1])
        return color


    def _get_matrix_bound(self, data):
        """
        :param data: the matrix coverage
        :type data: :class:`pandas.DataFrame` object 
        :return: the most right and left position of the coverage 
        :rtype: tuple of 2 int
        """
        for col in data.columns:
            try:
                int(col)
            except ValueError:
                continue
            else:
                break
        return int(col), int(data.columns[-1])


    def to_px(self):
        """
        tanslate the position of the mark relative to the reference in pixel.
        :return: the position of the mark in pixel.
        :rtype: positive int
        """
        return self.pos - self._min