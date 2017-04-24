.. heatmap:


=======
heatmap
=======

:func:`craw.heatmap.split_data` and  :func:`craw.heatmap.sort` work on data freshly parsed from coverage file.
That mean that the data contain the metadata (all columns which are not coverage scores
like chromosome, position strand , on so on)

The other functions sort normalization function work on pandas 2D DataFrame or numpy arrays containing only
scores of coverage. That means all metadata was removed (:func:`craw.heatmap.remove_metadata`).


sort
----
The is one public sort function which act as proxy for several private
sorting function.

normalisation
-------------

Several functions to normalize data.

The data can be normalize using min max of the whole data.
Or the min max is recalculated for each row.

in both case the formula is

zi = xi - min(x) / max(x) - min(x)

where x=(x1,...,xn) and zi is now your with normalized data.
in first case x is the whole matrix in 2nd is the row.

Normalization can be precede by 10 base log transformation.

.. note::
    In this case all 0 values are replace by 1 (10 base log is not define)


drawing heatmap
---------------

There are 2 way to generates figures, the first one is to generate a figures containing 2 heatmap for sense or antisense
with axis, legend on so on. But in this representation it's not possible to display a figure with no scaling out/in.
So the information of one pixel is not accessible. This representation is generate by :func:`craw.heatmap.draw_heatmap`
and use matplotlib.

The second representation is to produce raw image where one nucleotide (one position for one gene) is represent by one pixel
without any scale in/out.
In this representation there si not axis legend on so on it's only a raw image.

heatmap API reference
=====================

  .. automodule:: craw.heatmap
    :members:
    :private-members:
    :special-members:

