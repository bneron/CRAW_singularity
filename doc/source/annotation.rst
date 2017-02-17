.. _annotation:

==========
annotation
==========

The annotation module contains everything that is needed to parse annotation file and handle it.

AnnotationParser
----------------

The entry point to parse an annotation file is the :class:`craw.annotation.AnnotationParser`.
An annotation parser have two methods:

 * :meth:`craw.annotation.AnnotationParser.get_annotations` create a new type of Entry and iterate over the annotation file and for each line return a new
   instance of the newly :class:`craw.annotation.Entry` class it just create on the fly.
 * the other more technique give the maximum of nucleotides before and after the reference.
   It is needed to compute the size of the resulting matrix.

The force of this approach is to generate a new type of entry for each parsing. So it's very flexible and allow to fit with
most of annotation file. But for one file, all the parsing use the same *Entry* class so it ensure the coherence in data.


new_entry_type
--------------

Is a factory which generate a new subclass of :class:`craw.annotation.Entry` given the fields gather form the annotation file header
(first line non starting with #) and the columns semantic given by the user.
The first role of this factory is to check if all parameter given by user correspond ot header and do some coherence checking.
If everything seems Ok it generate on the fly a new subclass of :class:`craw.annotation.Entry`.

Entry Class
-----------

An Entry correspond to one line of the annotation file.

The Entry convert values if necessary (*strand* in a internal representation +/-, *position* in integer ...).
It also expose a generic api to access some fields whatever the named of the columns.


annotation API reference
========================

.. automodule:: craw.annotation
   :members:
   :private-members:
   :special-members:



