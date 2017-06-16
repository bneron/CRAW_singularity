.. _installation:

============
Installation
============


Requirements
============

For craw_coverage
-----------------

  - python > 3
  - pysam >= 0.9.1.4

For craw_htmp
-------------

  - python > 3
  - pysam >= 0.9.1.4
  - pandas >= 0.17.1
  - numpy >= 1.11.2
  - matplotlib >= 1.5.3
  - pillow >= 3.4.2


Installation
============

Installation from package
-------------------------

Using pip ::

    pip install craw

Do not forget to configure the `matplotlib` backend, specially if you use virtualenv.
Otherwise on some platform there won't any output.
See :ref:`matplotlibrc` for more explanation.

.. note::
    On MacOS install python > 3 from image on http://python.org . Then
    install craw using pip ::

        pip3 install craw

    craw will be installed in `/Library/Framework/Python.Framework/Version/3.6/`
    So if you want to use directly craw_coverage and craw_htmp just create a symbolic linc like this::

        ln -s /Library/Framework/Python.Framework/Version/3.6/bin/craw_coverage /usr/local/bin/craw_coverage
        ln -s /Library/Framework/Python.Framework/Version/3.6/bin/craw_htmp /usr/local/bin/craw_htmp

    The documentation (html and pdf) is located in /Library/Framework/Python.Framework/Version/3.6/share/craw/


Installation from repository
----------------------------

Clone the project and install with the setup.py ::

    git clone https://gitlab.pasteur.fr/bneron/craw.git
    cd craw
    python3 setup.py install

.. note::
    Instead of installing craw you can directly use the scripts from the repository.
    You can also use the package without installing it.
    To do this, you have to export the **CRAW_HOME** environment variable.
    `CRAW_HOME` must point to the `src` directory of the project.
    Then you can use `craw_coverage` and `craw_htmp` scripts located in `bin` directory.

This project is documented using `sphinx <http://www.sphinx-doc.org/en/stable/index.html>`_.
So if you use a clone, you have to generate the documentation from the source.

The project come from with some unit and functional tests.
to test if everything work fine.

    cd $CRAW_HOME
    python3 tests/run_tests.py -vvv

.. _matplotlibrc:

matplotlib configuration
------------------------

`matplotlib` is a python library to create graphics.
`craw_htmp` use this library to generate heat map.
The two parameters to configure for craw is:

    * the backend
    * figure.dpi

backend
"""""""

matplolib lay on low level graphic library of youre computer. This library will determine
the graphical formats manage by matplotlib and  then by craw_htmp. Most of backend handle `png`
but some library like `Qt` can handle 'jpeg', 'eps', `pdf` ...

In your `matplolibrc` file you must define the backend for instance to use Qt5 ::

    backend: qt5agg

An example of matplolibrc file and all supported backend is available here:
http://matplotlib.org/users/customizing.html#a-sample-matplotlibrc-file

figure.dpi
""""""""""

It's not an essential option but `matplolib` and `craw_htmp` will produce better graphic (on screen)
if you configure `matplotlib` to the native resolution of your screen.
To know the resolution of your screen you can visit the following page https://www.infobyip.com/detectmonitordpi.php
and report the resolution (for 1 inch) in `matplotlibrc` file like: ::

    figure.dpi: 96

For full explanation on how to configure matplotlib read
http://matplotlib.org/users/customizing.html#the-matplotlibrc-file.
