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


Installation
============

Installation from package
-------------------------

Using pip ::

    pip install craw-x.x.x.tar.gz

Do not forget to configure the `matplotlib` backend, specially if you use virtualenv.
see http://matplotlib.org/users/customizing.html#the-matplotlibrc-file for more explanation.

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
