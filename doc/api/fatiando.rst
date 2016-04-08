.. _fatiando:

========================
The ``fatiando`` package
========================

.. automodule:: fatiando
    :no-members:
    :no-inherited-members:


Geometric objects and meshes ``fatiando.mesher``
================================================
----

.. automodule:: fatiando.mesher
    :no-members:
    :no-inherited-members:

Classes
-------
.. currentmodule:: fatiando.mesher

.. autosummary::

    Polygon
    Square
    Prism
    PolygonalPrism
    Sphere
    Tesseroid
    SquareMesh
    PrismMesh
    PrismRelief
    TesseroidMesh
    PointGrid

Functions
----------
.. currentmodule:: fatiando.mesher

.. autosummary::

    extract
    vfilter
    vremove

Grids and irregularly spaced data ``fatiando.gridder``
=======================================================
----

.. automodule:: fatiando.gridder
    :no-members:
    :no-inherited-members:

Functions
----------
.. currentmodule:: fatiando.gridder

.. autosummary::

    regular
    scatter
    cut
    profile
    interp
    interp_at
    extrapolate_nans
    load_surfer
    spacing


Utility functions ``fatiando.utils``
====================================
----

.. automodule:: fatiando.utils
    :no-members:
    :no-inherited-members:

Classes
-------
.. currentmodule:: fatiando.utils

.. autosummary::

    SparseList

Functions
----------
.. currentmodule:: fatiando.utils

.. autosummary::

    normal
    gaussian
    gaussian2d
    safe_solve
    safe_dot
    safe_diagonal
    safe_inverse
    random_points
    circular_points
    connect_points
    si2mgal
    mgal2si
    si2eotvos
    eotvos2si
    si2nt
    nt2si
    sph2cart
    fromimage
    contaminate
    dircos
    ang2vec
    vecnorm
    vecmean
    vecstd
    sec2hms
    sec2year
    year2sec



.. toctree::
    :maxdepth: 1
    :hidden:

    gravmag.rst
    seismic.rst
    geothermal.rst
    mesher.rst
    gridder.rst
    vis.rst
    datasets.rst
    utils.rst
    constants.rst
    inversion.rst
