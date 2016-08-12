===============
AnechoDB_Access
===============

It is a library used to connect to a specific database, download data stored in it and make some simple calculation. The data are beam patterns obtained from measurements in two anechoic chambers saved on the database as HDF5 files.

The package is divided in two distinct modules called **connection** and **computation** that have different tasks.

*************
connection.py
*************
This module is a class with useful function to establish a connection with the chosen database. 
The database is structered in a way that each *beam* identifier (the data user is looking for) is linked to a *measurements* page with information about the measure and the links to  *projects* and *instruments* pages.
With **connection** is possible to find the desired beam identifier through a search in one of those pages and finally download the data that is converted from a .h5 file to a Python dict variable preserving the same structure

**************
computation.py
**************
This module has some function used to perform simple (but useful) calculation from the data previously obtained from **connection**. Till now this module has only two functions, one to compute mean and variance of the data and the other to perform normalization and centering of beam patterns, but more will be added in the future.

Example of usage
----------------
After being installed here is a classic way to use this package.

.. code-block:: python

    >>> c = share_belen.connection.Connection(Host)
    >>> i_m = c.search_meas_by_instruments('Instrument To Search')
    >>> #More than one measurement can be linked at the same instrument
    >>> i_b = c.search_beam_by_meas(i_1[0])
    >>> #More than one beam can be linked at the same measurement
    >>> b = c.get_beam_in_dict_by_id(i_b[0])
    >>> b_c = share_belen.computation.make_beam_meanvar(b)
    >>> b_c_2 = share_belen.computation.center_norm_beam(b_c)

Requirements
------------

* `Python <http://www.python.org>`_ (tested with version >=3.3)
* `h5py <http://www.h5py.org/>`_