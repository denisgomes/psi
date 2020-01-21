Pipe Stress Infinity
====================

Pipe Stress Infinity (PSI) is an engineering design and analysis software used
to evaluate the structural behavior and stresses of piping systems to a variety
of different codes and standards.

PSI has an active developer and user community. If you find a bug or a problem
with the documentation, please open an issue. Anyone is welcome to join our
discord server were a lot of the development discussion is going on. It's also
a great place to ask for help.

Some of the features of PSI are:

* Modeling of piping components such as Runs, Bends, Reducers, Valves, and
  Flanges
* Ability to specify section and material data properties for different
  cross-sections
* Access to a variety of different support types such as Anchors, GlobalX,
  GlobalY, and GlobalZ including non-linear support capability
* Assign loads such as Weight, Pressure, Thermal, Wind and Seismic among others
* Linear combinations of loads and combinations of loads
* Stress evaluation based on B31.1 power piping code
* Linear static analysis
* Movements, support reactions and internal force results
* And more...


Requirements
------------

PSI supports Python 3.5 and above. The following is a list of libraries that it
depends on:

* jinja2
* numpy
* pint
* tqdm


Installation
------------

To install PSI, type the following command in the terminal:

.. code:: sh

    $ pip install psi --user


Install From Source
-------------------

If you're reading this README from a source distribution, you can install PSI
with:

.. code:: sh

    $ python setup.py install --user

You can also install the latest development version direct from Github using:

.. code:: sh

    $ pip install https://github.com/denisgomes/psi/archive/master.zip

For local development install PSI in editable mode:

.. code:: sh

    # with pip
    $ pip install -e .

    # with setup.py
    $ python setup.py develop

There are no compilation steps during the installation; if you prefer, you can
simply add this directory to your PYTHONPATH and use PSI without installing it.
You can also copy PSI directly into your project folder.


Quickstart
----------

To start PSI in interactive mode just type:

.. code:: sh

    $ psi

PSI provides an interpreter with added functionality for creating and analyzing
piping systems using Python scripts. Type the input file shown below in your
favorite text editor and save it as *demo.inp*:

.. literalinclude:: ../examples/demo.inp
    :language: python

Now run the file above to get the displacements at the nodes and the reaction
force at the anchor:

.. code:: sh

    $ psi demo.py > demo.out        # run demo.py and redirect to demo.out

Inspect the demo.out file to view the output:

.. literalinclude:: ../examples/demo.out

To go directly into interactive mode after running the model, use the -i
switch:

.. code:: sh

    $ psi -i demo.py > demo.out     # start the PSI interpreter


Contribution
------------

There are many different ways to contribute. You can promote PSI, fix bugs,
participate on the mailing list, etc. Please read the documentation_ to find
out more ways to contribute to the Community and Enterprise Editions.


Support
-------

The Community Edition of PSI is open source and can be used free of charge. If
the software adds value to your life i.e. you use it to do commercial work for
example, considered donating to support PSI's ongoing development. Also please
read the documentation_ about the feature rich Enterprise Edition.


Building Docs
-------------

PSI's documentation_ is hosted on `Read the Docs <https://readthedocs.org>`_.
The PDF_ version and other popular formats are also available on this website.
To build the docs, type the following from the project root directory:

.. code:: sh

    $ cd docs
    $ make clean && make html


Building Website
----------------

The PSI website (i.e. this site) consists of the project README.rst file hosted
on GitHub via GitHub Pages.  It lives on a branch called *gh-pages* and is
viewable at the `project website  <https://denisgomes.github.io/psi>`_.

To build and upload the website, type the following from the project root
directory:

.. code:: sh

    $ cd www
    $ make clean && make html   # for testing
    $ make github               # to publish


Testing
-------

.. code:: sh

    $ cd tests
    $ tox


Contact
-------

PSI is developed by many individual volunteers, and there is no central point
of contact. If you have a question about developing with PSI, or you wish to
contribute, please join the mailing list or the discord server.

For license questions, please contact `Denis Gomes`_, the primary author. Also
check out the project links below.

* PSI on PyPI_
* PSI documentation_
* PSI discord_ server
* PSI `mailing list`_
* PSI `issue tracker`_

.. _PSI: https://denisgomes.github.io/psi
.. _PyPI: https://pypi.org/
.. _documentation: https://pipe-stress-infinity.readthedocs.io
.. _PDF: https://readthedocs.org/projects/pipe-stress-infinity/downloads/pdf/latest
.. _discord: https://discord.gg/xnHnwbD
.. _mailing list: https://groups.google.com/group/pipestressinfinity-users
.. _issue tracker: https://github.com/denisgomes/psi/issues
.. _Denis Gomes: denis.mp.gomes@gmail.com
