Pipe Stress Infinity
====================

Pipe Stress Infinity (PSI) is an engineering design and analysis software used
to evaluate the structural behavior and stresses of piping systems to a variety
of different codes and standards.

* PSI on PyPI_
* PSI documentation_
* PSI discord_ server
* PSI `mailing list`_
* PSI `issue tracker`_

PSI has an active developer and user community. If you find a bug or a problem
with the documentation, please open an issue. Anyone is welcome to join our
discord server were a lot of the development discussion is going on. It's also
a great place to ask for help.

Some of the features of PSI are:

* Linear static analysis
* Open source and free to use (consider supporting)


Requirements
------------

PSI supports Python 3.5 and above. The following is a list of libraries that
it depends on:

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

    $ psi demo.py > results.out     # run demo.py and redirect to results.out

Inspect the results.out file to view the output:

.. literalinclude:: ../examples/results.out

To go directly into interactive mode after running the model, use the -i
switch:

.. code:: sh

    $ psi -i demo.py > results.out  # start the PSI interpreter


Contribution
------------

Coming soon!


Support
-------

Coming soon!


Building Docs
-------------

PSI's documentation_ is hosted on the `Read the Docs
<https://readthedocs.org>`_ website.


Building Website
----------------

The PSI website (i.e. this site) consists of the project README.rst file hosted
on GitHub via GitHub Pages.  It lives on a branch called *gh-pages* and is
viewable at the `project website  <https://denisgomes.github.io/psi>`_.

To build and upload the website, type the following for the project root
directory:

.. code:: sh

    $ cd www
    $ make html     # for testing
    $ make github   # to publish


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

For license questions, please contact `Denis Gomes`_, the primary author.

.. _PSI: https://denisgomes.github.io/psi
.. _PyPI: https://pypi.org/
.. _documentation: https://pipe-stress-infinity.readthedocs.io
.. _discord: https://discord.gg/xnHnwbD
.. _mailing list: https://groups.google.com/group/pipestressinfinity-users
.. _issue tracker: https://github.com/denisgomes/psi/issues
.. _Denis Gomes: denis.mp.gomes@gmail.com
