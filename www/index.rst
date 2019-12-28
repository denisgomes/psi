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
favorite text editor and save it as **demo.py**:

.. code:: python

    #! /usr/bin

    # parameters
    L = 10 * 12

    # create top level
    mdl = Model('demo')

    # define properties
    pipe1 = Pipe.from_file('pipe1', '10', '40')
    mat1 = Material.from_file('mat1', 'A53A', 'B31.1')

    # create geometry
    pt10 = Point(10)
    run20 = Run(20, L)

    # assign supports
    anc1 = Anchor('A1', 10)
    anc1.apply([run20])

    # define loads for operating case 1
    w1 = Weight('W1', 1)
    p1 = Pressure('P1', 1, 250)

    # define a loadcase
    l1 = LoadCase('l1', 'ope', [w1, p1], [1, 1])

    # run the analysis
    mdl.analyze()

    # postprocess
    disp = Movements('r1', [l1])
    disp.to_screen()

Now run the file above to get the displacements at the nodes:

.. code:: sh

    $ psi demo.py       # run demo.py

To go directly into interactive mode after running the model, use the -i
switch:

.. code:: sh

    $ psi -i demo.py    # run demo.py and start interpreter


Contribution
------------

Soon to come!


Building Docs
-------------

Soon to come!


Testing
-------

Soon to come!


Contact
-------

PSI is developed by many individual volunteers, and there is no central point
of contact. If you have a question about developing with PSI, or you wish to
contribute, please join the mailing list or the discord server.

For license questions, please contact `Denis Gomes`_, the primary author.


.. _PyPI: https://pypi.org/
.. _documentation: https://readthedocs.com/
.. _discord: https://discord.gg/RZvjbAy
.. _mailing list: https://groups.google.com/group/pipestressinfinity-users
.. _issue tracker: https://github.com/denisgomes/psi/issues
.. _Denis Gomes: denis.mp.gomes@gmail.com