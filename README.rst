Pipe Stress Infinity
====================

Pipe Stress Infinity (PSI) is an engineering design and analysis software used
to evaluate the structural behavior and stresses of piping systems to a variety
of different codes and standards.

* PSI on [PyPI]
* PSI [documentation]
* PSI [discord] server
* PSI [mailing list]
* PSI [issue tracker]

PSI has an active developer and user community. If you find a bug or a problem
with the documentation, please open an issue. Anyone is welcome to join our
discord server were a lot of the development discussion is going on. It's also
a great place to ask for help.

Some of the features of PSI are:

* Linear static analysis


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

If you're reading this README from a source distribution, you can install
PSI with:

..code::sh

    $ python setup.py install --user

You can also install the latest development version direct from Github using:

pip install https://github.com/denisgomes/psi/archive/master.zip

For local development install PSI in editable mode:

# with pip
pip install -e .
# with setup.py
python setup.py develop

There are no compilation steps during the installation; if you prefer, you can
simply add this directory to your PYTHONPATH and use PSI without installing
it. You can also copy PSI directly into your project folder.


Quickstart
----------

To start PSI in interactive mode type:

.. code:: sh

    $ psi   # start psi interpreter

PSI is a python interpreter with added functionality which allows for creating
and analyzing piping systems. Write the input file shown below using a text
editor and save it as **demo.py**:

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

    # define loading
    w1 = Weight('W1')
    p1 = Pressure('P1', 250)

    # define a  loadcase
    l1 = LoadCase('l1', 'ope', [w1, p1])

    # run the analysis
    mdl.analyze()

    # postprocess
    disp = Movements('r1', [l1])
    disp.to_screen()

Now run the file above to get the displacements at the nodes:

.. code:: sh

    $ psi demo.py   # run demo.py

To go directly into interacive model after running the model, use:

..code:: sh

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

PSI is developed by many individual volunteers, and there is no central
point of contact. If you have a question about developing with PSI, or you
wish to contribute, please join the mailing list or the discord server.

For license issues, please contact Denis Gomes, the primary author.


[PyPI]
[documentation]
[discord]: https://discord.gg/RZvjbAy
[mailing list]: https://groups.google.com/group/pipestressinfinity-users
[issue tracker]: https://github.com/denisgomes/psi/issues
