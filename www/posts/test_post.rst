Testing
=======
.. thebe-button:: Start Jupyter Kernel


.. jupyter-execute::

    import numpy as np
    from matplotlib import pyplot
    %matplotlib inline

    x = np.linspace(1E-3, 2 * np.pi)

    pyplot.plot(x, np.sin(x) / x)
    pyplot.plot(x, np.cos(x))
    pyplot.grid()

.. image:: ./test_post/python.jpg
  :width: 200
  :alt: Alternative text

.. post:: Dec 27, 2021
   :tags: calc
   :category: Calc
   :author: me
