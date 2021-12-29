Testing
=======
.. thebe-button:: Interact!


.. jupyter-execute::

    %matplotlib widget
    import ipywidgets as widgets
    import matplotlib.pyplot as plt
    import numpy as np

    x = np.linspace(0,10)

    def sine_func(x, w, amp):
        return amp*np.sin(w*x)

    @widgets.interact(w=(0, 4, 0.25), amp=(0, 4, .1))
    def update(w = 1, amp = 1):
        plt.clf()
        plt.ylim(-4, 4)
        plt.plot(x, sine_func(x, w, amp))

.. image:: ./test_post/python.jpg
  :width: 200
  :alt: Alternative text

.. post:: Dec 27, 2021
   :tags: calc
   :category: Tools
   :author: me
