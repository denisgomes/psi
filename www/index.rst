===================
OpenPipeStress: PSI
===================

.. toctree::
   :maxdepth: 2
   :hidden:

   pages/about/about.rst

Piping Input
============

.. jupyter-execute::

    from pythreejs import *
    from IPython.display import display

    from ipywidgets import widgets
    from ipywidgets import HBox, VBox

    from ipysheet import sheet, cell


    mesh1 = Mesh(SphereBufferGeometry(20, 16, 16), MeshPhysicalMaterial(color='red'), position=[-20, 0, 0])
    mesh2 = Mesh(SphereBufferGeometry(20, 16, 16), MeshPhysicalMaterial(color='green'), position=[20, 0, 0])

    view_width = 640
    view_height = 400
    camera = CombinedCamera(position=[0, 0, 60], width=view_width, height=view_height)

    key_light = PointLight(position=[-100, 100, 100])
    ambient_light = AmbientLight(intensity=0.4)
    scene = Scene(children=[mesh1, mesh2, key_light, ambient_light, camera])
    renderer = Renderer(scene=scene, camera=camera, controls=[OrbitControls(controlling=camera)],
                        width=view_width, height=view_height)

    accordion = widgets.Accordion(children=[widgets.Text(),
                                  widgets.Text(),
                                  widgets.Text(),
                                  widgets.Text(),
                                  widgets.Text(),
                                  widgets.Text(),
                                  ])
    accordion.set_title(0, 'Points')
    accordion.set_title(1, 'Section')
    accordion.set_title(2, 'Material')
    accordion.set_title(3, 'Elements')
    accordion.set_title(4, 'Loads')
    accordion.set_title(5, 'Loadcases')

    vb = VBox((renderer, accordion))
    vb
