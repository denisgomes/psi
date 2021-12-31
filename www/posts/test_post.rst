Examples
========
.. jupyter-kernel:: python3

.. thebe-button:: Activate

.. jupyter-execute::

    import ipysheet
    sheet = ipysheet.sheet()
    sheet

.. jupyter-execute::

    from pythreejs import *
    from IPython.display import display


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
    display(renderer)


.. jupyter-execute::

    import plotly.graph_objs as go
    import plotly.offline as py

    import numpy as np
    from ipywidgets import interactive, HBox, VBox

    py.init_notebook_mode()

    x = y = np.arange(-5, 5, 0.1)
    yt = x[:, np.newaxis]
    z = np.cos(x * yt) + np.sin(x * yt) * 2

    f = go.FigureWidget(
        data=[
            go.Surface(z=z, x=x, y=y,
                    colorscale='Viridis')],
        layout=go.Layout(scene=go.layout.Scene(
            camera=go.layout.scene.Camera(
                up=dict(x=0, y=0, z=1),
                center=dict(x=0, y=0, z=0),
                eye=dict(x=1.25, y=1.25, z=1.25))
        ))
    )


    def update_z(frequency):
        f.data[0].z = np.cos(x * yt * frequency / 10.0) + np.sin(x * yt * frequency / 10.0) * 2


    freq_slider = interactive(update_z, frequency=(1, 50, 0.1))
    vb = VBox((f, freq_slider))
    vb.layout.align_items = 'center'
    vb


.. jupyter-execute::

    import plotly.graph_objects as go
    import numpy as np
    X, Y, Z = np.mgrid[-8:8:40j, -8:8:40j, -8:8:40j]
    values = np.sin(X*Y*Z) / (X*Y*Z)

    fig = go.Figure(data=go.Volume(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=values.flatten(),
        isomin=0.1,
        isomax=0.8,
        opacity=0.1, # needs to be small to see through all surfaces
        surface_count=17, # needs to be a large number for good volume rendering
        ))
    fig.show()

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

.. jupyter-execute::

    from mpl_toolkits import mplot3d
    import numpy as np
    import matplotlib.pyplot as plt
    from ipywidgets import widgets

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    z = np.linspace(0, 1, 100)
    x = z * np.sin(20 * z)
    y = z * np.cos(20 * z)
    ax.plot3D(x, y, z, 'gray')
    ax.set_title('3D line plot')
    fig

    accordion = widgets.Accordion(children=[widgets.IntSlider(),
                                  widgets.Text()])
    accordion.set_title(0, 'Slider')
    accordion.set_title(1, 'Text')
    accordion

.. jupyter-execute::

    from ipywidgets import Button, GridBox, Layout, ButtonStyle

    header  = Button(description='Header',
                    layout=Layout(width='auto', grid_area='header'),
                    style=ButtonStyle(button_color='lightblue'))
    main    = Button(description='Main',
                    layout=Layout(width='auto', grid_area='main'),
                    style=ButtonStyle(button_color='moccasin'))
    sidebar = Button(description='Sidebar',
                    layout=Layout(width='auto', grid_area='sidebar'),
                    style=ButtonStyle(button_color='salmon'))
    footer  = Button(description='Footer',
                    layout=Layout(width='auto', grid_area='footer'),
                    style=ButtonStyle(button_color='olive'))

    GridBox(children=[header, main, sidebar, footer],
            layout=Layout(
                width='50%',
                grid_template_rows='auto auto auto',
                grid_template_columns='25% 25% 25% 25%',
                grid_template_areas='''
                "header header header header"
                "main main . sidebar "
                "footer footer footer footer"
                ''')
        )

.. jupyter-execute::
    :hide-code:

    from ipywidgets import Layout, Button, Box, FloatText, Textarea, Dropdown, Label, IntSlider

    form_item_layout = Layout(
        display='flex',
        flex_flow='row',
        justify_content='space-between'
    )

    form_items = [
        Box([Label(value='Age of the captain'), IntSlider(min=40, max=60)], layout=form_item_layout),
        Box([Label(value='Egg style'),
            Dropdown(options=['Scrambled', 'Sunny side up', 'Over easy'])], layout=form_item_layout),
        Box([Label(value='Ship size'),
            FloatText()], layout=form_item_layout),
        Box([Label(value='Information'),
            Textarea()], layout=form_item_layout)
    ]

    form = Box(form_items, layout=Layout(
        display='flex',
        flex_flow='column',
        border='solid 2px',
        align_items='stretch',
        width='50%'
    ))
    form


.. jupyter-execute::

    from ipywidgets import Button, HBox, VBox

    words = ['correct', 'horse', 'battery', 'staple']
    items = [Button(description=w) for w in words]
    left_box = VBox([items[0], items[1]])
    right_box = VBox([items[2], items[3]])
    HBox([left_box, right_box])

.. jupyter-execute::

    from ipywidgets import Button, Layout

    b = Button(description='(50% width, 80px height) button',
            layout=Layout(width='50%', height='80px'))
    b

.. jupyter-execute::
    :hide-output:

    # importing mplot3d toolkits, numpy and matplotlib
    from mpl_toolkits import mplot3d
    import numpy as np
    import matplotlib.pyplot as plt

    fig = plt.figure()

    # syntax for 3-D projection
    ax = plt.axes(projection ='3d')

    # defining all 3 axes
    z = np.linspace(0, 1, 100)
    x = z * np.sin(25 * z)
    y = z * np.cos(25 * z)

    # plotting
    ax.plot3D(x, y, z, 'green')
    ax.set_title('3D line plot geeks for geeks')
    plt.show()

.. jupyter-execute::
    :hide-output:

    from IPython.display import display
    from pythreejs import (ParametricGeometry, Mesh, PerspectiveCamera, Scene,
                        MeshLambertMaterial, DirectionalLight, AmbientLight,
                        Renderer, OrbitControls, PerspectiveCamera)

    f = """
    function f(origu, origv, out) {
        // scale u and v to the ranges I want: [0, 2*pi]
        var u = 2*Math.PI*origu;
        var v = 2*Math.PI*origv;

        var x = Math.sin(u);
        var y = Math.cos(v);
        var z = Math.cos(u+v);

        out.set(x,y,z);
    }
    """
    surf_g = ParametricGeometry(func=f, slices=16, stacks=16)

    surf = Mesh(geometry=surf_g, material=MeshLambertMaterial(color='green', side='FrontSide'))
    surf2 = Mesh(geometry=surf_g, material=MeshLambertMaterial(color='yellow', side='BackSide'))
    c = PerspectiveCamera(position=[5, 5, 3], up=[0, 0, 1],
                        children=[DirectionalLight(color='white',
                                                    position=[3, 5, 1],
                                                    intensity=0.6)])
    scene = Scene(children=[surf, surf2, c, AmbientLight(intensity=0.5)])
    renderer = Renderer(camera=c, scene=scene, controls=[OrbitControls(controlling=c)], width=400, height=400)
    display(renderer)

.. jupyter-execute::

    from ipywidgets import interact, interactive, fixed, interact_manual
    import ipywidgets as widgets

    def f(x):
        return x

    interact(f, x=10);

.. jupyter-execute::
    :hide-output:

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

.. jupyter-execute::

    import ipywidgets as w
    from IPython.display import display

    a = w.IntSlider()
    b = w.IntText()
    w.jslink((a, 'value'), (b, 'value'))
    display(a, b)

.. image:: ./test_post/python.jpg
  :width: 200
  :alt: Alternative text

.. post:: Dec 27, 2021
   :tags: calc
   :category: Toolbox
   :author: me
