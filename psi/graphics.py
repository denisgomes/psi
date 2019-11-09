"""Display the model on a browser"""

from psi.elements import *


def render_model(model):
    for element in model.elements:
        if element.type == "Run":
            draw_run(element)


def draw_run(run):
    pt1 = run.from_point
    pt2 = run.to_point

    pos = pt1
    axis = vec(pt2 - pt1)
    radius = run.section.od / 2

    cylinder(pos=pos, axis=axis, radius=radius)
