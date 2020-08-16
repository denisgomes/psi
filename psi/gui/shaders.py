#
# Copyright Tristam Macdonald 2008.
#
# Distributed under the Boost Software License, Version 1.0
# (see http://www.boost.org/LICENSE_1_0.txt)
#

from pyglet.gl import *
from ctypes import *


class Shader:
    # vert, frag and geom take arrays of source strings
    # the arrays will be concattenated into one string by OpenGL
    def __init__(self, vert = [], frag = [], geom = []):
        # create the program handle
        self.handle = glCreateProgram()
        # we are not linked yet
        self.linked = False

        # create the vertex shader
        self.createShader(vert, GL_VERTEX_SHADER)
        # create the fragment shader
        self.createShader(frag, GL_FRAGMENT_SHADER)
        # the geometry shader will be the same, once pyglet supports the extension
        # self.createShader(frag, GL_GEOMETRY_SHADER_EXT)

        # attempt to link the program
        self.link()

    def createShader(self, strings, type):
        count = len(strings)
        # if we have no source code, ignore this shader
        if count < 1:
            return

        # create the shader handle
        shader = glCreateShader(type)

        # convert the source strings into a ctypes pointer-to-char array, and upload them
        # this is deep, dark, dangerous black magick - don't try stuff like this at home!
        src = (c_char_p * count)(*strings)
        glShaderSource(shader, count, cast(pointer(src), POINTER(POINTER(c_char))), None)

        # compile the shader
        glCompileShader(shader)

        temp = c_int(0)
        # retrieve the compile status
        glGetShaderiv(shader, GL_COMPILE_STATUS, byref(temp))

        # if compilation failed, print the log
        if not temp:
            # retrieve the log length
            glGetShaderiv(shader, GL_INFO_LOG_LENGTH, byref(temp))
            # create a buffer for the log
            buffer = create_string_buffer(temp.value)
            # retrieve the log text
            glGetShaderInfoLog(shader, temp, None, buffer)
            # print the log to the console
            print buffer.value
        else:
            # all is well, so attach the shader to the program
            glAttachShader(self.handle, shader);

    def link(self):
        # link the program
        glLinkProgram(self.handle)

        temp = c_int(0)
        # retrieve the link status
        glGetProgramiv(self.handle, GL_LINK_STATUS, byref(temp))

        # if linking failed, print the log
        if not temp:
            #   retrieve the log length
            glGetProgramiv(self.handle, GL_INFO_LOG_LENGTH, byref(temp))
            # create a buffer for the log
            buffer = create_string_buffer(temp.value)
            # retrieve the log text
            glGetProgramInfoLog(self.handle, temp, None, buffer)
            # print the log to the console
            print buffer.value
        else:
            # all is well, so we are linked
            self.linked = True

    def bind(self):
        # bind the program
        glUseProgram(self.handle)

    def unbind(self):
        # unbind whatever program is currently bound - not necessarily this program,
        # so this should probably be a class method instead
        glUseProgram(0)

    # upload a floating point uniform
    # this program must be currently bound
    def uniformf(self, name, *vals):
        # check there are 1-4 values
        if len(vals) in range(1, 5):
            # select the correct function
            { 1 : glUniform1f,
                2 : glUniform2f,
                3 : glUniform3f,
                4 : glUniform4f
                # retrieve the uniform location, and set
            }[len(vals)](glGetUniformLocation(self.handle, name), *vals)

    # upload an integer uniform
    # this program must be currently bound
    def uniformi(self, name, *vals):
        # check there are 1-4 values
        if len(vals) in range(1, 5):
            # select the correct function
            { 1 : glUniform1i,
                2 : glUniform2i,
                3 : glUniform3i,
                4 : glUniform4i
                # retrieve the uniform location, and set
            }[len(vals)](glGetUniformLocation(self.handle, name), *vals)

    # upload a uniform matrix
    # works with matrices stored as lists,
    # as well as euclid matrices
    def uniform_matrixf(self, name, mat):
        # obtian the uniform location
        loc = glGetUniformLocation(self.Handle, name)
        # uplaod the 4x4 floating point matrix
        glUniformMatrix4fv(loc, 1, False, (c_float * 16)(*mat))


edge_on = Shader(['''
varying vec4 vcolor;

float scale_factor = 0.05;

void main()
{
    vec4 vnormal;
    vec4 vposition;

    vcolor = gl_Color;  // pass vertex color to fragment

    // flipped normal in world csys
    vnormal = gl_ModelViewProjectionMatrix * vec4(gl_Normal, 0.0);

    // vertex world csys
    vposition = gl_ModelViewProjectionMatrix * gl_Vertex;

    gl_Position = vposition + scale_factor*vnormal;
}
'''], ['''
varying vec4 vcolor;

void main (void)
{
    if (gl_FrontFacing) discard; //discard front-facing
    else gl_FragColor = vec4(0.,0.,0.,1.); //draw rest black
}
'''])


# create our Phong Shader by Jerome GUINOT aka 'JeGX' - jegx [at] ozone3d [dot] net
# see http://www.ozone3d.net/tutorials/glsl_lighting_phong.php

phong = Shader(['''
varying vec3 normal, lightDir, eyeVec;
varying vec4 vcolor;

void main()
{
    normal = gl_NormalMatrix * gl_Normal;

    vec3 vVertex = vec3(gl_ModelViewMatrix * gl_Vertex);

    lightDir = vec3(gl_LightSource[0].position.xyz - vVertex);
    eyeVec = -vVertex;

    gl_Position = ftransform();


    vcolor = gl_Color;  // pass vertex color to fragment
}
'''], ['''
varying vec3 normal, lightDir, eyeVec;
varying vec4 vcolor;

void main (void)
{
    // ambient
    vec4 final_color =
    (gl_FrontLightModelProduct.sceneColor * gl_FrontMaterial.ambient) +
    (gl_LightSource[0].ambient * gl_FrontMaterial.ambient);

    vec3 N = normalize(normal);
    vec3 L = normalize(lightDir);

    float lambertTerm = dot(N,L);

    if(lambertTerm > 0.0)
    {
        // diffuse
        final_color += gl_LightSource[0].diffuse *
                       gl_FrontMaterial.diffuse *
                       lambertTerm;

        vec3 E = normalize(eyeVec);
        vec3 R = reflect(-L, N);

        // specular
        float specular = pow( max(dot(R, E), 0.0),
                         gl_FrontMaterial.shininess);
        final_color += gl_LightSource[0].specular *
                       gl_FrontMaterial.specular *
                       specular;
    }

    gl_FragColor = final_color * vcolor;
}
'''])

depth_clamp = Shader(['''
varying float z;
varying vec4 vcolor;
float near, far;

void main()
{
    gl_Position = ftransform();

    // transform z to window coordinates
    z = gl_Position.z / gl_Position.w;

    far = gl_DepthRange.far;
    near = gl_DepthRange.near;
    z = 0.5*z*(far-near) + 0.5*(far+near);

    // prevent z-clipping
    gl_Position.z = 0.0;

    vcolor = gl_Color;  // pass vertex color to fragment
}
'''], ['''
varying float z;
varying vec4 vcolor;

void main()
{
    gl_FragColor = vec4(vec3(vcolor), 1.0);
    gl_FragDepth = clamp(z, 0.0, 1.0);
}
'''])
