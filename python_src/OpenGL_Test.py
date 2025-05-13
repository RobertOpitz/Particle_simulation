#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 11:15:36 2019

@author: robertopitz
"""

import pygame
from pygame.locals import *
from OpenGL.GL import *
from OpenGL.GLU import *
from math import isnan
import struct
import argparse

verticies = (( 1, -1, -1),
             ( 1,  1, -1),
             (-1,  1, -1),
             (-1, -1, -1),
             ( 1, -1,  1),
             ( 1,  1,  1),
             (-1, -1,  1),
             (-1,  1,  1))

edges = ((0,1),
         (0,3),
         (0,4),
         (2,1),
         (2,3),
         (2,7),
         (6,3),
         (6,4),
         (6,7),
         (5,1),
         (5,4),
         (5,7))

# construct the argument parse and parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required = True,
                help = "path to input particle file")
args = vars(ap.parse_args())

def get_data(filePath, display, nb_of_bytes = 4):

    if nb_of_bytes == 4:
        float_format_descriptor = "f"
        integer_format_descriptor = "i"
    elif nb_of_bytes == 8:
        float_format_descriptor = "d"
        integer_format_descriptor = "l"

    def read(f, format_descriptor):
        byte = f.read(nb_of_bytes)
        if not byte:
            return None
        value = struct.unpack(format_descriptor, byte)[0]
        return value

    with open(filePath,"rb") as f:
        # nb records of first block
        a = read(f, integer_format_descriptor)
        print(a)
        # nb records of second block
        b = read(f, integer_format_descriptor)
        print(b)

        size_x = 0.5 * read(f, float_format_descriptor)
        size_y = 0.5 * read(f, float_format_descriptor)
        size_z = 0.5 * read(f, float_format_descriptor)
        size = [size_x, size_y, size_z]

        # jump to the important parts of the file
        for _ in range(a+b-5):
            byte = f.read(nb_of_bytes)

        data = []
        while True:
            p = read(f, float_format_descriptor)
            if isnan(p):
                break
            print('time point: ', p)

            p = read(f, integer_format_descriptor)
            print('entries: ', p)

            pos = []
            for _ in range(int(p/7)):
                # read particle id (long int) to actually skip it
                p = read(f, integer_format_descriptor)
                # read particle positions (x,y,z) as float
                temp = []
                for i in range(3): #['x', 'y', 'z']:
                    p = read(f, float_format_descriptor)
                    temp.append(p/size[i] - 1.0)
                pos.append(temp)
                # read particle velocities to actually skip this entries
                for _ in range(3):#['vx', 'vy', 'vz']:
                    p = read(f, float_format_descriptor)

            data.append(pos)

    return data

def Cube():
    glBegin(GL_LINES)
    for edge in edges:
        for vertex in edge:
            glVertex3fv(verticies[vertex])
    glEnd()

def Points(data):
    #glPointsize( 10 )
    glBegin(GL_POINTS)
    for point in data:
        glVertex3f(point[0],
                   point[1],
                   point[2])
    glEnd()

def main(display,data):
    pygame.init()
    pygame.display.set_mode(display, DOUBLEBUF|OPENGL)

    gluPerspective(30, (display[0]/display[1]), 0.1, 50.0 )

    glTranslatef(0.0, 0.0, -5.0)

    glRotatef(22.5, 0, 1, 0)

    #running = True
    while True:
        for i in range(len(data)):
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    pygame.quit()
                    quit()

            glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
            # rotate the view
            glRotatef(1, 0, 1, 0)
            # draw the cube. It represents the simulation area
            Cube()
            # draw all particles
            Points(data[i])

            pygame.display.flip()



display = (800,600)
data = get_data(args["input"], display)

main(display,data)
