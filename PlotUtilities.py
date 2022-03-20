# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 21:15:40 2022

@author: rubin
"""
#Plotting utilities

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.patches import Ellipse

def hline(y,color):
    plt.axhline(y=y,color = color)
    
def vline(x,color):
    plt.axhline(x=x,color = color)


def MakeRectangle(x,y,width,height, color, fill):
    rect = plt.Rectangle((x,y),width,height, color = color, fill = fill)
    return(rect)
                         
    
def MakeCircle(x,y,radius, color, fill):
    circ = plt.Circle((x,y),radius, color = color, fill = fill)
    return(circ)


def MakeEllipse(x,y,width,height,color,fill):
    ellip = Ellipse(x,y,width,heigh,color = color, fill = fill)
    return(ellip)

def PlotPatch(patch,limits):
    fig,ax = plt.subplots()
    ax.add_patch(patch)
    ax.set_xlim((-limits,limits))
    ax.set_ylim((-limits,limits))
    
def AddPatch(patch):
    ax = plt.gca()
    ax.add_patch(patch)
    
def MakeFigure(title,xlabel,ylabel,ticksize,fontsize,titlesize):
    fig,ax = plt.subplots()
    ax.tick_params(axis = 'both',which = 'major', labelsize = ticksize)
    ax.tick_params(axis = 'both',which = 'minor', labelsize = ticksize)
    plt.xlabel(xlabel,fontsize = fontsize)
    plt.ylabel(ylabel,fontsize = fontsize)
    plt.title(title,fontsize = titlesize)