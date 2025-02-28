# -*- coding: utf-8 -*-
"""
Created on Wed Aug 08 15:06:15 2018

@author: rh17872
"""

import scalebars as scale

def scale_bar(ax,sizex,sizey,labelx,labely,location,**kwargs) :
    scale.add_scalebar(ax,matchx =False, matchy = False, hidex = True, hidey=True, sizex = sizex, sizey = sizey, labelx = labelx, labely = labely, loc = location, **kwargs)