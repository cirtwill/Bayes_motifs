import sys
import os
import math
import random
from decimal import *
import numpy as np

#Pygrace libraries
from PyGrace.grace import Grace
from PyGrace.colors import ColorBrewerScheme
from PyGrace.dataset import SYMBOLS
from PyGrace.Extensions.panel import Panel,MultiPanelGrace
from PyGrace.Extensions.distribution import CDFGraph, PDFGraph
from PyGrace.Extensions.latex_string import LatexString, CONVERT
from PyGrace.drawing_objects import DrawText, DrawLine,DrawBox
from PyGrace.Extensions.colorbar import SolidRectangle, ColorBar
from PyGrace.axis import LINEAR_SCALE, LOGARITHMIC_SCALE
from PyGrace.Styles.el import ElGraph, ElLogColorBar

colors=ColorBrewerScheme('Blues', reverse=True)  # The blue is very beautiful but maybe harder to see.
colors.add_color(173,132,191,'Purple')
colors.add_color(34,90,34,'Green')
# colors.add_color(120,120,120,'grey')
# colors.add_color(255,125,125,'lightish_red')
# colors.add_color(200,200,200,'lightgrey')

def format_graph(graph,graphtype):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  graph.world.xmin=0
  graph.world.xmax=0.51
  graph.xaxis.tick.major=.1

  graph.world.ymin=0
  graph.world.ymax=.6000001
  graph.yaxis.tick.major=.2
  graph.yaxis.label.configure(text='Proportion',char_size=1,just=2,place='normal')
  graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=.75)


  graph.xaxis.ticklabel.configure(format='decimal',prec=2,char_size=.75)
  graph.xaxis.label.configure(text='Proportion basal resources',char_size=1,just=2,place='normal')

  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.4,minor_size=.5,place='both',major_linewidth=.5,minor_linewidth=1)
  graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.4,minor_size=.5,place='both',major_linewidth=.5,minor_linewidth=1)

  return graph

def total_dataset(nettype):
  if nettype=='chain':
    I=0.29775
    b=-0.08041
    c=-0.49806
    bc=0.05543
  elif nettype=='apparent':
    I=0.443506
    b=0.202398
    c=-0.596003
    bc=0.257850
  elif nettype=='direct':
    I=0.184305
    b=0.075896
    c=0.199148
    bc=-2.845068
  elif nettype=='omnivory':
    I=0.074440
    b=-0.197886
    c=0.894912
    bc=2.531783

  datadict={}
  for C in [0.02,0.06,0.1,0.14,0.18]:
    dataset=[]
    for B in [0.02,0.1,0.2,0.3,0.4,0.5]:
      y=I+B*b+C*c+B*C*bc
      dataset.append((B,y))
    datadict[C]=dataset

  return datadict

def populate_graph(graph,nettype):
  i=1
  datadict=total_dataset(nettype)
  for C in sorted(datadict):
    i+=1
    dats2=graph.add_dataset(datadict[C])
    dats2.symbol.configure(shape=0,color=i,fill_color=i)
    dats2.line.configure(color=i)

    dats1=graph.add_dataset(datadict[C])
    dats1.symbol.configure(color=i,fill_color=i)
    dats1.line.linestyle=0

    if lm=='apparent':
      dats1.legend=str(C)
  if lm=='apparent':
    graph.add_drawing_object(DrawText,text='Connectance',char_size=.75,x=.4,y=.37,just=2,loctype='world')
    graph.legend.configure(box_linestyle=0,char_size=.75,loc=(0.35,0.35),loctype='world')


  return graph

###############################################################################################
###############################################################################################
#
#           Assemble the plots
#
###############################################################################################
###############################################################################################

grace=MultiPanelGrace(colors=colors)
grace.add_label_scheme('dummy',['Apparent competition','Direct competition','Omnivory','3-species chain','Season length','F','G','H'])
grace.set_label_scheme('dummy')
for lm in ['apparent','direct','omnivory','chain']:
  graph=grace.add_graph(Panel)
  graph=format_graph(graph,'proportion')
  graph=populate_graph(graph,lm)
  graph.panel_label.configure(placement='ouc',char_size=1,dx=.03,dy=.01)

grace.multi(rows=2,cols=2,vgap=.03,hgap=.04)
grace.hide_redundant_labels()
grace.set_row_xaxislabel(label='Connectance',row=1,colspan=(None,None),char_size=1,perpendicular_offset=.05,just=2)
grace.set_col_yaxislabel(label='Proportion of total motifs',col=0,rowspan=(None,None),char_size=1,perpendicular_offset=.05,just=2)
# grace.set_row_xaxislabel(label='Day of year',row=2,colspan=(None,None),char_size=1,perpendicular_offset=.05)
for graph in grace.graphs:
  print graph.get_view()
grace.graphs[0].set_view(0.15,0.45,0.53,0.70)
grace.graphs[1].set_view(0.57,0.45,0.95,0.70)
grace.graphs[2].set_view(0.15,0.15,0.53,0.40)
grace.graphs[3].set_view(0.57,0.15,0.95,0.40)
grace.write_file('../manuscript/figures/motif_proportion_basal_lms.eps')

