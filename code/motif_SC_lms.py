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

colors=ColorBrewerScheme('Blues')  # The blue is very beautiful but maybe harder to see.
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
  graph.world.xmax=0.21
  graph.xaxis.tick.major=.05

  graph.world.ymin=0
  if graphtype=='total':
    graph.world.ymax=50000
    graph.yaxis.tick.major=10000
    graph.yaxis.label.configure(text='Total motifs',char_size=1,just=2,place='normal')
    graph.yaxis.ticklabel.configure(format='scientific',prec=1,char_size=.75)
  else:
    graph.world.ymax=.60000001
    graph.yaxis.tick.major=.2
    graph.yaxis.label.configure(text='Proportion',char_size=1,just=2,place='normal')
    graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=.75)


  graph.xaxis.ticklabel.configure(format='decimal',prec=2,char_size=.75)
  graph.xaxis.label.configure(text='Connectance',char_size=1,just=2,place='normal')

  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.4,minor_size=.5,place='both',major_linewidth=.5,minor_linewidth=1)
  graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.4,minor_size=.5,place='both',major_linewidth=.5,minor_linewidth=1)

  return graph

def total_dataset(nettype):
  if nettype=='total':
    I=-2.846e+03
    s=4.662e+01
    c=-2.122e+05
    sc=4.716e+03
  elif nettype=='chain':
    I=0.2617064
    s=0.0001775
    c=-0.1566627
    sc=-0.0032675
  elif nettype=='apparent':
    I=5.409e-01
    s=-5.058e-04
    c=-1.001e+00
    sc=2.590e-03
  elif nettype=='direct':
    I=1.926e-01
    s=-3.563e-05
    c=-5.556e-01
    sc=5.917e-03
  elif nettype=='omnivory':
    I=4.824e-03
    s=3.639e-04
    c=1.714e+00
    sc=-5.239e-03

  datadict={}
  for S in [50,60,70,80,90,100]:
    dataset=[]
    for C in [0.02,0.06,0.1,0.14,0.18]:
      y=I+S*s+C*c+S*C*sc
      dataset.append((C,y))
    datadict[S]=dataset

  return datadict

def populate_graph(graph,nettype):
  if nettype=='total':
    datadict=total_dataset('total')
    for S in sorted(datadict):
      dats=graph.add_dataset(datadict[S])
      dats.symbol.configure(color=S/10,fill_color=S/10)
      dats.line.configure(color=S/10)
      dats.legend=str(S)

    graph.add_drawing_object(DrawText,text='Size',char_size=.75,x=.027,y=40500,just=2,loctype='world')
    graph.legend.configure(box_linestyle=0,char_size=.75,loc=(0.02,40000),loctype='world')
  else:
    datadict=total_dataset(nettype)
    for S in sorted(datadict):
      dats2=graph.add_dataset(datadict[S])
      dats2.symbol.configure(shape=0,color=S/10,fill_color=S/10)
      dats2.line.configure(color=S/10)

      dats1=graph.add_dataset(datadict[S])
      dats1.symbol.configure(color=S/10,fill_color=S/10)
      dats1.line.linestyle=0

      if lm=='omnivory':
        dats1.legend=str(S)
    if lm=='omnivory':
      graph.add_drawing_object(DrawText,text='Size',char_size=.75,x=.027,y=.53,just=2,loctype='world')
      graph.legend.configure(box_linestyle=0,char_size=.75,loc=(0.005,0.52),loctype='world')


  return graph

###############################################################################################
###############################################################################################
#
#           Assemble the plots
#
###############################################################################################
###############################################################################################


grace=MultiPanelGrace(colors=colors)
# grace.add_label_scheme('dummy',['Start of flowering','Mass flowering','Inferred peak','End of flowering','Season length','F','G','H'])
# grace.set_label_scheme('dummy')
for lm in ['mot_total']:
  graph=grace.add_graph(Panel)
  graph=format_graph(graph,'total')
  graph=populate_graph(graph,'total')
  graph.panel_label.configure(placement='iul',char_size=1,dx=.03,dy=.03)

# grace.multi(rows=3,cols=2,vgap=.09,hgap=.04)
# grace.hide_redundant_labels()
# grace.set_row_xaxislabel(label='Degree (Number of interaction partners per species)',row=1,colspan=(None,None),char_size=1,perpendicular_offset=.05)
# grace.set_row_xaxislabel(label='Day of year',row=2,colspan=(None,None),char_size=1,perpendicular_offset=.05)
# for graph in grace.graphs:
  # print graph.get_view()
graph.set_view(0.15,0.15,0.95,0.65)
grace.write_file('../manuscript/figures/motif_context_lms.eps')


grace=MultiPanelGrace(colors=colors)
grace.add_label_scheme('dummy',['Apparent competition','Direct competition','Omnivory','3-species chain','Season length','F','G','H'])
grace.set_label_scheme('dummy')
for lm in ['apparent','direct','omnivory','chain']:
  graph=grace.add_graph(Panel)
  graph=format_graph(graph,'proportion')
  graph=populate_graph(graph,lm)
  graph.panel_label.configure(placement='ouc',char_size=1,dx=.03,dy=.01)

grace.multi(rows=2,cols=2,vgap=.05,hgap=.04)
grace.hide_redundant_labels()
grace.set_col_yaxislabel(label='Proportion of total motifs',col=0,rowspan=(None,None),char_size=1,perpendicular_offset=.06,just=2)
# grace.set_row_xaxislabel(label='Day of year',row=2,colspan=(None,None),char_size=1,perpendicular_offset=.05)
for graph in grace.graphs:
  print graph.get_view()
grace.graphs[0].set_view(0.15, 0.6010588235294116,0.53,0.95)
grace.graphs[1].set_view(0.57, 0.6010588235294116,0.95,0.95)
grace.graphs[2].set_view(0.15,0.18211764705882338,0.53,0.5310588235294116)
grace.graphs[3].set_view(0.57,0.18211764705882338,0.95,0.5310588235294116)
# grace.add_drawing_object(DrawText,text='Connectance',loctype='view',x=.55,y=0.08,char_size=1,just=2)
grace.write_file('../manuscript/figures/motif_proportion_lms.eps')

