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

def read_dispfile(dispfile):
  dispdict=[]

  f=open(dispfile,'r')
  for line in f:
    if line.split()[0]!='"Distance"':
      dist=float(line.split()[1][1:-1])
      S=int(line.split()[3][1:-1])
      C=float(line.split()[4][1:-1])
      dispdict.append((S,C,dist))
  f.close()
  return dispdict

def format_graph(graph,graphtype):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  graph.world.xmin=0
  graph.world.xmax=0.21
  graph.xaxis.tick.major=.05

  graph.world.ymin=0
  if graphtype=='prop':
    graph.world.ymax=0.25
    graph.yaxis.tick.major=.05
    graph.yaxis.ticklabel.configure(format='decimal',prec=2,char_size=.75)
  else:
    graph.world.ymax=0.4
    graph.yaxis.tick.major=0.1
    graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=.75)

  graph.xaxis.ticklabel.configure(format='decimal',prec=2,char_size=.75)
  graph.yaxis.label.configure(text='Motif profile dispersion',char_size=1,just=2,place='normal')
  graph.xaxis.label.configure(text='Connectance',char_size=1,just=2,place='normal')

  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.4,minor_size=.5,place='both',major_linewidth=.5,minor_linewidth=1)
  graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.4,minor_size=.5,place='both',major_linewidth=.5,minor_linewidth=1)

  return graph

def populate_graph(graph,dispdict,graphtype):
  for C in [0.02,0.06,0.1,0.14,0.18]:
    line=graph.add_dataset([(C,0),(C,1)])
    line.symbol.shape=0
    line.linewidth=.25

  for (S,C,disp) in dispdict:
    data=graph.add_dataset([(C+(float(S-70)/5000),disp)])
    data.symbol.configure(fill_color=S/10,color=S/10,fill_pattern=0)

  for S in [50,60,70,80,90,100]:
    temp=graph.add_dataset([(S,S)])
    temp.line.linestyle=0
    temp.symbol.configure(fill_color=S/10,color=S/10,fill_pattern=0)
    temp.legend=str(S)

  if graphtype=='heatmap':
    graph.add_drawing_object(DrawText,text='Size',char_size=.75,x=.115,y=.38,just=2,loctype='world')
    graph.legend.configure(box_linestyle=0,char_size=.75,loc=(0.105,0.375),loctype='world')
  else:
    graph.add_drawing_object(DrawText,text='Size',char_size=.75,x=.115,y=.235,just=2,loctype='world')
    graph.legend.configure(box_linestyle=0,char_size=.75,loc=(0.105,0.23),loctype='world')

  return graph

###############################################################################################
###############################################################################################
#
#           Assemble the plots
#
###############################################################################################
###############################################################################################

# dispfile='stat_analysis/motif_variability_SC.tsv'
# dispdict=read_dispfile(dispfile) # S, C, dict

# grace=MultiPanelGrace(colors=colors)
# # grace.add_label_scheme('dummy',['Start of flowering','Mass flowering','Inferred peak','End of flowering','Season length','F','G','H'])
# # grace.set_label_scheme('dummy')

# graph=grace.add_graph(Panel)
# graph=format_graph(graph,'heatmap')
# graph=populate_graph(graph,dispdict,'heatmap')
# graph.panel_label.configure(placement='iul',char_size=1,dx=.03,dy=.03)

# # grace.multi(rows=3,cols=2,vgap=.09,hgap=.04)
# # grace.hide_redundant_labels()
# # grace.set_row_xaxislabel(label='Degree (Number of interaction partners per species)',row=1,colspan=(None,None),char_size=1,perpendicular_offset=.05)
# # grace.set_row_xaxislabel(label='Day of year',row=2,colspan=(None,None),char_size=1,perpendicular_offset=.05)
# # for graph in grace.graphs:
#   # print graph.get_view()
# graph.set_view(0.15,0.15,0.95,0.65)
# grace.write_file('../manuscript/figures/motif_profile_dispersion.eps')


dispfile='stat_analysis/proportion_variability_SC.tsv'
dispdict=read_dispfile(dispfile) # S, C, dict

grace=MultiPanelGrace(colors=colors)

graph=grace.add_graph(Panel)
graph=format_graph(graph,'prop')
graph=populate_graph(graph,dispdict,'prop')
graph.panel_label.configure(placement='iul',char_size=1,dx=.03,dy=.03)

# grace.multi(rows=3,cols=2,vgap=.09,hgap=.04)
# grace.hide_redundant_labels()
# grace.set_row_xaxislabel(label='Degree (Number of interaction partners per species)',row=1,colspan=(None,None),char_size=1,perpendicular_offset=.05)
# grace.set_row_xaxislabel(label='Day of year',row=2,colspan=(None,None),char_size=1,perpendicular_offset=.05)
# for graph in grace.graphs:
  # print graph.get_view()
graph.set_view(0.15,0.15,0.95,0.65)
grace.write_file('../manuscript/figures/proportion_profile_dispersion.eps')

