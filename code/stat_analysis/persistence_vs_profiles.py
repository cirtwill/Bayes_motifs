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

colors=ColorBrewerScheme('Greys',n=10)  # The blue is very beautiful but maybe harder to see.
# Match Anna's ggplot colours
colors.add_color(158,101,184,'Rpurple')
colors.add_color(145,148,215,'Rblue')
colors.add_color(42,120,142,'Rgreen1')
colors.add_color(34,168,132,'Rgreen2')
colors.add_color(122,209,81,'Rgreen3')
colors.add_color(253,231,37,'Ryellow')
colors.add_color(254,224,139,'Ryellow2')
colors.add_color(253,174,97,'Rorange1')
colors.add_color(244,109,67,'Rorange2')
colors.add_color(213,62,79,'Rred1')
colors.add_color(158,1,66,'Rred2')

mot_ranges={
  'omnivory':(0.005952381,0.345547537),
  'apparent':(0.2499348,0.7647059),
  'direct':(0.05010438,0.35878963),
  'chain':(0.1088861,0.4125000),
}

scales={ # center, scale
  'omnivory':(0.1641884,0.07895101),
  'apparent':(0.4222605,0.06446181),
  'direct':(0.1787041,0.0407462),
  'chain':(0.234847,0.04233961),
  'Disturbance':(0.3,0.1366298)
}

lms={ # Intercept, motif, disturbance, interaction
  'omnivory':(0.5743430145,-0.0140077658,-0.1385163991,-0.0005586564),
  'apparent':(0.574343015,0.011206293,-0.138516399,0.000589999),
  'direct':(0.5743430145,0.0033850297,-0.1385163991,0.0001809503),
  'chain':(5.743430e-01,5.801241e-03,-1.385164e-01,-3.067919e-05)  
}

motif_names={'omnivory':'Omnivory','apparent':'Apparent competition','direct':'Direct competition','chain':'Three-species chain'}

# Each point will be mean of one network
def read_datafile(datafile):
  netprops={}
  f=open(datafile,'r')
  for line in f:
    if line.split()[0]!='"ID"':
      ID=line.split()[0]
      chain=float(line.split()[1])
      apparent=float(line.split()[2])
      direct=float(line.split()[3])
      omnivory=float(line.split()[4])
      persistence=float(line.split()[5])

      netprops[ID]={'persisence':persistence,
      'chain':chain,'dircomp':direct,'omnivory':omnivory,'appcomp':apparent}

  f.close()
  return netprops

def format_graph(graph,simple):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  graph.world.xmin=0
  graph.world.ymin=0

  if simple=='profile':
    graph.xaxis.label.text='Proportion omnivory'
    graph.yaxis.label.text='Proportion other motif'
    graph.xaxis.tick.major=.1
    graph.yaxis.tick.major=.2
    graph.world.xmax=0.35
    graph.world.ymax=.6000001
  elif simple=='persistence':
    graph.world.ymin=.1
    graph.world.ymax=.9
    graph.world.xmax=.8
    graph.xaxis.tick.major=.2
    graph.yaxis.tick.major=.2
    graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=.75)
    graph.yaxis.label.text='Mean persistence'
    graph.xaxis.label.text='Proportion of network profile'


  graph.xaxis.ticklabel.configure(format='decimal',prec=1,char_size=1)
  graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=1)
  graph.xaxis.label.configure(char_size=1.5,just=2,place='normal')

  graph.yaxis.label.configure(char_size=1.5,just=2,place='normal')
  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.4,place='both',major_linewidth=.5,minor_linewidth=1)
  graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.4,place='both',major_linewidth=.5,minor_linewidth=1)


  return graph

def populate_persgraph(graph,motif):

  for dist in [0.1,0.5]:
    if motif=='apparent':
      j='Rpurple'
      sty=2
    elif motif=='chain':
      j='Rblue'
      sty=3
    elif motif=='direct':
      j='Rgreen2'
      sty=8
    else:
      j='Rorange1'
      sty=1

    # for dist in [0.1,0.18,0.26,0.34,0.42,0.5]:
    scaldist=(dist-scales['Disturbance'][0])/scales['Disturbance'][1]

    dats=[]

    minx=mot_ranges[motif][0]
    maxx=mot_ranges[motif][1]
    for x in [minx,maxx]:
      scalx=(x-scales[motif][0])/scales[motif][1]
      y=lms[motif][0]+scalx*lms[motif][1]+scaldist*lms[motif][2]+scalx*scaldist*lms[motif][3]

      dats.append((x,y))

    data=graph.add_dataset(dats)
    data.symbol.shape=0
    data.line.configure(linestyle=sty,linewidth=4,color=j)

    if dist==0.1:
      data.legend=motif_names[motif]

  return graph

###############################################################################################
###############################################################################################
#
#           Assemble the plots
#
###############################################################################################
###############################################################################################

# # Far too many species to see anything. Going to apply stats.
# datafile='motif_proportions_persistence.tsv'
# netprops=read_datafile(datafile)

grace=Grace(colors=colors)

graph2=grace.add_graph()
graph2=format_graph(graph2,'persistence')
for mot in ['omnivory','chain','apparent','direct']:
  graph2=populate_persgraph(graph2,mot)

graph2.add_drawing_object(DrawText,x=0.5,y=0.75,loctype='world',text='Lowest disturbance',just=0,char_size=1)
graph2.add_drawing_object(DrawText,x=0.5,y=0.34,loctype='world',text='Highest disturbance',just=0,char_size=1)

graph2.legend.configure(char_size=1,loc=(0.05,0.685),loctype='world',box_linestyle=0,fill_pattern=0)

grace.write_file('../../manuscript/figures/persistence_motif_profiles.eps')


