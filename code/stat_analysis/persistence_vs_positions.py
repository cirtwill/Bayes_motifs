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
colors.add_color(68,1,84,'Rpurple')
colors.add_color(65,68,135,'Rblue')
colors.add_color(42,120,142,'Rgreen1')
colors.add_color(34,168,132,'Rgreen2')
colors.add_color(122,209,81,'Rgreen3')
colors.add_color(253,231,37,'Ryellow')
colors.add_color(254,224,139,'Ryellow2')
colors.add_color(253,174,97,'Rorange1')
colors.add_color(244,109,67,'Rorange2')
colors.add_color(213,62,79,'Rred1')
colors.add_color(158,1,66,'Rred2')

codes={'Three-species chain':'p.12','Omnivory':'p.38','Direct competition':'p.36','Apparent competition':'p.6'}

# Each point will be mean of one network
def read_datafile(datafile):
  netprops={'betas':{},'pvals':{}}
  f=open(datafile,'r')
  for line in f:
    if line.split()[0]!='"Position"':
      position=line.split()[0][1:-1]
      betas=[float(b[1:-1]) for b in line.split()[1:5]]
      pvals=[float(p[1:-1]) for p in line.split()[5:]]
      netprops['betas'][position]=betas
      netprops['pvals'][position]=pvals
  f.close()
  return netprops

def format_graph(graph,simple,motif):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  graph.world.xmin=0
  graph.world.ymin=0

  if simple=='persistence':
    graph.world.ymin=.1
    graph.world.ymax=.9
    graph.world.xmax=.8
    graph.xaxis.tick.major=.2
    graph.yaxis.tick.major=.1
    graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=.75)
    graph.yaxis.label.text='Mean persistence'
    graph.xaxis.label.text='Proportion of '+motif+' in network profile'

    graph.xaxis.ticklabel.configure(format='decimal',prec=1,char_size=.75)
    graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=.75)
    graph.xaxis.label.configure(char_size=1,just=2,place='normal')

    graph.yaxis.label.configure(char_size=1,just=2,place='normal')
    graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.4,place='both',major_linewidth=.5,minor_linewidth=1)
    graph.yaxis.tick.configure(onoff='on',minor_ticks=1,major_size=.4,minor_size=.5,place='both',major_linewidth=.5,minor_linewidth=1)

  elif simple=='dummy':
    graph.yaxis.bar.configure(linewidth=0,linestyle=0,color=0)
    graph.xaxis.bar.configure(linewidth=0,linestyle=0,color=0)
    graph.frame.configure(linewidth=0,linestyle=0,color=0)
    graph.xaxis.tick.onoff='off'
    graph.yaxis.tick.onoff='off'
    graph.yaxis.label.char_size=0
    graph.xaxis.label.char_size=0
    graph.xaxis.ticklabel.char_size=0
    graph.yaxis.ticklabel.char_size=0

  return graph

def populate_persgraph(graph,motif,datadict):
  positions=[key for key in datadict['betas'] if codes[motif] in key]

  for pos in positions:
    npred=int(pos.split('..')[1])
    nprey=int(pos.split('..')[2].split('.')[0])
    j=12
    if motif=='Omnivory':
      if nprey==0:
        linestyle=2
      elif nprey==1:
        linestyle=3
      else:
        linestyle=1
    elif motif=='Three-species chain':
      if nprey==0:
        linestyle=2
      if npred==0:
        linestyle=1
      if nprey==1 and npred==1:
        linestyle=3
    else: # both competition motifs can use same code
      if npred==0:
        linestyle=1
      else:
        linestyle=2


    # betas go intercept, disturbance, position, interaction
    lms=datadict['betas'][pos]
    for dist in [0.1,0.26,0.42]:
      dats=[]

      for x in [0,0.2,0.4,0.6,0.8]:
        y=lms[0]+dist*lms[1]+x*lms[2]+x*dist*lms[3]

        dats.append((x,y))

      data=graph.add_dataset(dats)
      data.symbol.configure(shape=0,color=j,fill_color=j)
      data.line.configure(linestyle=linestyle,linewdith=1.5,color=j)

      if linestyle==1:
        data.legend=str(dist)        
      j+=2

  if motif=='Three-species chain':
    data=graph.add_dataset([(0,0),(100,0)])
    data.symbol.configure(shape=0)
    data.legend='Top'

    data=graph.add_dataset([(0,0),(100,0)])
    data.symbol.shape=0
    data.legend='Middle'
    data.line.linestyle=3

    data=graph.add_dataset([(0,0),(100,0)])
    data.symbol.shape=0
    data.legend='Bottom'
    data.line.linestyle=2

  elif motif=='Omnivory':
    data=graph.add_dataset([(0,0),(100,0)])
    data.symbol.configure(shape=0)
    data.legend='Top'

    data=graph.add_dataset([(0,0),(100,0)])
    data.symbol.shape=0
    data.legend='Middle'
    data.line.linestyle=3

    data=graph.add_dataset([(0,0),(100,0)])
    data.symbol.shape=0
    data.legend='Bottom'
    data.line.linestyle=2

  else:
    data=graph.add_dataset([(0,0),(100,0)])
    data.symbol.configure(shape=0)
    data.legend='Top'

    data=graph.add_dataset([(0,0),(100,0)])
    data.symbol.shape=0
    data.legend='Bottom'
    data.line.linestyle=2


  graph.legend.configure(box_linestyle=0,char_size=.75,loctype='world',loc=(0.05,0.35))
  return graph

###############################################################################################
###############################################################################################
#
#           Assemble the plots
#
###############################################################################################
###############################################################################################

# # Far too many species to see anything. Going to apply stats.
datafile='persistence_vs_positions.tsv'
predictors=read_datafile(datafile)

for mot in ['Omnivory','Three-species chain','Apparent competition','Direct competition']:

  grace=Grace(colors=colors)
  graph2=grace.add_graph()
  graph2=format_graph(graph2,'persistence',mot)
  graph2=populate_persgraph(graph2,mot,predictors)
  # graph2.panel_label.configure(placement='ouc',char_size=1,dx=.0,dy=.01)

  grace.write_file('../../manuscript/figures/persistence_positions_'+mot.split()[0]+'.eps')

