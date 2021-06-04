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

def format_graph(graph,simple):
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
    graph.xaxis.label.text='Proportion of network profile'

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
    if nprey==0:
      linestyle=2
    elif nprey==1:
      linestyle=3
    else:
      linestyle=1

    if npred==0:
      shap=3
    elif npred==1:
      shap=1
    else:
      shap=2

    # betas go intercept, disturbance, position, interaction
    lms=datadict['betas'][pos]
    for dist in [0.1,0.26,0.42]:
      dats=[]

      for x in [0,0.2,0.4,0.6,0.8]:
        y=lms[0]+dist*lms[1]+x*lms[2]+x*dist*lms[3]

        dats.append((x,y))

      data=graph.add_dataset(dats)
      data.symbol.configure(shape=shap,color=j,fill_color=j)
      data.line.configure(linestyle=linestyle,linewdith=1.5,color=j)

      if motif=='Direct competition' and nprey==1:
        data.legend=str(dist)        
      j+=2

  if motif=='Three-species chain':
    data=graph.add_dataset([(0,0),(100,0)])
    data.symbol.shape=0
    data.line.configure(linestyle=2,linewidth=1.5,color=1)
    data.legend='0'

    data=graph.add_dataset([(0,0),(100,0)])
    data.symbol.shape=0
    data.line.configure(linestyle=3,linewidth=1.5,color=1)
    data.legend='1'

    data=graph.add_dataset([(0,0),(100,0)])
    data.symbol.shape=0
    data.line.configure(linestyle=1,linewidth=1.5,color=1)
    data.legend='2'

  if motif=='Omnivory':
    data=graph.add_dataset([(0,0),(100,0)])
    data.symbol.configure(shape=3)
    data.legend='0'

    data=graph.add_dataset([(0,0),(100,0)])
    data.symbol.shape=1
    data.legend='1'

    data=graph.add_dataset([(0,0),(100,0)])
    data.symbol.shape=2
    data.legend='2'


  if motif=='Three-species chain':
    graph.add_drawing_object(DrawText,text='No. prey',x=0.85,y=0.695,loctype='world',char_size=.5)
    graph.legend.configure(char_size=.5,loc=(0.85,0.675),loctype='world',box_linestyle=0,fill_pattern=0)
  elif motif=='Omnivory':
    graph.add_drawing_object(DrawText,text='No. preds',x=1.75,y=0.295,loctype='world',char_size=.5)
    graph.legend.configure(char_size=.5,loc=(1.75,0.275),loctype='world',box_linestyle=0,fill_pattern=0)
  elif motif=='Direct competition':
    graph.legend.configure(char_size=.5,loc=(0.85,0.675),loctype='world',box_linestyle=0,fill_pattern=0)
    graph.add_drawing_object(DrawText,text='Basal species',x=0.85,y=0.755,loctype='world',char_size=.5)
    graph.add_drawing_object(DrawText,text='extinction',x=0.85,y=0.725,loctype='world',char_size=.5)
    graph.add_drawing_object(DrawText,text='probability',x=0.85,y=0.695,loctype='world',char_size=.5)
  return graph

def populate_graph(graph,nprey,datadict):
  motifs=datadict['betas'].keys()
  positions=[key for key in motifs if int(key.split('..')[2].split('.')[0])==nprey]

  for pos in positions:
    motif=[c for c in codes if pos.split('..')[0]==codes[c]][0]
    npred=int(pos.split('..')[1])
    nprey=int(pos.split('..')[2].split('.')[0])
    if nprey==0:
      linestyle=2
    elif nprey==1:
      linestyle=3
    else:
      linestyle=1

    if npred==0:
      shap=3
    elif npred==1:
      shap=1
    else:
      shap=2

    if motif=='Omnivory':
      j=13
    elif motif=='Three-species chain':
      j=15
    elif motif=='Apparent competition':
      j=16
    else:
      j=17

    # betas go intercept, disturbance, position, interaction
    lms=datadict['betas'][pos]
    for dist in [0.1,0.26,0.5]:
      dats=[]

      for x in [0,0.2,0.4,0.6,0.8]:
        y=lms[0]+dist*lms[1]+x*lms[2]+x*dist*lms[3]

        dats.append((x,y))

      data=graph.add_dataset(dats)
      data.symbol.configure(shape=shap,color=j,fill_color=j)
      data.line.configure(linestyle=linestyle,linewdith=1.5,color=j)

      if nprey==0 and dist==0.1:
        data.legend=motif        

  if nprey==2:
    data=graph.add_dataset([(0,0),(100,0)])
    data.symbol.configure(shape=3)
    data.legend='0'

    data=graph.add_dataset([(0,0),(100,0)])
    data.symbol.shape=1
    data.legend='1'

    data=graph.add_dataset([(0,0),(100,0)])
    data.symbol.shape=2
    data.legend='2'


  if nprey==0:
    # graph.add_drawing_object(DrawText,text='Motif',x=0.1,y=0.3,loctype='world',char_size=.5)
    graph.legend.configure(char_size=.5,loc=(0.45,0.35),loctype='world',box_linestyle=0,fill_pattern=0)
  elif nprey==2:
    graph.add_drawing_object(DrawText,text='No. preds',x=.95,y=0.795,loctype='world',char_size=.5)
    graph.legend.configure(char_size=.5,loc=(.95,0.75),loctype='world',box_linestyle=0,fill_pattern=0)
  # elif motif=='Direct competition':
  #   graph.legend.configure(char_size=.5,loc=(0.85,0.675),loctype='world',box_linestyle=0,fill_pattern=0)
  #   graph.add_drawing_object(DrawText,text='Basal species',x=0.85,y=0.755,loctype='world',char_size=.5)
  #   graph.add_drawing_object(DrawText,text='extinction',x=0.85,y=0.725,loctype='world',char_size=.5)
  #   graph.add_drawing_object(DrawText,text='probability',x=0.85,y=0.695,loctype='world',char_size=.5)

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

grace=MultiPanelGrace(colors=colors)
grace.add_label_scheme('dummy',['Omnivory','Three-species chain','Apparent competition','Direct competition'])
grace.set_label_scheme('dummy')

# Organized by motif
for mot in ['Omnivory','Three-species chain','Apparent competition','Direct competition']:
  graph2=grace.add_graph(Panel)
  graph2=format_graph(graph2,'persistence')
  graph2=populate_persgraph(graph2,mot,predictors)
  graph2.panel_label.configure(placement='ouc',char_size=1,dx=.0,dy=.01)

grace.multi(rows=2,cols=2,vgap=.07,hgap=.07)
grace.hide_redundant_labels()

grace.write_file('../../manuscript/figures/persistence_positions_bymotif.eps')

grace=MultiPanelGrace(colors=colors)
grace.add_label_scheme('dummy',['0 prey','1 prey','2 prey',''])
grace.set_label_scheme('dummy')

# Organized by no. prey
for prey in [0,1,2]:
  graph2=grace.add_graph(Panel)
  graph2=format_graph(graph2,'persistence')
  graph2=populate_graph(graph2,prey,predictors)
  graph2.panel_label.configure(placement='ouc',char_size=1,dx=.0,dy=.01)

graph=grace.add_graph(Panel)
graph=format_graph(graph,'dummy')

grace.multi(rows=2,cols=2,vgap=.07,hgap=.07)
grace.hide_redundant_labels()

grace.write_file('../../manuscript/figures/persistence_positions_byprey.eps')

