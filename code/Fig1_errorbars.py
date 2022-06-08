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
colors.add_color(94,79,162,'Rpurple')
colors.add_color(50,136,189,'Rblue')
colors.add_color(102,194,165,'Rgreen1')
colors.add_color(171,221,164,'Rgreen2')
colors.add_color(230,245,152,'Rgreen3')
colors.add_color(255,255,191,'Ryellow1')
colors.add_color(254,224,139,'Ryellow2')
colors.add_color(253,174,97,'Rorange1')
colors.add_color(244,109,67,'Rorange2')
colors.add_color(213,62,79,'Rred1')
colors.add_color(158,1,66,'Rred2')
colors.add_color(204,199,242,'Gpurple')
colors.add_color(210,246,249,'Gblue')
colors.add_color(192,244,235,'Ggreen')
colors.add_color(211,250,204,'Ggreen2')
colors.add_color(240,255,222,'Ggreen3')
colors.add_color(245,245,221,'Gyellow')


# colors.add_color(120,120,120,'grey')
# colors.add_color(255,125,125,'lightish_red')
# colors.add_color(200,200,200,'lightgrey')

motnames={'apparent':'Apparent comp.','direct':'Diect comp.','chain':'Three-sp. chain','omnivory':'Omnivory'}

scaledict={'dist':(0.3,0.1366261),'apparent':(0.4051135,0.1461447),'direct':(0.1923225,0.1143021),
  'chain':(0.2489432,0.1272357),'omnivory':(0.1536209,0.1007478)}


def read_lmfile(lmfile):
  lmdict={'(Intercept)':[],'motif':[],'scale(Disturbance)':[],'inter':[]} # estimate, SD

  f=open(lmfile,'r')
  for line in f:
    if line.split()[0]!='"Estimate"':
      pred=line.split()[0].split('"')[1]
      if pred!='(Intercept)' and 'Disturbance' not in pred:
        pred='motif'
      if ':' in pred:
        pred='inter'
      est=float(line.split()[1])
      se=float(line.split()[2])
      try:
        lmdict[pred]=((est,1.96*se))
      except:
        print pred
  f.close()
  return lmdict

def format_graph(graph,motif):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  graph.world.xmin=0
  graph.world.xmax=1
  graph.xaxis.tick.major=.25
  graph.xaxis.ticklabel.configure(format='decimal',prec=2,char_size=.5)

  graph.world.ymin=0.3
  graph.world.ymax=.9
  graph.yaxis.tick.major=.1
  graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=.5)
  graph.yaxis.label.text='Probability of persistence'

  graph.yaxis.label.configure(char_size=.75,just=2,place='normal')
  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.4,place='normal',major_linewidth=.5,minor_linewidth=1)
  graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.4,minor_size=.5,place='both',major_linewidth=.5,minor_linewidth=1)


  upbox=graph.add_dataset([(0,0.9),(1,0.9)])
  upbox.symbol.shape=0
  upbox.line.linestyle=0
  upbox.fill.configure(type=2,color=5)

  dbox=graph.add_dataset([(0,0.85),(1,0.85)])
  dbox.symbol.shape=0
  dbox.line.linestyle=0
  dbox.fill.configure(type=2,color=0)

  return graph

def populate_graph(graph,motif,lmdict):
  j=12
  for basal_p in [0.1,0.18,0.26,0.34,0.42,0.5]:
    p=(basal_p-scaledict['dist'][0])/scaledict['dist'][1]

    dats=[]
    lower=[]
    upper=[]
    if motif!='omnivory':
      xmax=100
    else:
      xmax=55
    for x in range(0,xmax):
      skex=(float(x)/100-scaledict[motif][0])/scaledict[motif][1]
      base=lmdict['(Intercept)'][0]+p*lmdict['scale(Disturbance)'][0]
      maxbase=lmdict['(Intercept)'][0]+lmdict['(Intercept)'][1]+p*(lmdict['scale(Disturbance)'][0]+lmdict['scale(Disturbance)'][1])
      minbase=lmdict['(Intercept)'][0]-lmdict['(Intercept)'][1]+p*(lmdict['scale(Disturbance)'][0]+lmdict['scale(Disturbance)'][1])

      mod=skex*lmdict['motif'][0]+skex*p*lmdict['inter'][0]
      maxmod=skex*(lmdict['motif'][0]+lmdict['motif'][1])+skex*p*(lmdict['inter'][0]+lmdict['inter'][1])
      minmod=skex*(lmdict['motif'][0]-lmdict['motif'][1])+skex*p*(lmdict['inter'][0]-lmdict['inter'][1])

      y=base+mod
      maxy=maxbase+maxmod
      miny=minbase+minmod

      dats.append((float(x)/100,y))
      upper.append((float(x)/100,maxy))
      lower.append((float(x)/100,miny))

    uperr=graph.add_dataset(upper)
    uperr.symbol.shape=0
    uperr.line.linestyle=0
    uperr.fill.configure(type=2,color=j+11,pattern=1)

    derr=graph.add_dataset(lower)
    derr.symbol.shape=0
    derr.line.linestyle=0
    derr.fill.configure(type=2,color=0)


    data=graph.add_dataset(dats)
    data.symbol.shape=0
    data.line.configure(linestyle=1,linewidth=1,color=j)

    j+=1

    if motif=='omnivory':
      data.legend=str(basal_p)
  if motif=='omnivory':
    graph.legend.configure(loc=(.95,.6),loctype='view',box_linestyle=0,box_fill=0,char_size=.5)
    graph.add_drawing_object(DrawText,text='Basal species',char_size=.5,x=.95,y=.628,just=0)
    graph.add_drawing_object(DrawText,text='extinction',char_size=.5,x=.95,y=.615,just=0)
    graph.add_drawing_object(DrawText,text='probability',char_size=.5,x=.95,y=.602,just=0)
  return graph


###############################################################################################
###############################################################################################
#
#           Assemble the plots
#
###############################################################################################
###############################################################################################

# # Far too many species to see anything. Going to apply stats.

grace=MultiPanelGrace(colors=colors)
grace.add_label_scheme('dummy',['Apparent comp.','Three-sp. chain','Direct comp.','Omnivory'])
grace.set_label_scheme('dummy')
for motif in ['apparent','chain','direct','omnivory']:
  lmfile='stat_analysis/'+motif+'_lm.tsv'
  lmdict=read_lmfile(lmfile)
  graph=grace.add_graph(Panel)
  graph=format_graph(graph,motif)
  graph=populate_graph(graph,motif,lmdict)
  graph.panel_label.configure(placement='iuc',char_size=.75,dx=.0,dy=.02)
      
    
grace.multi(rows=1,cols=4,vgap=0,hgap=.01)
grace.hide_redundant_labels()
grace.graphs[0].set_view(0.1, 0.15, 0.28, 0.85)
grace.graphs[1].set_view(0.32, 0.15, 0.5, 0.85)
grace.graphs[2].set_view(0.54, 0.15, 0.72, 0.85)
grace.graphs[3].set_view(0.76, 0.15, 0.94, 0.85)
grace.set_row_xaxislabel(colspan=(None,None),row=0,label='Proportion of role',char_size=.75,perpendicular_offset=.05)

grace.write_file('../manuscript/figures/persistence_vs_motifs.eps')

