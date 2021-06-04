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

mot_ranges={
  'omnivory':(0.005952381,0.345547537),
  'apparent':(0.2499348,0.7647059),
  'direct':(0.05010438,0.35878963),
  'chain':(0.1088861,0.4125000),
}

scales={ # center, scale
  'omnivory':(0.169822,0.07561787),
  'apparent':(0.4169211,0.06134048),
  'direct':(0.1805743,0.04025804),
  'chain':(0.2326826,0.04136131),
  'Disturbance':(0.3,0.1366261)
}

lms={ # Intercept, motif, disturbance, interaction
  'omnivory':(5.723e-01,-1.328e-02,-1.393e-01,-4.886e-04),
  'apparent':(5.723e-01,1.053e-02,-1.393e-01,4.781e-04),
  'direct':(5.723e-01,2.618e-03,-1.393e-01,3.844e-04),
  'chain':(0.5723309,0.0061244,-0.1393351,-0.0001900)  
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

def binner(netprops,binby):
  if binby=='omnivory':
    bindict={'chain':[],'dircomp':[],'appcomp':[]}
    step=.01
    minni=0
    for i in range(0,100):
      binmin=minni+step*i
      binmed=minni+step*i+step/2
      binmax=minni+step*(i+1)

      for motif in bindict:
        motpoints=[]
        for ID in netprops:
          if netprops[ID]['omnivory']>=binmin and netprops[ID]['omnivory']<binmax:
            motpoints.append(netprops[ID][motif])
        if motpoints!=[]:
          bindict[motif].append((binmed,np.mean(motpoints),np.std(motpoints)/len(motpoints)))
  else:
    bindict={}

  return bindict

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


  return graph

def populate_graph(graph,bindict):
  j=2
  motnames={'appcomp':'Apparent competition','chain':'3-species chain','dircomp':'Direct competition'}
  for motif in ['appcomp','chain','dircomp']:
    data=graph.add_dataset(bindict[motif],type='xydy')
    data.symbol.shape=0
      
    data.line.configure(linestyle=1,linewdith=1.5,color=j)
    data.errorbar.configure(riser_color=j,color=j,riser_linewidth=1,size=.5,linewidth=1)

    # data.legend=motnames[motif]
    # graph.legend.configure(char_size=.75,loc=(0.2,0.58),loctype='world',box_linestyle=0,fill_pattern=0)
    if j==11:
      j=13
    elif j==2:
      j=11

  return graph

def populate_persgraph(graph,motif):

  j=12
  for dist in [0.1,0.18,0.26,0.34,0.42,0.5]:

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
    data.line.configure(linestyle=1,linewdith=1.5,color=j)


    if motif=='omnivory':
      data.legend=str(dist)
    j+=1

  if motif=='omnivory':
    graph.add_drawing_object(DrawText,text='Basal species',x=0.5,y=0.755,loctype='world',char_size=.5)
    graph.add_drawing_object(DrawText,text='extinction',x=0.5,y=0.725,loctype='world',char_size=.5)
    graph.add_drawing_object(DrawText,text='probability',x=0.5,y=0.695,loctype='world',char_size=.5)
    graph.legend.configure(char_size=.5,loc=(0.5,0.685),loctype='world',box_linestyle=0,fill_pattern=0)
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

grace=MultiPanelGrace(colors=colors)
grace.add_label_scheme('dummy',['Omnivory','Three-species chain','Apparent competition','Direct competition'])
grace.set_label_scheme('dummy')

for mot in ['omnivory','chain','apparent','direct']:
  graph2=grace.add_graph(Panel)
  graph2=format_graph(graph2,'persistence')
  graph2=populate_persgraph(graph2,mot)
  graph2.panel_label.configure(placement='ouc',char_size=1,dx=.0,dy=.01)
    
# graph=grace.add_graph(Panel)
# graph=format_graph(graph,'profile')
# omnibin=binner(netprops,'omnivory')
# graph=populate_graph(graph,omnibin)
# graph.panel_label.configure(placement='iul',char_size=1,dx=.03,dy=.03)

grace.multi(rows=2,cols=2,vgap=.07,hgap=.07)
grace.hide_redundant_labels()
# grace.set_row_xaxislabel(label='in-degree (number of prey)',row=0,colspan=(None,None),char_size=1,perpendicular_offset=.05)
# grace.set_row_xaxislabel(label='Trophic level (STL)',row=1,colspan=(None,None),char_size=1,perpendicular_offset=.05)
# grace.set_col_yaxislabel(label='Count of motif',col=0,rowspan=(None,None),char_size=1,perpendicular_offset=.07)
# grace.set_col_yaxislabel(label='Proportion of motif role',col=1,rowspan=(None,None),char_size=1,perpendicular_offset=.07)

grace.write_file('../../manuscript/figures/persistence_motif_profiles.eps')

