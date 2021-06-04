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

# colors.add_color(120,120,120,'grey')
# colors.add_color(255,125,125,'lightish_red')
# colors.add_color(200,200,200,'lightgrey')

# Each point will be mean of one network
def read_datafile(datafile):
  speciesprops={}
  f=open(datafile,'r')
  for line in f:
    if line.split()[0]!='Size':
      S=int(line.split()[0])
      C=float(line.split()[1])
      basal_p=line.split()[2]
      network=line.split()[3]
      species=line.split()[4]
      ID=str(S)+'_'+str(C)+'_'+network+'_'+species
      pers=line.split()[5]
      in_deg=float(line.split()[6])
      out_deg=float(line.split()[7])
      STL=int(line.split()[8])
      PATL=float(line.split()[9])
      chain=int(line.split()[10])
      dircomp=int(line.split()[11])
      omnivory=int(line.split()[12])
      appcomp=int(line.split()[13])
      total=chain+dircomp+omnivory+appcomp
      speciesprops[ID]={'TL':STL,'Deg':in_deg,
      'Raw':{'chain':chain,'dircomp':dircomp,'omnivory':omnivory,'appcomp':appcomp},
      'Prop':{'chain':float(chain)/float(total),'dircomp':float(dircomp)/float(total),
         'omnivory':float(omnivory)/float(total),'appcomp':float(appcomp)/float(total)}}

  f.close()

  return speciesprops

def read_lmfile(lmfile):
  lmdict={'Deg':{'Count':{},'Prop':{}},'TL':{'Count':{},'Prop':{}}}

  f=open(lmfile,'r')
  for line in f:
    if line.split()[0]!='"Predctor"':
      pred=line.split()[1][1:-1]
      role=line.split()[2][1:-1]
      motif=line.split()[3][1:-1]
      intercept=float(line.split()[4][1:-1])
      intercept_SD=float(line.split()[5][1:-1])
      slope=float(line.split()[6][1:-1])
      slope_SD=float(line.split()[7][1:-1])

      lmdict[pred][role][motif]=((intercept,intercept_SD,slope,slope_SD))
  f.close()

  return lmdict

def read_persfiles(Degfile,TLfile):
  persdict={'Deg':{},'TL':{}}
  f=open(Degfile,'r')
  for line in f:
    if line.split()[0]!='"Estimate"':
      persdict['Deg'][line.split()[0][1:-1]]=((float(line.split()[1]),float(line.split()[2])))
  f.close()

  f=open(TLfile,'r')
  for line in f:
    if line.split()[0]!='"Estimate"':
      persdict['TL'][line.split()[0][1:-1]]=((float(line.split()[1]),float(line.split()[2])))
  f.close()

  return persdict

def format_graph(graph,simple,roletype):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  graph.world.xmin=0
  graph.world.ymin=0

  if simple=='Deg':
    graph.world.xmax=100
    graph.xaxis.tick.major=25
    graph.xaxis.label.text='In-degree (number of prey)'
  elif simple=='TL':
    graph.world.xmax=5
    graph.xaxis.tick.major=1
    graph.xaxis.label.text='Trophic level (STL)'

  if roletype=='Prop':
    graph.world.ymax=.600000001
    graph.yaxis.tick.major=.2
    graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=.75)
    graph.yaxis.label.text='Proportion of motif role'
  elif roletype=='Count':
    graph.world.ymax=600
    graph.yaxis.tick.major=200
    graph.yaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)
    graph.yaxis.label.text='Count of motif'
  elif roletype=='persistence':
    graph.world.ymax=1
    graph.yaxis.tick.major=.2
    graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=.75)
    graph.yaxis.label.text='Persistence'


  graph.xaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)
  graph.xaxis.label.configure(char_size=1,just=2,place='normal')

  graph.yaxis.label.configure(char_size=1,just=2,place='normal')
  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.4,place='both',major_linewidth=.5,minor_linewidth=1)
  graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.4,minor_size=.5,place='both',major_linewidth=.5,minor_linewidth=1)

  return graph

def populate_graph(graph,minidict,simple,roletype):
  j=19
  for motif in ['Apparent','Chain','Omnivory','Direct']:
    if j==23:
      j=1
    dats=[]
    lower=[]
    upper=[]
    for x in range(0,100):
      dats.append((x,minidict[motif][0]+x*minidict[motif][2]))
      lower.append((x,minidict[motif][0]-minidict[motif][1]+x*(minidict[motif][2]-minidict[motif][3])))
      upper.append((x,minidict[motif][0]+minidict[motif][1]+x*(minidict[motif][2]+minidict[motif][3])))

    upply=graph.add_dataset(upper)
    upply.symbol.shape=0
    upply.line.configure(linestyle=1,linewidth=1,color=j)
    # upply.fill.configure(pattern=4,type=2,color=j-1,)

    lowly=graph.add_dataset(lower)
    lowly.symbol.shape=0
    lowly.line.configure(linestyle=1,linewidth=1,color=j)
    # lowly.fill.configure(type=2,pattern=0,color=0)

    data=graph.add_dataset(dats)
    data.symbol.shape=0
    data.line.configure(linestyle=1,linewdith=1.5,color=j)

    if simple=='TL':
      data.legend=motif
    j+=1

  if simple=='TL':
    graph.add_drawing_object(DrawText,text='Motif',x=5.25,y=0.47,loctype='world',char_size=.75) 
    graph.legend.configure(char_size=.75,loc=(5.25,0.450),loctype='world',box_linestyle=0,fill_pattern=0)

  return graph

def populate_persgraph(graph,persdict,simple):
  if simple=='Deg':
    key='in_Degree'
  else:
    key='STL'

  j=12
  for basal_p in [0,0.2,0.4,0.6,0.8,1]:
    p=0.1+0.4*basal_p

    dats=[]
    # lower=[]
    # upper=[]
    for x in range(0,100):
      xcomp=persdict['(Intercept)'][0]+x*persdict[key][0]
      bcomp=basal_p*persdict['Basal_p'][0]+x*basal_p*persdict[key+':Basal_p'][0]
      dats.append((x,xcomp+bcomp))

    data=graph.add_dataset(dats)
    data.symbol.shape=0
    data.line.configure(linestyle=1,linewdith=1.5,color=j)

    j+=1
    if simple=='TL':
      data.legend=str(p)

  if simple=='TL':
    graph.add_drawing_object(DrawText,text='Basal species',x=5.25,y=0.58,loctype='world',char_size=.75) 
    graph.add_drawing_object(DrawText,text='extinction',x=5.25,y=0.52,loctype='world',char_size=.75) 
    graph.add_drawing_object(DrawText,text='probability',x=5.25,y=0.47,loctype='world',char_size=.75) 
    graph.legend.configure(char_size=.75,loc=(5.25,0.45),loctype='world',box_linestyle=0,fill_pattern=0)


  return graph

###############################################################################################
###############################################################################################
#
#           Assemble the plots
#
###############################################################################################
###############################################################################################

# # Far too many species to see anything. Going to apply stats.
# datafile='../data/3sp_roles_participation.tsv'
# speciesprops=read_datafile(datafile)

lmfile='stat_analysis/roles_vs_TL_Deg.tsv'
lmdict=read_lmfile(lmfile)
persdict=read_persfiles('stat_analysis/persistence_vs_Deg.tsv','stat_analysis/persistence_vs_TL.tsv')

grace=MultiPanelGrace(colors=colors)
for graphtype in ['Persistence','Proportion']:
  for simple in ['Deg','TL']:
    if graphtype=='Proportion':
      graph=grace.add_graph(Panel)
      graph=format_graph(graph,simple,'Prop')
      graph=populate_graph(graph,lmdict[simple]['Prop'],simple,'Prop')
      graph.panel_label.configure(placement='iul',char_size=1,dx=.03,dy=.03)
    elif graphtype=='Persistence':
      graph=grace.add_graph(Panel)
      graph=format_graph(graph,simple,'persistence')
      graph=populate_persgraph(graph,persdict[simple],simple)
      graph.panel_label.configure(placement='iul',char_size=1,dx=.03,dy=.03)
      
    
grace.multi(rows=2,cols=2,vgap=.05,hgap=.05)
grace.hide_redundant_labels()
# grace.set_row_xaxislabel(label='in-degree (number of prey)',row=0,colspan=(None,None),char_size=1,perpendicular_offset=.05)
# grace.set_row_xaxislabel(label='Trophic level (STL)',row=1,colspan=(None,None),char_size=1,perpendicular_offset=.05)
# grace.set_col_yaxislabel(label='Count of motif',col=0,rowspan=(None,None),char_size=1,perpendicular_offset=.07)
# grace.set_col_yaxislabel(label='Proportion of motif role',col=1,rowspan=(None,None),char_size=1,perpendicular_offset=.07)

grace.write_file('../manuscript/figures/roles_vs_TL.eps')

