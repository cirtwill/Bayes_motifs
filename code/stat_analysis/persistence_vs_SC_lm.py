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

colors=ColorBrewerScheme('Spectral',n=13)  # The blue is very beautiful but maybe harder to see.
colors.add_color(173,132,191,'Purple')
colors.add_color(34,90,34,'Green')
# colors.add_color(120,120,120,'grey')
# colors.add_color(255,125,125,'lightish_red')
# colors.add_color(200,200,200,'lightgrey')


scales={
  'Size':(79.1882,16.62225),
  'Connectance':(0.1044983,0.05559927),
  'Disturbance':(0.3,0.1366261),
}

lms={
  'Intercept': 5.723e-01,
  'Size':-5.148e-04,
  'Connectance':-1.293e-02,
  'Disturbance':-1.393e-01,
  'Size:Connectance':-5.245e-04,
  'Size:Disturbance':-3.485e-03,
  'Connectance:Disturbance':-6.998e-05,
  'Size:Connectance:Disturbance':-2.693e-04 }

def format_graph(graph,simple):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1


  graph.world.ymin=0
  graph.world.ymax=1
  graph.yaxis.tick.major=.2
  graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=.75)
  graph.yaxis.label.text='Mean persistence'
  if simple=='Size':
    graph.xaxis.label.text='Network size'
    graph.world.xmin=45
    graph.world.xmax=105
    graph.xaxis.tick.major=10
    graph.xaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)
  else:
    graph.xaxis.label.text='Connectance'
    graph.world.xmin=0
    graph.world.xmax=.2
    graph.xaxis.tick.major=.05
    graph.xaxis.ticklabel.configure(format='decimal',prec=2,char_size=.75)

  graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=.75)
  graph.xaxis.label.configure(char_size=1,just=2,place='normal')

  graph.yaxis.label.configure(char_size=1,just=2,place='normal')
  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.4,place='both',major_linewidth=.5,minor_linewidth=1)
  graph.yaxis.tick.configure(onoff='on',minor_ticks=1,major_size=.4,minor_size=.5,place='both',major_linewidth=.5,minor_linewidth=1)

  return graph

def populate_persgraph(graph,simple):

  if simple=='Size':
    levels=[50,60,70,80,90,100]
    altlevels=[0.02,0.1,0.18]
    altpred='Connectance'
  else:
    levels=[0.02,0.06,0.1,0.14,0.18]
    altlevels=[50,75,100]
    altpred='Size'


  for Disturbance in [0.1,0.3,0.5]:
    scalbp=(Disturbance-scales['Disturbance'][0])/scales['Disturbance'][1]
    base=lms['Intercept']+scalbp*lms['Disturbance']
    for al in altlevels:
      dats=[]
      scalalt=(al-scales[altpred][0])/scales[altpred][1]
      base2=base+scalalt*lms[altpred]+scalbp*scalalt*lms[altpred+':Disturbance']
      for level in levels:
        scalmain=(level-scales[simple][0])/scales[simple][1]
        base3=base2+scalmain*lms[simple]+scalbp*scalmain*lms[simple+':Disturbance']
        final=base3+scalmain*scalalt*lms['Size:Connectance']+scalmain*scalalt*scalbp*lms['Size:Connectance:Disturbance']

        dats.append((level,final))

      data=graph.add_dataset(dats)
      if al==altlevels[0]:
        sty=2
      elif al==altlevels[1]:
        sty=1
      else:
        sty=5
      if Disturbance==0.1:
        col=2
      elif Disturbance==0.5:
        col=14
      else:
        col=6

        print al, col
      data.symbol.shape=0
      data.line.configure(linestyle=sty,linewdith=1.5,color=col)

      if simple=='Size' and al==altlevels[1]:
        data.legend='p(Basal extinct)='+str(Disturbance)
      # Doesn't work - add DrawText.
      # if simple=='Size' and Basal_p==0.5:
      #   data.legend='C='+str(al)        

  if simple=='Size':
    # graph.add_drawing_object(DrawText,text='p(Basal)',x=107,y=0.825,loctype='world',char_size=.5)
    graph.legend.configure(char_size=.5,loc=(50,0.3),loctype='world',box_linestyle=0,fill_pattern=0)
    # graph.add_drawing_object(DrawText,text='Connectance',x=107,y=0.5,loctype='world',char_size=.5)

    # dummy1=graph.add_dataset([(10000,100),(1000000,40)])
    # dummy1.symbol.shape=0
    # dummy1.line.configure(linestyle=2,linewidth=1.5)
    # dummy1.legend='C=0.02'
    # dummy2=graph.add_dataset([(10000,100),(1000000,40)])
    # dummy2.symbol.shape=0
    # dummy2.line.configure(linestyle=1,linewidth=1.5)
    # dummy2.legend='C=0.1'
    # dummy3=graph.add_dataset([(10000,100),(1000000,40)])
    # dummy3.symbol.shape=0
    # dummy3.line.configure(linestyle=3,linewidth=1.5)
    # dummy3.legend='C=0.18'

  return graph

###############################################################################################
###############################################################################################
#
#           Assemble the plots
#
###############################################################################################
###############################################################################################

grace=MultiPanelGrace(colors=colors)

for simple in ['Size','Connectance']:
  graph2=grace.add_graph(Panel)
  graph2=format_graph(graph2,simple)
  graph2=populate_persgraph(graph2,simple)
  graph2.panel_label.configure(placement='iul',char_size=1,dx=.03,dy=.03)
    

grace.multi(rows=1,cols=2,vgap=.09,hgap=.04)
grace.hide_redundant_labels()
# grace.set_row_xaxislabel(label='in-degree (number of prey)',row=0,colspan=(None,None),char_size=1,perpendicular_offset=.05)
# grace.set_row_xaxislabel(label='Trophic level (STL)',row=1,colspan=(None,None),char_size=1,perpendicular_offset=.05)
# grace.set_col_yaxislabel(label='Count of motif',col=0,rowspan=(None,None),char_size=1,perpendicular_offset=.07)
# grace.set_col_yaxislabel(label='Proportion of motif role',col=1,rowspan=(None,None),char_size=1,perpendicular_offset=.07)

grace.write_file('../../manuscript/figures/persistence_vs_SC_lm.eps')


