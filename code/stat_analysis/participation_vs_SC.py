import sys
import os
import math

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
  'omnivory':(0,0.54917),
  'apparent':(0,1),
  'direct':(0,1),
  'chain':(0,1),
}

motif_names={'Omnivory':'Omnivory','Apparent':'Apparent competition','Direct':'Direct competition','Chain':'Three-species chain'}

# Each point will be mean of one network
def read_datafile(datafile):
  netprops={}
  f=open(datafile,'r')
  for line in f:
    if line.split()[0]!='"V1"':
      motif=line.split()[1][1:-1]
      intercept=float(line.split()[2][1:-1])
      pBasal=float(line.split()[3][1:-1])
      C=float(line.split()[4][1:-1])
      interaction=float(line.split()[5][1:-1])

      netprops[motif]=((intercept,pBasal,C,interaction))

  f.close()
  return netprops

def format_graph(graph,simple):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  graph.xaxis.label.configure(char_size=1,just=2,place='normal')
  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.4,place='both',major_linewidth=.5,minor_linewidth=1)
  graph.world.xmin=0
  if simple=='C':
    graph.world.xmax=.2
    graph.xaxis.tick.major=.04
    graph.xaxis.label.text='Connectance'
    graph.xaxis.ticklabel.configure(format='decimal',prec=2,char_size=.75)
  else:
    graph.world.xmin=45
    graph.world.xmax=105
    graph.xaxis.tick.major=10
    graph.xaxis.label.text='Network size'
    graph.xaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)

  graph.world.ymin=0
  graph.world.ymax=.6
  graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=.75)
  graph.yaxis.label.configure(text='Mean persistence',char_size=1,just=2,place='normal')
  graph.yaxis.tick.configure(major=.2,onoff='on',minor_ticks=0,major_size=.4,minor_size=.5,place='both',major_linewidth=.5,minor_linewidth=1)

  return graph

def populate_persgraph(graph,xprop,val,netprops):
  print xprop, val
  for motif in netprops:
    if motif=='Apparent':
      j=12
      sty=2
    elif motif=='Chain':
      j=14
      sty=3
    elif motif=='Direct':
      j=15
      sty=5
    else:
      j=18
      sty=1
    betas=netprops[motif]
    dats=[]
    if xprop=='S':
      for x in [50,60,70,80,90,100]:
        y=betas[0]+x*betas[1]+val*betas[2]+x*val*betas[3]
        logity=math.exp(y)/(1+math.exp(y))
        dats.append((x,logity))
    elif xprop=='C':
      for x in [0.02,0.06,0.1,0.14,0.18]:
        y=betas[0]+x*betas[2]+val*betas[1]+x*val*betas[3]      
        logity=math.exp(y)/(1+math.exp(y))
        dats.append((x,logity))

    data=graph.add_dataset(dats)
    data.symbol.shape=0
    data.line.configure(linestyle=sty,linewdith=1.5,color=j)

    if xprop=='C' and val==70:
      data.legend=motif_names[motif]

  graph.legend.configure(box_linestyle=0,box_fill_pattern=0,char_size=.75,loc=(0.21,0.5),loctype='world')
  return graph

###############################################################################################
###############################################################################################
#
#           Assemble the plots
#
###############################################################################################
###############################################################################################

# # Far too many species to see anything. Going to apply stats.
datafile='roles_vs_SC.tsv'
netprops=read_datafile(datafile)

grace=MultiPanelGrace(colors=colors)
grace.add_label_scheme('dummy',['C=0.1','S=70'])
grace.set_label_scheme('dummy')


for xprop in ['S','C']:
  if xprop=='S':
    val=0.1
  elif xprop=='C':
    val=70
  graph2=grace.add_graph(Panel)
  graph2=format_graph(graph2,xprop)
  graph2=populate_persgraph(graph2,xprop,val,netprops)
  graph2.panel_label.configure(placement='iul',char_size=.75,dx=.02,dy=.02)


grace.multi(rows=1,cols=2,vgap=.06,hgap=.06)
grace.hide_redundant_labels()
grace.set_col_yaxislabel(label='Proportion of motif participation',col=0,rowspan=(None,None),char_size=1,perpendicular_offset=.07)

grace.write_file('../../manuscript/figures/participation_vs_SC.eps')


