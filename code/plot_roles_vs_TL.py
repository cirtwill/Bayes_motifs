import sys
import os
import math
import participation_vs_SC as SC

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

colors=ColorBrewerScheme('Greys',n=10)  
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

colors.add_color(228,26,28,'omnivory')
colors.add_color(55,126,184,'direct') 
colors.add_color(255,127,0,'apparent') 
colors.add_color(152,78,163,'chain') 

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
    graph.xaxis.tick.major=50
    graph.xaxis.label.text='In-degree (number of prey)'
  elif simple=='TL':
    graph.world.xmax=5
    graph.xaxis.tick.major=2
    graph.xaxis.label.text='Trophic level (STL)'

  if roletype=='Prop':
    graph.world.ymax=.600000001
    graph.yaxis.tick.major=.2
    graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=.5)
    graph.yaxis.label.text='Proportion of motif in role'


  graph.xaxis.ticklabel.configure(format='decimal',prec=0,char_size=.5)
  graph.xaxis.label.configure(char_size=.75,just=2,place='normal')

  graph.yaxis.label.configure(char_size=.75,just=2,place='normal')
  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.4,place='both',major_linewidth=.5,minor_linewidth=1)
  graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.4,minor_size=.5,place='both',major_linewidth=.5,minor_linewidth=1)

  return graph

def populate_graph(graph,minidict,simple,roletype):
  for motif in ['Apparent','Chain','Omnivory','Direct']:
    if motif=='Apparent':
      sty=2
    elif motif=='Chain':
      sty=3
    elif motif=='Direct':
      sty=5
    else:
      sty=1
    dats=[]
    lower=[]
    upper=[]
    for x in range(0,100):
      y=minidict[motif][0]+x*minidict[motif][2]
      logity=math.exp(y)/(1+math.exp(y))
      dats.append((x,logity))
      lowy=minidict[motif][0]-minidict[motif][1]+x*(minidict[motif][2]-minidict[motif][3])
      logitlow=math.exp(lowy)/(1+math.exp(lowy))
      lower.append((x,logitlow))
      upy=minidict[motif][0]+minidict[motif][1]+x*(minidict[motif][2]+minidict[motif][3])
      logitup=math.exp(upy)/(1+math.exp(upy))
      upper.append((x,logitup))

    # upply=graph.add_dataset(upper)
    # upply.symbol.shape=0
    # upply.line.configure(linestyle=sty,linewidth=1,color=motif.lower())

    # lowly=graph.add_dataset(lower)
    # lowly.symbol.shape=0
    # lowly.line.configure(linestyle=sty,linewidth=1,color=motif.lower())

    data=graph.add_dataset(dats)
    data.symbol.shape=0
    data.line.configure(linestyle=sty,linewidth=3.5,color=motif.lower())

    if simple=='TL':
      data.legend=motif

  if simple=='TL':
    graph.add_drawing_object(DrawText,text='Motif',x=3.15,y=-.21,loctype='world',char_size=.5) 
    graph.legend.configure(char_size=.5,loc=(3,-.22),loctype='world',box_linestyle=0,fill_pattern=0)

  return graph

###############################################################################################
###############################################################################################
#
#           Assemble the plots
#
###############################################################################################
###############################################################################################

# # Far too many species to see anything. Going to apply stats.
lmfile='stat_analysis/roles_vs_TL_Deg.tsv'
lmdict=read_lmfile(lmfile)
persdict=read_persfiles('stat_analysis/persistence_vs_Deg_norandom.tsv','stat_analysis/persistence_vs_TL.tsv')

grace=MultiPanelGrace(colors=colors)
grace.add_label_scheme('dummy',['A','B','C (S=50)','D (C=0.20)'])
# 'D (C=0.02)','E (S=100)','F (C=0.20)'])
grace.set_label_scheme('dummy')

graphtype='Proportion'
for simple in ['Deg','TL']:
  if graphtype=='Proportion':
    graph=grace.add_graph(Panel)
    graph=format_graph(graph,simple,'Prop')
    graph=populate_graph(graph,lmdict[simple]['Prop'],simple,'Prop')
    graph.panel_label.configure(placement='iul',char_size=.5,dx=.02,dy=.02)

datafile='stat_analysis/roles_vs_SC.tsv'
netprops=SC.read_datafile(datafile)
SEfile='stat_analysis/roles_vs_SC_SE.tsv'
netSE=SC.read_datafile(SEfile)

# for sel in ['lo','high']:
for xprop in ['C','S']:
  if xprop=='S':
    sel='high'
    if sel=='lo':
      val=0.02
    else:
      val=0.2
  elif xprop=='C':
    sel='lo'
    if sel=='lo':
      val=50
    else:
      val=100
  graph2=grace.add_graph(Panel)
  graph2=SC.format_graph(graph2,xprop)
  graph2=SC.populate_persgraph(graph2,xprop,val,netprops,netSE)
  graph2.panel_label.configure(placement='iul',char_size=.5,dx=.02,dy=.02)


grace.multi(rows=2,cols=2,vgap=.07,hgap=.04)
grace.hide_redundant_labels()
grace.set_col_yaxislabel(col=0,rowspan=(None,None),label='Proportion of motif in participation vector',char_size=1,perpendicular_offset=.05)
for graph in grace.graphs:
  print(graph.get_view())
grace.graphs[0].set_view(0.10, 0.55, 0.5, 0.85)
grace.graphs[1].set_view(0.55, 0.55, 0.95, 0.85)
grace.graphs[2].set_view(0.10, 0.17, 0.5, 0.47)
grace.graphs[3].set_view(0.55, 0.17, 0.95, 0.47)

grace.write_file('../manuscript/figures/roles_vs_TL.eps')

grace=MultiPanelGrace(colors=colors)
grace.add_label_scheme('dummy',['A (S=50)','B (S=70)','C (S=100)','D (C=0.02)', 'E (C=0.1)', 'F (C=0.2)'])
grace.set_label_scheme('dummy')
for xprop in ['C','S']:
  for sel in ['lo','med','high']:
    if xprop=='S':
      if sel=='lo':
        val=0.02
      elif sel=='med':
        val=0.1
      else:
        val=0.2
    elif xprop=='C':
      if sel=='lo':
        val=50
      elif sel=='med':
        val=70
      else:
        val=100
    graph2=grace.add_graph(Panel)
    graph2=SC.format_graph(graph2,xprop)
    graph2=SC.populate_persgraph(graph2,xprop,val,netprops,netSE)
    graph2.panel_label.configure(placement='iul',char_size=.5,dx=.02,dy=.02)


grace.multi(rows=2,cols=3,vgap=.07,hgap=.04)
grace.hide_redundant_labels()
grace.set_col_yaxislabel(col=0,rowspan=(None,None),label='Proportion of motif in participation vector',char_size=1,perpendicular_offset=.05)
grace.set_row_xaxislabel(row=0,colspan=(None,None),label='Connectance',char_size=1,perpendicular_offset=.04)
grace.set_row_xaxislabel(row=1,colspan=(None,None),label='Network Size',char_size=1,perpendicular_offset=.04)

grace.write_file('../manuscript/figures/roles_vs_SC_all.eps')
