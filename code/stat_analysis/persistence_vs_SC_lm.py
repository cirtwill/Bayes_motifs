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


scales={
  'Size':(79.1882,16.62225),
  'Connectance':(0.1044983,0.05559927),
  'Disturbance':(0.3,0.1366261),
  'Basal':(0.1601632,0.07986507)
}

lms={
  'Intercept': 0.3259301,
  'Size':-0.0007211,
  'Connectance':-0.0611595,
  'Disturbance':-0.5959284,
  'Size:Connectance':-0.0044764,
  'Size:Disturbance':-0.0162922,
  'Connectance:Disturbance':0.0045203,
  'Size:Connectance:Disturbance':-0.0013456 }

Blms={'Intercept':5.777e-01,
      'Basal':2.841e-02,
      'Size':5.309e-03,
      'Connectance':6.482e-03,
      'Disturbance':-1.452e-01,
      'Basal:Size':1.704e-03,
      'Basal:Connectance':6.680e-03,
      'Size:Connectance':2.302e-03,
      'Basal:Disturbance':-1.116e-02,
      'Size:Disturbance':-6.699e-03,
      'Connectance:Disturbance':-6.235e-03,
      'Basal:Size:Connectance':6.934e-04,
      'Basal:Size:Disturbance':-1.587e-03,
      'Basal:Connectance:Disturbance':-7.291e-03,
      'Size:Connectance:Disturbance':-2.105e-03,
      'Basal:Size:Connectance:Disturbance':-1.445e-03
}

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

def format_graph(graph,simple):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  graph.world.ymin=0
  graph.world.ymax=1
  graph.yaxis.tick.major=.2
  graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=.75)
  graph.xaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)

  if simple=='Deg':
    graph.xaxis.label.text='In-degree (number of prey)'
    graph.world.xmax=100
    graph.xaxis.tick.major=25
  elif simple=='TL':
    graph.xaxis.label.text='Trophic level (STL)'
    graph.world.xmax=5
    graph.xaxis.tick.major=1
  elif simple=='Size':
    graph.xaxis.label.text='Network size'
    graph.world.xmin=45
    graph.world.xmax=105
    graph.xaxis.tick.major=10
  elif simple=='Connectance':
    graph.xaxis.label.text='Connectance'
    graph.world.xmin=0
    graph.world.xmax=.2
    graph.xaxis.tick.major=.05
    graph.xaxis.ticklabel.configure(format='decimal',prec=2,char_size=.75)

  graph.xaxis.label.configure(char_size=1,just=2,place='normal')
  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.4,place='both',major_linewidth=.5,minor_linewidth=1)
  graph.yaxis.tick.configure(onoff='on',minor_ticks=1,major_size=.4,minor_size=.5,place='both',major_linewidth=.5,minor_linewidth=1)


  return graph

def populate_persgraph(graph,simple):

  if simple=='Size':
    levels=[50,60,70,80,90,100]
    altlevels=[0.02,0.1,0.18]
    altpred='Connectance'
  elif simple=='Connectance':
    levels=[0.02,0.06,0.1,0.14,0.18]
    altlevels=[50,75,100]
    altpred='Size'

  j=13
  for Disturbance in [0.1,0.18,0.26,0.34,0.42,0.5]:
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
        logy=math.exp(final)/(1+math.exp(final))
        dats.append((level,logy))

      data=graph.add_dataset(dats)
      if al==altlevels[0]:
        sty=2
      elif al==altlevels[1]:
        sty=1
      else:
        sty=5


      data.symbol.shape=0
      data.line.configure(linestyle=sty,linewdith=1.5,color=j)

      if simple=='Size' and Disturbance==0.1:
        data.legend='C='+str(al)        
        graph.legend.configure(char_size=.5,loc=(50,0.3),loctype='world',box_linestyle=0,fill_pattern=0)
      elif simple=='Connectance' and Disturbance==0.1:
        data.legend='Size='+str(al)
        graph.legend.configure(char_size=.5,loc=(.0175,0.3),loctype='world',box_linestyle=0,fill_pattern=0)
    j+=1

  return graph


def populate_specgraph(graph,persdict,simple):
  if simple=='Deg':
    key='scale(in_Degree)'
    scal=12.51637
    cent=11.95232
  else:
    key='scale(STL)'
    scal=0.4803683
    cent=2.274325

  print persdict.keys()

  j=13
  for basal_p in [0.1,0.18,0.26,0.34,0.42,0.5]:
    p=(basal_p-0.3)/0.1366261

    dats=[]
    # lower=[]
    # upper=[]
    for x in range(0,100):
      skex=(x-cent)/scal
      xcomp=persdict['(Intercept)'][0]+skex*persdict[key][0]
      bcomp=p*persdict['scale(Disturbance)'][0]+skex*p*persdict[key+':scale(Disturbance)'][0]
      y=xcomp+bcomp
      logity=math.exp(y)/(1+math.exp(y))

      dats.append((x,logity))

    data=graph.add_dataset(dats)
    data.symbol.shape=0
    data.line.configure(linestyle=1,linewidth=3.5,color=j)

    if simple=='TL':
      data.legend=str(basal_p)
    j+=1

  if simple=='TL':
    graph.add_drawing_object(DrawText,text='Basal species',x=5.25,y=0.71,loctype='world',char_size=.75) 
    graph.add_drawing_object(DrawText,text='extinction',x=5.25,y=0.65,loctype='world',char_size=.75) 
    graph.add_drawing_object(DrawText,text='probability',x=5.25,y=0.58,loctype='world',char_size=.75) 
    graph.legend.configure(char_size=.75,loc=(5.25,0.55),loctype='world',box_linestyle=0,fill_pattern=0)


  return graph


###############################################################################################
###############################################################################################
#
#           Assemble the plots
#
###############################################################################################
###############################################################################################
persdict=read_persfiles('persistence_vs_Deg_norandom.tsv','persistence_vs_TL.tsv')


grace=MultiPanelGrace(colors=colors)

for simple in ['Connectance','Deg','Size','TL']:
  graph2=grace.add_graph(Panel)
  graph2=format_graph(graph2,simple)
  if simple in ['Size','Connectance']:
    graph2=populate_persgraph(graph2,simple)
  else:
    graph2=populate_specgraph(graph2,persdict[simple],simple)
  graph2.panel_label.configure(placement='iul',char_size=1,dx=.03,dy=.03)
    

grace.multi(rows=2,cols=2,vgap=.09,hgap=.09)
grace.hide_redundant_labels()
grace.set_col_yaxislabel(label='Mean consumer probability of persistence',col=0,rowspan=(None,None),char_size=1.5,perpendicular_offset=.07)
grace.set_col_yaxislabel(label='Consumer probability of persistence',col=1,rowspan=(None,None),char_size=1.5,perpendicular_offset=.04)

grace.write_file('../../manuscript/figures/persistence_vs_SC_lm.eps')
