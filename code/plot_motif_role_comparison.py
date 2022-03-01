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

colors=ColorBrewerScheme('Spectral', reverse=True,n=18)  # The blue is very beautiful but maybe harder to see.
colors.add_color(173,132,191,'Purple')
colors.add_color(34,90,34,'Green')
# colors.add_color(120,120,120,'grey')
# colors.add_color(255,125,125,'lightish_red')
# colors.add_color(200,200,200,'lightgrey')

# Want to use AC, DC as x's, chain, omni as y's

def read_motif_roles(motfile):
  profiles={}
  f=open(motfile,'r')
  for line in f:
    if line.split()[0]!='Size':
      S=int(line.split()[0])
      C=float(line.split()[1])
      basal_p=float(line.split()[2])
      ID=line.split()[3]
      if ID not in profiles:
        profiles[ID]={'AC':[],'DC':[],'chain':[],'omnivory':[],'S':S,'C':C}
      sp=line.split()[4]
      persis=line.split()[5]
      deg=float(line.split()[6])
      out_deg=line.split()[7]
      STL=line.split()[8]
      PATL=line.split()[9]
      chain=int(line.split()[10])
      DC=int(line.split()[11])
      omnivory=int(line.split()[12])
      AC=int(line.split()[13])
      total=chain+DC+omnivory+AC
      if basal_p==1.0 and deg>0: # Only want one per network
        profiles[ID]['AC'].append(float(AC)/float(total))
        profiles[ID]['DC'].append(float(DC)/float(total))
        profiles[ID]['chain'].append(float(chain)/float(total))
        profiles[ID]['omnivory'].append(float(omnivory)/float(total))
  f.close()

  return profiles

def read_empirical_motif_roles(motfile):
  profiles={}
  f=open(motfile,'r')
  for line in f:
    if line.split()[0]!='Size':
      S=int(line.split()[0])
      C=float(line.split()[1])
      basal_p=float(line.split()[2])
      ID=line.split()[3]
      if ID not in profiles:
        profiles[ID]={'AC':[],'DC':[],'chain':[],'omnivory':[],'S':S,'C':C}
      sp=line.split()[4]
      persis=line.split()[5]
      deg=float(line.split()[6])
      out_deg=line.split()[7]
      STL=line.split()[8]
      chain=float(line.split()[9])
      DC=float(line.split()[10])
      omnivory=float(line.split()[11])
      AC=float(line.split()[12])
      total=chain+DC+omnivory+AC
      if basal_p==1.0 and deg>0: # Only want one per network
        profiles[ID]['AC'].append(float(AC)/float(total))
        profiles[ID]['DC'].append(float(DC)/float(total))
        profiles[ID]['chain'].append(float(chain)/float(total))
        profiles[ID]['omnivory'].append(float(omnivory)/float(total))
  f.close()

  return profiles

def format_graph(graph,xtype,ytype):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  if xtype=='DC':
    graph.world.xmin=0.05
  else:
    graph.world.xmin=0.25

  graph.world.xmax=.600000001
  graph.xaxis.tick.major=.1
  graph.xaxis.ticklabel.configure(format='decimal',prec=1,char_size=.75)
  graph.xaxis.label.configure(text='Frequency of '+xtype,char_size=1,just=2,place='normal')
  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.4,minor_size=.5,place='both',major_linewidth=.5,minor_linewidth=1)


  graph.world.ymin=0
  graph.world.ymax=.400000001
  graph.yaxis.tick.major=.1
  graph.yaxis.label.configure(text='Frequency of '+str.capitalize(ytype),char_size=1,just=2,place='normal')
  graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=.75)
  graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.4,minor_size=.5,place='both',major_linewidth=.5,minor_linewidth=1)

  return graph

def populate_graph(graph,xtype,ytype,sim_profiles,emp_profiles):

  for ID in sim_profiles:
    xmean=np.mean(sim_profiles[ID][xtype])
    xse=np.std(sim_profiles[ID][xtype])/math.sqrt(len(sim_profiles[ID][xtype]))
    C=sim_profiles[ID]['C']
    ymean=np.mean(sim_profiles[ID][ytype])
    yse=np.std(sim_profiles[ID][ytype])/math.sqrt(len(sim_profiles[ID][ytype]))
    # print int(C*100)
    sims=graph.add_dataset([(xmean,ymean,xse,yse)],type='xydxdy')
    sims.line.linestyle=0
    sims.symbol.configure(fill_pattern=1,shape=1,color=int(100*C),fill_color=int(100*C),size=.25)
    sims.errorbar.configure(riser_linewidth=.5,linewidth=.5,size=.25,color=int(100*C))

  for ID in emp_profiles:
    xmean=np.mean(emp_profiles[ID][xtype])
    xse=np.std(emp_profiles[ID][xtype])/math.sqrt(len(emp_profiles[ID][xtype]))
    C=emp_profiles[ID]['C']
    ymean=np.mean(emp_profiles[ID][ytype])
    yse=np.std(emp_profiles[ID][ytype])/math.sqrt(len(emp_profiles[ID][ytype]))
    print ID
    # print int(C*100)
    ems=graph.add_dataset([(xmean,ymean,xse,yse)],type='xydxdy')
    ems.line.linestyle=0
    ems.symbol.configure(fill_pattern=1,shape=2,fill_color=int(100*C),color=int(100*C))
    ems.errorbar.configure(riser_linewidth=.5,linewidth=.5,size=.25,color=int(100*C))

    # graph.add_drawing_object(DrawText,text=str.capitalize(ID.split('_')[0]),char_size=.55,x=emp_profiles[xtype][ytype][ID][0][0],
    #   y=emp_profiles[xtype][ytype][ID][0][1],loctype='world',just=2)


  return graph

def add_legend(graph):
  for C in [0.02,0.06,0.1,0.14,0.18]:
    dat=graph.add_dataset([(100,100)])
    dat.line.linestyle=0
    dat.symbol.configure(size=.5,shape=1,color=int(100*C),fill_color=int(100*C))
    dat.legend=str(C)

  graph.legend.configure(loc=(0.45,0.3),loctype='world',char_size=.5,box_linestyle=0,fill_pattern=0)
  graph.add_drawing_object(DrawText,text='Connectance',x=0.45,y=0.305,loctype='world',char_size=.6)

  return graph

###############################################################################################
###############################################################################################
#
#           Assemble the plots
#
###############################################################################################
###############################################################################################

sim_profiles=read_motif_roles('../data/3sp_roles_participation.tsv')
emp_profiles=read_empirical_motif_roles('../data/empirical_3sp_roles_participation_nonlinear.tsv')

grace=MultiPanelGrace(colors=colors)
grace.add_label_scheme('dummy',['A','B','C','D','Season length','F','G','H'])
grace.set_label_scheme('dummy')
for x in ['AC','DC']:
  for y in ['chain','omnivory']:
    graph=grace.add_graph(Panel)
    graph=format_graph(graph,x,y)
    graph=populate_graph(graph,x,y,sim_profiles,emp_profiles)
    graph.panel_label.configure(placement='iul',char_size=1,dx=.03,dy=.01)
    if x=='DC' and y=='omnivory':
      graph=add_legend(graph)

grace.multi(rows=2,cols=2,vgap=.07,hgap=.07)
grace.hide_redundant_labels()
# grace.set_col_yaxislabel(label='Frequency of three-species chain',col=0,rowspan=(None,None),char_size=1,perpendicular_offset=.07,just=2)
# grace.set_col_yaxislabel(label='Frequency of omnivory',col=1,rowspan=(None,None),char_size=1,perpendicular_offset=.03,just=2)
# grace.set_row_xaxislabel(label='Frequency of direct competition',row=1,colspan=(None,None),char_size=1,perpendicular_offset=.05,just=2)
# grace.set_row_xaxislabel(label='Day of year',row=2,colspan=(None,None),char_size=1,perpendicular_offset=.05)
# for graph in grace.graphs:
#   print graph.get_view()
# grace.graphs[0].set_view(0.15,0.45,0.53,0.70)
# grace.graphs[1].set_view(0.57,0.45,0.95,0.70)
# grace.graphs[2].set_view(0.15,0.15,0.53,0.40)
# grace.graphs[3].set_view(0.57,0.15,0.95,0.40)
grace.write_file('../manuscript/figures/empirical_mean_roles.eps')
