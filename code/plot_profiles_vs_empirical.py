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
      if S not in profiles:
        profiles[S]={}
      C=float(line.split()[1])
      if C not in profiles[S]:
        profiles[S][C]={'AC':[],'DC':[],'chain':[],'omnivory':[]}
      basal_p=float(line.split()[2])
      ID=line.split()[3]
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
      if basal_p==1.0 and deg>0 and S>50: # Only want one per network
        profiles[S][C]['AC'].append(float(AC)/float(total))
        profiles[S][C]['DC'].append(float(DC)/float(total))
        profiles[S][C]['chain'].append(float(chain)/float(total))
        profiles[S][C]['omnivory'].append(float(omnivory)/float(total))
  f.close()

  return profiles

def read_motif_profiles(motfile):
  profiles={}
  f=open(motfile,'r')
  for line in f:
    if line.split()[0]!='Size':
      S=int(line.split()[0])
      if S not in profiles:
        profiles[S]={}
      C=float(line.split()[1])
      if C not in profiles[S]:
        profiles[S][C]={'AC':[],'DC':[],'chain':[],'omnivory':[]}
      network=line.split()[2]
      apparent=int(line.split()[3])
      chain=int(line.split()[4])
      direct=int(line.split()[5])
      omni=int(line.split()[6])
      total=float(line.split()[7])

      profiles[S][C]['AC'].append(float(apparent)/total)
      profiles[S][C]['DC'].append(float(direct)/total)
      profiles[S][C]['chain'].append(float(chain)/total)
      profiles[S][C]['omnivory'].append(float(omni)/total)

  f.close()

  return profiles

def read_empirical_profiles(motfile):
  profiles={}
  f=open(motfile,'r')
  for line in f:
    if line.split()[0]!='Size':
      S=int(line.split()[0])
      C=float(line.split()[1])
      network=line.split()[2]
      if S>50:
        profiles[network]={'AC':[],'DC':[],'chain':[],'omnivory':[]}
      apparent=int(line.split()[3])
      chain=int(line.split()[4])
      direct=int(line.split()[5])
      omni=int(line.split()[6])
      total=float(line.split()[7])

      # if float(chain)/total<0.12:
      #   print network, S, C, float(chain)/total

      if S>50:
        profiles[network]['AC'].append(float(apparent)/total)
        profiles[network]['DC'].append(float(direct)/total)
        profiles[network]['chain'].append(float(chain)/total)
        profiles[network]['omnivory'].append(float(omni)/total)

  f.close()

  return profiles

def read_empirical_roles(motfile):
  profiles={}
  f=open(motfile,'r')
  for line in f:
    if line.split()[0]!='Size':
      S=int(line.split()[0])
      C=float(line.split()[1])
      basal_p=float(line.split()[2])
      network=line.split()[3]
      if network not in profiles:
        profiles[network]={'AC':[],'DC':[],'chain':[],'omnivory':[]}
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
      if basal_p==1.0 and deg>0 and total>0:
        if network in ['WEB246', 'WEB244', 'WEB245', 'WEB242', 'WEB240', 'WEB241', 'WEB181', 'WEB182', 'WEB210', 'WEB239', 'WEB238', 'WEB214', 'WEB217', 'WEB216', 'WEB233', 'WEB232', 'WEB230', 'WEB237', 'WEB236', 'WEB234', 'WEB359', 'WEB263', 'WEB219', 'WEB218', 'WEB229', 'WEB321', 'WEB320', 'WEB323', 'WEB322', 'WEB228', 'WEB324', 'WEB224', 'WEB225', 'WEB227', 'WEB222', 'WEB223']:
          profiles[network]['AC'].append(float(AC)/float(total))
          profiles[network]['DC'].append(float(DC)/float(total))
          profiles[network]['chain'].append(float(chain)/float(total))
          profiles[network]['omnivory'].append(float(omnivory)/float(total))
  f.close()


  return profiles

def format_graph(graph,value):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  graph.world.ymin=0
  if value=='profile':
    graph.world.ymax=.8
    graph.yaxis.tick.major=.1
  else:
    graph.world.ymax=1
    graph.yaxis.tick.major=.2
  graph.world.xmin=0
  graph.world.xmax=5

  graph.xaxis.label.configure(text='Motif',char_size=1,just=2,place='normal')
  graph.xaxis.tick.set_spec_ticks(major_ticks=[1,2,3,4],tick_labels=['AC','DC','Chain','Omni'],minor_ticks=[])
  graph.xaxis.ticklabel.configure(char_size=.75)
  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.4,minor_size=.5,place='both',major_linewidth=.5,minor_linewidth=1)

  graph.yaxis.label.configure(text='Proportion',char_size=1,just=2,place='normal')
  graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=.75)
  graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.4,minor_size=.5,place='both',major_linewidth=.5,minor_linewidth=1)

  return graph

def populate_graph(graph,sim_profiles,emp_profiles):

  motnames={'AC':1,'DC':2,'chain':3,'omnivory':4}

  for mot in ['AC','DC','chain','omnivory']:
    x=motnames[mot]
    ally=[]
    for S in sim_profiles:
      for C in sim_profiles[S]:
        for m in sim_profiles[S][C][mot]:
          ally.append(m)
    mean=np.mean(ally)
    sd=np.std(ally)
    # print sorted(ally), len(ally)
    lowQ=sorted(ally)[749]
    miny=sorted(ally)[0]
    highQ=sorted(ally)[-750]
    maxy=sorted(ally)[-1]

    # graph.add_drawing_object(DrawBox,lowleft=(x-.3,lowQ),upright=(x+.3,highQ),loctype='world',fill_pattern=0)
    graph.add_drawing_object(DrawBox,lowleft=(x-.25,miny),upright=(x+.25,maxy),loctype='world',fill_pattern=0,color=1,linewidth=1)
    dot=graph.add_dataset([(x-.35,mean,sd)],type='xydy')

  for net in sorted(emp_profiles,reverse=True):
    for mot in ['AC','DC','chain','omnivory']:
      y=emp_profiles[net][mot][0]
      x=motnames[mot]+random.random()/4-0.15
      dat=graph.add_dataset([(x,y)])
      dat.line.linestyle=0
      if net in ['WEB320','WEB321','WEB322','WEB323','WEB324']:
        dat.symbol.configure(color=18,fill_pattern=0,linewidth=1)
        if net=='WEB320' and mot=='AC':
          dat.legend='California coast'
      elif net in ['WEB273','WEB274','WEB275','WEB276']:
        dat.symbol.configure(color=6,fill_pattern=0,linewidth=1)
        if net=='WEB273' and mot=='AC':
          dat.legend='Alaskan stream'
      elif net=='WEB357':
        dat.symbol.configure(color=1,fill_pattern=1,linewidth=1)
        if mot=='AC':
          dat.legend='Benguela'
      elif net in ['WEB214','WEB215','WEB216','WEB217']:
        dat.symbol.configure(shape=2,color=4,fill_pattern=0,linewidth=1)
        if mot=='AC' and net=='WEB214':
          dat.legend='US Stream'
      elif net in ['WEB218','WEB219','WEB220','WEB221','WEB222','WEB223',
        'WEB224','WEB225','WEB226','WEB227','WEB228','WEB229','WEB230',
        'WEB231','WEB232','WEB233','WEB234','WEB235','WEB236','WEB237',
        'WEB238','WEB239','WEB240','WEB241','WEB242','WEB243','WEB244',
        'WEB245','WEB246','WEB247']:
        dat.symbol.configure(shape=2,color=2,fill_pattern=0,linewidth=1)
        if net=='WEB242' and mot=='AC':
          dat.legend='NZ Stream'
      else:
        if mot=='AC':
          print net
        dat.symbol.configure(shape=2,color=14,fill_pattern=0,linewidth=1)
        if net=='WEB181' and mot=='AC':
          dat.legend='Other invert-dominated'

  graph.legend.configure(char_size=.5,loc=(3.0,0.75),loctype='world',box_linestyle=0)

  # for ID in sim_profiles:
  #   xmean=np.mean(sim_profiles[ID][xtype])
  #   xse=np.std(sim_profiles[ID][xtype])/math.sqrt(len(sim_profiles[ID][xtype]))
  #   C=sim_profiles[ID]['C']
  #   ymean=np.mean(sim_profiles[ID][ytype])
  #   yse=np.std(sim_profiles[ID][ytype])/math.sqrt(len(sim_profiles[ID][ytype]))
  #   # print int(C*100)
  #   sims=graph.add_dataset([(xmean,ymean,xse,yse)],type='xydxdy')
  #   sims.line.linestyle=0
  #   sims.symbol.configure(fill_pattern=1,shape=1,color=int(100*C),fill_color=int(100*C),size=.25)
  #   sims.errorbar.configure(riser_linewidth=.5,linewidth=.5,size=.25,color=int(100*C))

  # for ID in emp_profiles:
  #   xmean=np.mean(emp_profiles[ID][xtype])
  #   xse=np.std(emp_profiles[ID][xtype])/math.sqrt(len(emp_profiles[ID][xtype]))
  #   C=emp_profiles[ID]['C']
  #   ymean=np.mean(emp_profiles[ID][ytype])
  #   yse=np.std(emp_profiles[ID][ytype])/math.sqrt(len(emp_profiles[ID][ytype]))
  #   print ID
  #   # print int(C*100)
  #   ems=graph.add_dataset([(xmean,ymean,xse,yse)],type='xydxdy')
  #   ems.line.linestyle=0
  #   ems.symbol.configure(fill_pattern=1,shape=2,fill_color=int(100*C),color=int(100*C))
  #   ems.errorbar.configure(riser_linewidth=.5,linewidth=.5,size=.25,color=int(100*C))

  #   # graph.add_drawing_object(DrawText,text=str.capitalize(ID.split('_')[0]),char_size=.55,x=emp_profiles[xtype][ytype][ID][0][0],
  #   #   y=emp_profiles[xtype][ytype][ID][0][1],loctype='world',just=2)


  return graph

def populate_pargraph(graph,sim_profiles,emp_profiles):

  motnames={'AC':1,'DC':2,'chain':3,'omnivory':4}

  for mot in ['AC','DC','chain','omnivory']:
    x=motnames[mot]
    ally=[]
    for S in sim_profiles:
      for C in sim_profiles[S]:
        for m in sim_profiles[S][C][mot]:
          ally.append(m)
    mean=np.mean(ally)
    sd=np.std(ally)
    # print sorted(ally), len(ally)
    lowQ=sorted(ally)[749]
    miny=sorted(ally)[0]
    highQ=sorted(ally)[-750]
    maxy=sorted(ally)[-1]

    # graph.add_drawing_object(DrawBox,lowleft=(x-.3,lowQ),upright=(x+.3,highQ),loctype='world',fill_pattern=0)
    graph.add_drawing_object(DrawBox,lowleft=(x-.25,miny),upright=(x+.25,maxy),loctype='world',fill_pattern=0,color=1,linewidth=1)
    dot=graph.add_dataset([(x-.35,mean,sd)],type='xydy')

  for net in sorted(emp_profiles,reverse=True):
    for mot in ['AC','DC','chain','omnivory']:
      for y in emp_profiles[net][mot]:
        x=motnames[mot]+random.random()/2-.25
        dat=graph.add_dataset([(x,y)])
        dat.line.linestyle=0
        if net in ['WEB320','WEB321','WEB322','WEB323','WEB324']:
          dat.symbol.configure(color=18,fill_pattern=0,linewidth=1,size=.05)
        elif net in ['WEB214','WEB215','WEB216','WEB217']:
          dat.symbol.configure(shape=2,color=4,fill_pattern=0,linewidth=1,size=.05)
        elif net in ['WEB218','WEB219','WEB220','WEB221','WEB222','WEB223',
          'WEB224','WEB225','WEB226','WEB227','WEB228','WEB229','WEB230',
          'WEB231','WEB232','WEB233','WEB234','WEB235','WEB236','WEB237',
          'WEB238','WEB239','WEB240','WEB241','WEB242','WEB243','WEB244',
          'WEB245','WEB246','WEB247']:
          dat.symbol.configure(shape=2,color=2,fill_pattern=0,linewidth=1,size=.05)
        else:
          dat.symbol.configure(shape=2,color=14,fill_pattern=0,linewidth=1,size=.05)

 
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

emp_profiles=read_empirical_profiles('../data/3sp_global_verts_motif_profiles.tsv')
emp_roles=read_empirical_roles('../data/global_verts_3sp_roles_participation_nonlinear.tsv')
sim_profiles=read_motif_profiles('../data/3sp_motif_profiles.tsv')
sim_roles=read_motif_roles('../data/3sp_roles_participation.tsv')

grace=MultiPanelGrace(colors=colors)
grace.add_label_scheme('dummy',['A: Motif profiles','B: Motif participation','C','D','Season length','F','G','H'])
grace.set_label_scheme('dummy')
for value in ['profile','participation']:
  graph=grace.add_graph(Panel)
  graph=format_graph(graph,value)
  graph.panel_label.configure(placement='ouc',char_size=.75,dx=.03,dy=.01)
  if value=='profile':
    graph=populate_graph(graph,sim_profiles,emp_profiles)
  else:
    graph=populate_pargraph(graph,sim_roles,emp_roles)

grace.multi(rows=1,cols=2,vgap=.07,hgap=.05)
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
grace.write_file('../manuscript/figures/motif_profiles_participation_vs_empirical.eps')

