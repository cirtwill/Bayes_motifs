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
  'omnivory':(0,0.5491713),
  'apparent':(0,1),
  'direct':(0,1),
  'chain':(0,1),
}

scales={ # center, scale
  'prop_chain':(0.2489432,0.1272357),
  'prop_apparent':(0.4051135,0.1461447),
  'prop_direct':(0.1923225,0.1143021),
  'prop_omni':(0.1536209,0.1007478),
  'Disturbance':(0.3,0.1366261),
  'Size':(79.1882,16.62225),
  'Connectance':(0.1044983,0.05559927),
  'Intercept':(0,1)
}

motif_names={'omnivory':'Omnivory','apparent':'Apparent competition','direct':'Direct competition','chain':'Three-species chain'}

# Each point will be mean of one network
def read_coefs(datafile):
  netprops={}
  f=open(datafile,'r')
  key=1
  for line in f:
    if line.split()[0]!='"Estimate"':
      ID=line.split()[0].split('"')[1]
      ID=ID.split(':')
      goodID=set()
      for term in ID:
        if '(' in term:
          newterm=term.split('(')[1].split(')')[0]
          goodID.add(newterm)

      effect=float(line.split()[1])
      netprops[key]={'terms':goodID,'coef':effect}
      key+=1

  f.close()
  return netprops

def format_graph(graph,S,C):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  graph.world.xmin=0
  graph.world.ymin=0

  graph.world.ymin=0
  graph.world.ymax=0.9
  graph.world.xmax=1
  graph.xaxis.tick.major=.3
  graph.yaxis.tick.major=.2
  if C==0.2:
    graph.yaxis.label.configure(text='S='+str(S),place='opposite',char_size=.75)
  if S==50:
    graph.xaxis.label.configure(text='C='+str(C),place='opposite',char_size=.75)

  # graph.xaxis.label.text='Proportion of motif in network profile'

  if S==100:
    graph.xaxis.ticklabel.configure(format='decimal',prec=1,char_size=.5)
  else:
    graph.xaxis.ticklabel.configure(format='decimal',prec=1,char_size=0)
  if C==0.02:
    graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=.5)
  else:
    graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=0)    

  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.4,place='both',major_linewidth=.5,minor_linewidth=1)
  graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.4,place='both',major_linewidth=.5,minor_linewidth=1)
  graph.panel_label.configure(placement='iul',char_size=.5,just=2,dy=.01,dx=.01)

  return graph

def populate_persgraph(graph,coefs,S,C):

  minx=int(mot_ranges['apparent'][0]*100)
  maxx=int(mot_ranges['apparent'][1]*100)

  # print(coefs)
  # sys.exit()
  j=13
  for dist in [0.1,0.18,0.26,0.34,0.42,0.5]:
  # for dist in [0.5,0.42,0.34,0.26,0.18,0.1]:
    dats=[]
    sty=1
    scalDist=(dist-scales['Disturbance'][0])/scales['Disturbance'][1]
    scalS=(S-scales['Size'][0])/scales['Size'][1]
    scalC=(C-scales['Connectance'][0])/scales['Connectance'][1]

    # print(coefs[17])
    for x in range(minx,maxx):
      scalx=(float(x)/100-scales['prop_apparent'][0])/scales['prop_apparent'][1]
      y=coefs[1]['coef'] # Global intercept
      y+=coefs[2]['coef']*scalx # Apparent
      y+=coefs[3]['coef']*scalDist # Disturbance
      y+=coefs[4]['coef']*scalS # Size
      y+=coefs[5]['coef']*scalC # Connectance
      y+=coefs[6]['coef']*scalDist*scalx # Disturbance, apparent
      y+=coefs[7]['coef']*scalx*scalS # apparent, size
      y+=coefs[8]['coef']*scalDist*scalS # Disturbance, size
      y+=coefs[9]['coef']*scalx*scalC # apparent, Connectance
      y+=coefs[10]['coef']*scalDist*scalC # Connectance, Disturbance
      y+=coefs[11]['coef']*scalS*scalC # Connectance, Size
      y+=coefs[12]['coef']*scalS*scalDist*scalx # apparent, Disturbance, Size
      y+=coefs[13]['coef']*scalC*scalDist*scalx # apparent, Disturbance, Connectance
      y+=coefs[14]['coef']*scalC*scalS*scalx # apparent, Size, Connectance 
      y+=coefs[15]['coef']*scalC*scalS*scalDist # Disturbance, Size, Connectance 
      y+=coefs[16]['coef']*scalC*scalS*scalDist*scalx # Disturbance, Size, Connectance, apparent 

      logity=math.exp(y)/(1+math.exp(y))

      dats.append((float(x)/100,logity))

    data=graph.add_dataset(dats)
    data.symbol.shape=0
    data.line.configure(linestyle=sty,linewidth=2,color=j)

    if S==50 and C==0.02:
      data.legend=str(dist)
    j+=1
  return graph

###############################################################################################
###############################################################################################
#
#           Assemble the plots
#
###############################################################################################
###############################################################################################

datfile='apparent_lm_SC.tsv'
coefs=read_coefs(datfile)

# # Far too many species to see anything. Going to apply stats.
# datafile='motif_proportions_persistence.tsv'
# netprops=read_datafile(datafile)

grace=MultiPanelGrace(colors=colors)
# grace.add_label_scheme('dummy',['Apparent comp.','Three-sp. chain','Direct comp.','Omnivory'])
# grace.set_label_scheme('dummy')

for S in [50,100]:
  for C in [0.02,0.1,0.2]:
    graph2=grace.add_graph(Panel)
    graph2=format_graph(graph2,S,C)
    graph2=populate_persgraph(graph2,coefs,S,C)

# graph2.add_drawing_object(DrawText,x=0.5,y=0.75,loctype='world',text='Lowest disturbance',just=0,char_size=1)

# grace.graphs[0].add_drawing_object(DrawText,x=0.1,y=0.20,loctype='world',text='Basal species',just=0,char_size=.5)
# grace.graphs[0].add_drawing_object(DrawText,x=0.1,y=0.18,loctype='world',text='extinction',just=0,char_size=.5)
# grace.graphs[0].add_drawing_object(DrawText,x=0.1,y=0.16,loctype='world',text='probability',just=0,char_size=.5)
grace.graphs[0].legend.configure(char_size=.5,loc=(0.05,0.5),loctype='world',box_linestyle=0,fill_pattern=0)
grace.multi(rows=2,cols=3,hgap=.05,vgap=.05)
# grace.hide_redundant_labels()
# grace.graphs[0].set_view(0.15,0.15,0.33,0.95)
# grace.graphs[1].set_view(0.35,0.15,0.53,0.95)
# grace.graphs[2].set_view(0.55,0.15,0.73,0.95)
# grace.graphs[3].set_view(0.75,0.15,0.93,0.95)
grace.set_row_xaxislabel(row=1,colspan=(None,None),label='Proportion of apparent competition motif in role',char_size=1,perpendicular_offset=.05)
grace.set_col_yaxislabel(col=0,rowspan=(None,None),label='Probability of persistence',char_size=1,perpendicular_offset=.05)

grace.write_file('../../manuscript/figures/persistence_apparent_detailpers.eps')


