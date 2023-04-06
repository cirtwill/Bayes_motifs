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
  'chain':(0.2489432,0.1272357),
  'apparent':(0.4051135,0.1461447),
  'direct':(0.1923225,0.1143021),
  'omnivory':(0.1536209,0.1007478),
  'Disturbance':(0.3,0.1366261)
}

lms={ # Intercept, motif, disturbance, interaction
  'chain':(0.316878,-0.009774,-0.598433,-0.015098),  
  'apparent':(0.316934,-0.002628,-0.598478,0.019535),
  'direct':(0.316977,0.030266,-0.598517,-0.001665),
  'omnivory':(0.316886,-0.018064,-0.598426,-0.007762)
}
# Disturbance effect probably wrong... should get approx. evenly spaced curves

# Empty dict to fake adding SE for error bars
zerodict={ # Intercept, motif, disturbance, interaction
  'chain':(0,0,0,0),
  'apparent':(0,0,0,0),
  'direct':(0,0,0,0),
  'omnivory':(0,0,0,0)
}

# SE for upper error bars
plusSE={ # Intercept, motif, disturbance, interaction
  'chain':(0.001997,0.002001,0.002057,0.002067),  
  'apparent':(0.001997,0.001994,0.002057,0.002054),
  'direct':(0.001997,0.002005,0.002058,0.002066),
  'omnivory':(0.001997,0.001996,0.002057,0.002058)
}

# negative SE for lower error bars
minusSE={ # Intercept, motif, disturbance, interaction
  'chain':(-0.001997,-0.002001,-0.002057,-0.002067),  
  'apparent':(-0.001997,-0.001994,-0.002057,-0.002054),
  'direct':(-0.001997,-0.002005,-0.002058,-0.002066),
  'omnivory':(-0.001997,-0.001996,-0.002057,-0.002058)
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
    graph.world.ymin=0
    graph.world.ymax=0.9
    graph.world.xmax=1
    graph.xaxis.tick.major=.3
    graph.yaxis.tick.major=.2
    graph.yaxis.label.text='Probability of persistence'
    graph.xaxis.label.text='Proportion of motif in network profile'


  graph.xaxis.ticklabel.configure(format='decimal',prec=1,char_size=.5)
  graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=.5)
  graph.xaxis.label.configure(char_size=1.5,just=2,place='normal')

  graph.yaxis.label.configure(char_size=1,just=2,place='normal')
  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.4,place='both',major_linewidth=.5,minor_linewidth=1)
  graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.4,place='both',major_linewidth=.5,minor_linewidth=1)
  graph.panel_label.configure(placement='iuc',char_size=.75,just=2,dy=.01)

  return graph

def populate_persgraph(graph,motif):

  j=13
  for dist in [0.1,0.18,0.26,0.34,0.42,0.5]:
  # for dist in [0.5,0.42,0.34,0.26,0.18,0.1]:

    sty=1
    scaldist=(dist-scales['Disturbance'][0])/scales['Disturbance'][1]

    # Will need to add error bars... 

    for lin in ['max','min','mean']:
      dats=[]
      if lin=='mean':
        offset=zerodict
      elif lin=='max':
        offset=plusSE
      else:
        offset=minusSE

      minx=int(mot_ranges[motif][0]*100)
      maxx=int(mot_ranges[motif][1]*100)

      # 5SE is a silly error bar, just using to test plotting.
      for x in range(minx,maxx):
        scalx=(float(x)/100-scales[motif][0])/scales[motif][1]
        y=lms[motif][0]+2*offset[motif][0]
        y+=scalx*(lms[motif][1]+2*offset[motif][1])
        y+=scaldist*(lms[motif][2]+2*offset[motif][2])
        y+=scalx*scaldist*(lms[motif][3]+2*offset[motif][3])

        logity=math.exp(y)/(1+math.exp(y))

        dats.append((float(x)/100,logity))

      data=graph.add_dataset(dats)
      data.symbol.shape=0
      if lin=='mean':
        data.line.configure(linestyle=sty,linewidth=2,color=j)
      else:
        data.line.linestyle=0
        data.fill.configure(color=j,pattern=4,type=2)
        if lin=='min':
          data.fill.color=0

    if motif=='apparent':
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

# # Far too many species to see anything. Going to apply stats.
# datafile='motif_proportions_persistence.tsv'
# netprops=read_datafile(datafile)

grace=MultiPanelGrace(colors=colors)
grace.add_label_scheme('dummy',['Apparent comp.','Three-sp. chain','Direct comp.','Omnivory'])
grace.set_label_scheme('dummy')

for mot in ['apparent','chain','direct','omnivory']:
  graph2=grace.add_graph(Panel)
  graph2=format_graph(graph2,'persistence')
  graph2=populate_persgraph(graph2,mot)

# graph2.add_drawing_object(DrawText,x=0.5,y=0.75,loctype='world',text='Lowest disturbance',just=0,char_size=1)

grace.graphs[0].add_drawing_object(DrawText,x=0.1,y=0.20,loctype='world',text='Basal species',just=0,char_size=.5)
grace.graphs[0].add_drawing_object(DrawText,x=0.1,y=0.18,loctype='world',text='extinction',just=0,char_size=.5)
grace.graphs[0].add_drawing_object(DrawText,x=0.1,y=0.16,loctype='world',text='probability',just=0,char_size=.5)
grace.graphs[0].legend.configure(char_size=.5,loc=(0.1,0.15),loctype='world',box_linestyle=0,fill_pattern=0)
grace.multi(rows=1,cols=4)
grace.hide_redundant_labels()
grace.graphs[0].set_view(0.15,0.15,0.33,0.95)
grace.graphs[1].set_view(0.35,0.15,0.53,0.95)
grace.graphs[2].set_view(0.55,0.15,0.73,0.95)
grace.graphs[3].set_view(0.75,0.15,0.93,0.95)
grace.set_row_xaxislabel(row=0,colspan=(None,None),label='Proportion of motif in role',char_size=1,perpendicular_offset=.05)

grace.write_file('../../manuscript/figures/persistence_motif_participation.eps')


