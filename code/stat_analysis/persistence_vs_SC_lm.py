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
  elif simple=='Connectance':
    graph.xaxis.label.text='Connectance'
    graph.world.xmin=0
    graph.world.xmax=.2
    graph.xaxis.tick.major=.05
    graph.xaxis.ticklabel.configure(format='decimal',prec=2,char_size=.75)
  else:
    graph.xaxis.label.text='Proportion basal'
    graph.world.xmin=0
    graph.world.xmax=.5
    graph.xaxis.tick.major=.1
    graph.xaxis.ticklabel.configure(format='decimal',prec=1,char_size=.75)

  graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=.75)
  graph.xaxis.label.configure(char_size=1,just=2,place='normal')

  graph.yaxis.label.configure(char_size=1,just=2,place='normal')
  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.4,place='both',major_linewidth=.5,minor_linewidth=1)
  graph.yaxis.tick.configure(onoff='on',minor_ticks=1,major_size=.4,minor_size=.5,place='both',major_linewidth=.5,minor_linewidth=1)

  if simple=='dummy':
    graph.yaxis.bar.linestyle=0
    graph.xaxis.bar.linestyle=0
    graph.frame.linestyle=0
    graph.yaxis.tick.onoff='off'
    graph.xaxis.tick.onoff='off'
    graph.xaxis.label.char_size=0
    graph.yaxis.label.char_size=0
    graph.xaxis.ticklabel.char_size=0
    graph.yaxis.ticklabel.char_size=0
    graph.panel_label.char_size=0

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
        logy=math.exp(final)/(1+math.exp(final))
        dats.append((level,logy))

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

def populate_graph(graph,simple,dist):

  scalD=(dist-scales['Disturbance'][0])/scales['Disturbance'][1]
  leveldict={
    'Size':[50,60,70,80,90,100],
    'Connectance':[0.02,0.06,0.1,0.14,0.18],
    'Basal':[0.02,0.1,0.2,0.3,0.4,0.5]
  }

  if simple=='Size':
    alt1='Basal'
    alt2='Connectance'
  elif simple=='Connectance':
    alt1='Basal'
    alt2='Size'
  else:
    alt1='Size'
    alt2='Connectance'

  for pred1 in [leveldict[alt1][0],leveldict[alt1][1]]:
    scal1=(pred1-scales[alt1][0])/scales[alt1][1]
    base=lms['Intercept']+scalD*Blms['Disturbance']+scal1*Blms[alt1]+scalD*scal1*Blms[alt1+':Disturbance']

    for pred2 in [leveldict[alt2][0],leveldict[alt2][1]]:
      scal2=(pred2-scales[alt2][0])/scales[alt2][1]
      base2=scal2*Blms[alt2]+scalD*scal2*Blms[alt2+':Disturbance']+scal1*scal2*Blms[alt1+':'+alt2]+scalD*scal1*scal2*Blms[alt1+':'+alt2+':Disturbance']

      dats=[]
      for level in leveldict[simple]:
        scal3=(level-scales[simple][0])/scales[simple][1]
        if simple=='Size':
          base3=scal3*Blms[simple]+scalD*scal3*Blms[simple+':Disturbance']+scal1*scal3*Blms[alt1+':'+simple]+scal2*scal3*Blms[simple+':'+alt2]
          base4=scalD*scal1*scal3*Blms[alt1+':'+simple+':Disturbance']+scalD*scal2*scal3*Blms[simple+':'+alt2+':Disturbance']
        elif simple=='Connectance':
          base3=scal3*Blms[simple]+scalD*scal3*Blms[simple+':Disturbance']+scal1*scal3*Blms[alt1+':'+simple]+scal2*scal3*Blms[alt2+':'+simple]
          base4=scalD*scal1*scal3*Blms[alt1+':'+simple+':Disturbance']+scalD*scal2*scal3*Blms[alt2+':'+simple+':Disturbance']
        else:
          base3=scal3*Blms[simple]+scalD*scal3*Blms[simple+':Disturbance']+scal1*scal3*Blms[simple+':'+alt1]+scal2*scal3*Blms[simple+':'+alt2]
          base4=scalD*scal1*scal3*Blms[simple+':'+alt1+':Disturbance']+scalD*scal2*scal3*Blms[simple+':'+alt2+':Disturbance']

        final=base+base2+base3+base4+scal1*scal2*scal3*Blms['Basal:Size:Connectance']+scalD*scal1*scal2*scal3*Blms['Basal:Size:Connectance:Disturbance']
        dats.append((level,final))

      data=graph.add_dataset(dats)
      if pred1==leveldict[alt1][0]:
        sty=2
      else:
        sty=1
      if pred2==leveldict[alt2][0]:
        col=2
      else:
        col=14
      print col
      data.symbol.shape=0
      data.line.configure(linestyle=sty,linewdith=1.5,color=col)

  #     if simple=='Size' and al==altlevels[1]:
  #       data.legend='p(Basal extinct)='+str(Disturbance)
  #     # Doesn't work - add DrawText.
  #     # if simple=='Size' and Basal_p==0.5:
  #     #   data.legend='C='+str(al)        

  # if simple=='Size':
  #   # graph.add_drawing_object(DrawText,text='p(Basal)',x=107,y=0.825,loctype='world',char_size=.5)
  #   graph.legend.configure(char_size=.5,loc=(50,0.3),loctype='world',box_linestyle=0,fill_pattern=0)
  #   # graph.add_drawing_object(DrawText,text='Connectance',x=107,y=0.5,loctype='world',char_size=.5)

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

## I don't think we use this figure
# grace=MultiPanelGrace(colors=colors)

# for simple in ['Size','Connectance','Basal']:
#   for dist in [0.1, 0.5]:
#     graph2=grace.add_graph(Panel)
#     graph2=format_graph(graph2,simple)
#     if simple=='Size':
#       graph2.add_drawing_object(DrawText,text='Disturbance: '+str(dist),char_size=0.75,x=75,y=1.1,loctype='world',just=2)
#     graph2=populate_graph(graph2,simple,dist)
#     graph2.panel_label.configure(placement='iul',char_size=.75,dx=.02,dy=.02)
    

# # dummy=grace.add_graph(Panel)
# # dummy=format_graph(dummy,'dummy')

# grace.multi(rows=3,cols=2,vgap=.09,hgap=.07)
# grace.hide_redundant_labels()
# grace.set_row_xaxislabel(label='Network size',row=0,colspan=(None,None),char_size=1,perpendicular_offset=.05)
# grace.set_row_xaxislabel(label='Connectance',row=1,colspan=(None,None),char_size=1,perpendicular_offset=.05)
# grace.set_row_xaxislabel(label='Proportion basal',row=2,colspan=(None,None),char_size=1,perpendicular_offset=.05)
# grace.set_col_yaxislabel(label='Mean Persistence',col=0,rowspan=(None,None),char_size=1,perpendicular_offset=.07)
# # grace.set_col_yaxislabel(label='High disturbance',col=1,rowspan=(None,None),char_size=1,perpendicular_offset=.07)

# grace.write_file('../../manuscript/figures/persistence_vs_BSC_lm.eps')

