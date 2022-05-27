import sys
import os
import pymfinder as py

def main():

    webdir='../data/empirical/global_verts/edgelists/'
    roledir='../data/roles/empirical/global_verts/'
    motifdir='../data/motif_frequencies/empirical/global_verts/'
    pardir='../data/motif_participation/empirical/global_verts/'

    for web in os.listdir(webdir):
      infile=webdir+web
      print web

      size=3
      motfile=motifdir+web.split('.tsv')[0]+'.motifs_'+str(size)
      parfile=pardir+web.split('.tsv')[0]+'.participation_'+str(size)
      rolefile=roledir+web.split('.tsv')[0]+'.roles_'+str(size)

      result=py.motif_roles(infile,motifsize=size,stoufferIDs=False,allroles=True)
      motifs=str(result).split('node')[0]
      roles='node'+str(result).split('node')[1]

      presult=py.motif_participation(infile,motifsize=size,stoufferIDs=False,allmotifs=True)
      partic=str(presult).split('node')[1]

      e=open(parfile,'w')
      e.write(partic)
      e.close()

      f=open(rolefile,'w')
      f.write(roles)
      f.close()

      g=open(motfile,'w')
      g.write(motifs)
      g.close()

      # size=4
      # parfile2=pardir+web.split('.tsv')[0]+'.participation_'+str(size)          
      # presult=py.motif_participation(infile,motifsize=size,stoufferIDs=False,allmotifs=True)
      # partic=str(presult).split('node')[1]

      # e=open(parfile2,'w')
      # e.write(partic)
      # e.close()

if __name__ == '__main__':
  main()
