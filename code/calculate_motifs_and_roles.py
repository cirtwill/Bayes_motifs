import sys
import os
import pymfinder as py

def main():

    webdir='../data/edgelists/acyclic/'
    motifdir='../data/motif_frequencies/acyclic/'
    roledir='../data/roles/acyclic/'
    pardir='../data/motif_participation/acyclic/'

    for s in sorted(os.listdir(webdir)):
      try:
        os.mkdir(roledir+s)
      except:
        pass
      try:
        os.mkdir(motifdir+s)
      except:
        pass
      try:
        os.mkdir(pardir+s)
      except:
        pass
      for c in sorted(os.listdir(webdir+s)):
        print s,c
        try:
          os.mkdir(roledir+s+'/'+c)
        except:
          pass
        try:
          os.mkdir(motifdir+s+'/'+c)
        except:
          pass
        try:
          os.mkdir(pardir+s+'/'+c)
        except:
          pass
        for web in os.listdir(webdir+s+'/'+c+'/'):
          infile=webdir+s+'/'+c+'/'+web

          size=3
          rolefile=roledir+s+'/'+c+'/'+web.split('.tsv')[0]+'.roles_'+str(size)
          motfile=motifdir+s+'/'+c+'/'+web.split('.tsv')[0]+'.motifs_'+str(size)
          parfile=pardir+s+'/'+c+'/'+web.split('.tsv')[0]+'.participation_'+str(size)

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

          size=4
          parfile2=pardir+s+'/'+c+'/'+web.split('.tsv')[0]+'.participation_'+str(size)          
          presult=py.motif_participation(infile,motifsize=size,stoufferIDs=False,allmotifs=True)
          partic=str(presult).split('node')[1]

          e=open(parfile2,'w')
          e.write(partic)
          e.close()

if __name__ == '__main__':
  main()
