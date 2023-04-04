import os
import sys

# igraph expects predators on rows/

def read_matrix(infile):
  species=set()
  linkdict={}
  f=open(infile,'r')
  for line in f:
    comsplit=line.split(',')
    if 'WEB239' in infile:
      print comsplit
    if comsplit[-5:]==['','','','','\r\n']: # Ignore headers
      pass
    else:
      prey=line.split('\t')[0]
      links=line.split('\t')[1:]
      if '0' in links or 0 in links:
        species.add(prey)
        for i in range(0,len(links)):
          pred=preds[i].split('\n')[0]
          if 'WEB239' in infile:
            print pred
          species.add(pred)
          if float(links[i].split('\n')[0])>0:
            try:
              linkdict[pred].append(prey)
            except KeyError:
              linkdict[pred]=[prey]
      else:
        preds=links
  f.close()

  if 'WEB239' in infile:
    print species

  return species, linkdict

def write_matrix(species,linkdict,net):

  f=open('../data/empirical/global_verts/matrix/'+net.split('_mod')[0]+'.tsv','w')
  f.write('\t'+'\t'.join(sorted(species))+'\n')
  for pred in sorted(species):
    f.write(pred)
    if pred in linkdict:
      for prey in sorted(species):
        if prey in linkdict[pred]:
          f.write('\t1')
        else:
          f.write('\t0')
      f.write('\n')
    else:
      for spec in species:
        f.write('\t0')
      f.write('\n')      

  f.close()



def main():

  # indir='../data/empirical/global_verts/original/'
  indir='../data/empirical/global_inverts/original/'

  nets=os.listdir(indir)
  i=1
  for net in [n for n in nets if '_mod.csv' in n]:
    print net,i
    species, linkdict=read_matrix(indir+net)
    write_matrix(species, linkdict, net)
    i+=1
if __name__ == '__main__':
  main()
