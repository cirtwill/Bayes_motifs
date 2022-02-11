import os
import sys

# igraph expects predators on rows/

def read_matrix(infile):
  species=set()
  linkdict={}
  f=open(infile,'r')
  for line in f:
    pred=line.split()[0]
    prey=line.split()[1]
    if pred not in linkdict:
      linkdict[pred]=[prey]
    else:
      linkdict[pred].append(prey)
    species.add(pred)
    species.add(prey)
  f.close()

  return species, linkdict

def write_matrix(species,linkdict,net):

  f=open('../data/empirical/matrix/'+net.split('.')[0]+'.tsv','w')
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

  indir='../data/empirical/original/'
  nets=os.listdir(indir)

  for net in nets:
    species, linkdict=read_matrix(indir+net)
    write_matrix(species, linkdict, net)


if __name__ == '__main__':
  main()
