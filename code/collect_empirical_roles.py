import os
import sys

def CSreader(roledir): # Read the edgelists, calculate S and C, make a dict.
  Cdict={}
  for web in os.listdir(roledir):
    S=set()
    L=set()
    f=open(roledir+web,'r')
    for line in f:
      (pred,prey)=line.split()
      S.add(pred)
      S.add(prey)
      L.add((pred,prey))
    f.close()
    C=float(len(L))/float(len(S)**2)
    Cdict[web.split('.tsv')[0]]=((len(S),C))
  return Cdict

# Then bin the species according to ... things, like in Benno's paper.
# Then do some stats and send to Anna.
def read_persistence_file(persfile,role_dict):
  persdict={}
  TLdict={}
  degdict={}
  bpdict={}
  for webfile in role_dict:
    web=webfile.split('.')[0]
    persdict[web]={}
    TLdict[web]={}
    degdict[web]={}
    bpdict[web]={}
    for bp in [0.0,0.2,0.4,0.6,0.8,1.0]:
      persdict[web][bp]={} 
      TLdict[web][bp]={}
      degdict[web][bp]={}

  f=open(persfile,'r')
  for line in f:
    if line.split()[0]!='"per"':
      per=float(line.split()[1])
      sp='"sp'+''.join(line.split()[2].split('"'))+'"'
      indeg=int(line.split()[3])
      outdeg=int(line.split()[4])
      try:
        STL=int(line.split()[5])
      except ValueError:
        STL=10000
      basal_p=float(line.split()[6])
      web=line.split()[7].split('.')[0].split('"')[1]
      persdict[web][basal_p][sp]=(str(per))
      degdict[web][basal_p][sp]=(str(indeg),str(outdeg))
      TLdict[web][basal_p][sp]=(str(STL))
  f.close()

  return TLdict,degdict,persdict


# partic_dict takes s, c, network, species, size, motif: count
def read_participation_file(pardir):
  par_dict={}
  full_motifs={3:set(),4:set()}
  size=3
  for webfile in os.listdir(pardir):
    web=webfile.split('.')[0]
    par_dict[web]={}
    f=open(pardir+webfile,'r')
    fixed_header=['motif']
    for line in f:
      if len(line.split())==13:
        for mot in line.split():
          fixed_header.append(mot)
          full_motifs[size].add(mot)
      elif len(line.split())>13:
        spnm='"sp'+''.join(line.split()[0].split('"'))+'"'
        par_dict[web][spnm]={}
        par_dict[web][spnm][size]={}
        for j in range(1,len(fixed_header)):
          par_dict[web][spnm][size][fixed_header[j]]=line.split()[j]  
    f.close()
  return par_dict, full_motifs


# Creates a dictionary with position counts listed by species, connectance, network, species.
def read_role_file(roledir):
  role_dict={}
  full_positions=set()
  for webfile in os.listdir(roledir):
    web=webfile.split('.')[0]
    role_dict[web]={}
    if '.roles_3' in webfile:
      net=roledir+webfile
      f=open(net,'r')
      for line in f:
        if len(line.split())>0:
          if line.split()[0]=='node':
            header=line.split('(')
            fixed_header=['node']
            for pos in header[1:]:
              position='('+pos[:-1]
              fixed_header.append(position)
              full_positions.add(position)
          else:
            sp='"sp'+''.join(line.split()[0].split('"'))+'"'
            role_dict[web][sp]={}
            for j in range(1,len(fixed_header)):
              role_dict[web][sp][fixed_header[j]]=line.split()[j]

      f.close()
  return role_dict, full_positions


def sort_nonzero_positions(partic_dict,full_motifs):
  good_motifs={3:set(),4:set()}
  motif_counts={}
  for size in full_motifs:
    for motif in full_motifs[size]:
      motif_counts[motif]=set()
  for size in full_motifs:
    for motif in full_motifs[size]:
      for web in partic_dict:
        for sp in partic_dict[web]:
          motif_counts[motif].add(partic_dict[web][sp][size][motif])
      if motif_counts[motif]!=set(['0']):
        good_motifs[size].add(motif)
  return good_motifs


def sort_nonzero_roles(role_dict,full_positions):
  good_motifs=set()
  position_counts={}
  for motif in full_positions:
    position_counts[motif]=set()
    for web in role_dict:
      for sp in role_dict[web]:
        position_counts[motif].add(role_dict[web][sp][motif])
    if position_counts[motif]!=set(['0']):
      good_motifs.add(motif)
  return good_motifs


def collect_and_print_3sp(outfile,role_dict,Cdict,good_positions,partic_dict,good_motifs,degdict,TLdict,persdict):

  nomotif_species=set()
  f=open(outfile,'w')
  f.write('Size\tConnectance\tBasal_p\tNetwork\tSpecies\tPersistence\tin_Degree\tout_Degree\tSTL\tm')
  f.write('\tm'.join(sorted(good_motifs[3]))+'\tp')
  f.write('\tp'.join(sorted(good_positions))+'\n')
  for web in persdict:
    print web
    for bp in [0.0,0.2,0.4,0.6,0.8,1.0]:
      for sp in persdict[web][bp]: # Seems to be only non-basal species?
        # Particdict doesn't include all species - some may not be in 3sp motifs
        # Especially in webs
        f.write(str(Cdict[web][0])+'\t'+str(Cdict[web][1])+'\t')
        f.write('\t'.join([str(bp),web,sp]))
        f.write('\t'+persdict[web][bp][sp])
        f.write('\t'+'\t'.join(degdict[web][bp][sp]))
        f.write('\t'+str(TLdict[web][bp][sp]))
        for motif in sorted(good_motifs[3]):
          if sp in partic_dict[web]:
            f.write('\t'+str(partic_dict[web][sp][3][motif]))
          else:
            f.write('\t0')
            nomotif_species.add((web,sp))
        for pos in sorted(good_positions):
          if sp in partic_dict[web]:
            f.write('\t'+str(role_dict[web][sp][pos]))
          else:
            f.write('\t0')
        f.write('\n')
  f.close()
  for web in sorted(persdict):
    missing=[]
    for (w,s) in nomotif_species:
      if w==web:
        missing.append(s)
    print web, len(missing), "Not present in 3sp motifs"
  print len(persdict.keys())
  return

def main():
  Cdict=CSreader('../data/empirical/global_verts/edgelists/')
  # role_dict takes s, c, network, species, position: count
  # full_positions is a full set of positions
  role_dict,full_positions=read_role_file('../data/roles/empirical/global_verts/')
  good_positions=sort_nonzero_roles(role_dict,full_positions)
  print 'roles read'
  # partic_dict takes s, c, network, species, size, motif: count
  # full_motifs is a dict for size: motif
  partic_dict,full_motifs=read_participation_file('../data/motif_participation/empirical/global_verts/')
  good_motifs=sort_nonzero_positions(partic_dict,full_motifs)
  print 'participation read'
  # TLdict,degdict=read_TL_file('../data/TL/',role_dict) # Role dict to give a reference of i's
  # print 'TLs and degrees read'


  linear_TLdict,linear_degdict,linear_perdict=read_persistence_file('../data/global_verts_linear.tsv',partic_dict)
  print 'persistence read'
  collect_and_print_3sp('../data/global_verts_3sp_roles_participation_linear.tsv',role_dict,Cdict,good_positions,partic_dict,good_motifs,linear_degdict,linear_TLdict,linear_perdict)
  print '3sp role files complete'

  nonlinear_TLdict,nonlinear_degdict,nonlinear_perdict=read_persistence_file('../data/global_verts_nonlinear.tsv',partic_dict)
  print 'persistence read'
  collect_and_print_3sp('../data/global_verts_3sp_roles_participation_nonlinear.tsv',role_dict,Cdict,good_positions,partic_dict,good_motifs,nonlinear_degdict,nonlinear_TLdict,nonlinear_perdict)
  print 'nonlinear 3sp role files complete'
      

if __name__ == '__main__':
  main()
