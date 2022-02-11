import os
import sys

listdir='../data/edgelists/acyclic/'
Annas_to_mine={}
for S in os.listdir(listdir):
  Annas_to_mine[S]={}
  for C in os.listdir(listdir+S):
    i=1
    Annas_to_mine[S][C]={}
    myfils=os.listdir(listdir+S+'/'+C)
    for filname in sorted(myfils):
      webnum=filname.split('_')[-1].split('.tsv')[0]
      Annas_to_mine[S][C][str(i)]=webnum
      i=i+1

# Need to get Anna to redo "all_new" with web numbers from new webs, not sequentially.

# Also have motif frequencies per web, may be interesting to compare to S, C?

# Then bin the species according to ... things, like in Benno's paper.
# Then do some stats and send to Anna.
def read_persistence_file(persfile,role_dict):
  persdict={}
  TLdict={}
  degdict={}
  bpdict={}
  for s in role_dict:
    persdict[s]={}
    TLdict[s]={}
    degdict[s]={}
    bpdict[s]={}
    for c in role_dict[s]:
      persdict[s][c]={} 
      TLdict[s][c]={}
      degdict[s][c]={}
      for bp in [0.0,0.2,0.4,0.6,0.8,1.0]:
        persdict[s][c][bp]={} 
        TLdict[s][c][bp]={}
        degdict[s][c][bp]={}
        for i in role_dict[s][c]:
          persdict[s][c][bp][i]={}
          TLdict[s][c][bp][i]={}
          degdict[s][c][bp][i]={} 

  f=open(persfile,'r')
  for line in f:
    if line.split()[0]!='"sp_no"':
      sp='"sp'+line.split()[1][1:-1]+'"'
      indeg=int(line.split()[2])
      outdeg=int(line.split()[3])
      STL=int(line.split()[4])
      PATL=float(line.split()[5])
      per=float(line.split()[6])
      basal_p=float(line.split()[7])
      S=str(int(line.split()[8]))
      C=str(float(line.split()[9]))
      web=str(int(line.split()[10]))
      Aweb=Annas_to_mine[S][C][web]
      persdict[int(S)][C][basal_p][Aweb][sp]=(str(per))
      degdict[int(S)][C][basal_p][Aweb][sp]=(str(indeg),str(outdeg))
      TLdict[int(S)][C][basal_p][Aweb][sp]=(str(STL),str(PATL))

  f.close()

  return TLdict,degdict,persdict


# partic_dict takes s, c, network, species, size, motif: count
def read_participation_file(pardir):
  par_dict={}
  full_motifs={3:set(),4:set()}
  for s in range(50,110,10):
    par_dict[s]={}
    for c in os.listdir(pardir+'/'+str(s)):
      par_dict[s][c]={}
      for netfile in os.listdir(pardir+'/'+str(s)+'/'+str(c)+'/'):
        i=netfile.split('.participation')[0].split('_')[2]
        if i not in par_dict[s][c]:
          par_dict[s][c][i]={}
        for sp in range(1,s+1):
          if '"sp'+str(sp)+'"' not in par_dict[s][c][i]:
            par_dict[s][c][i]['"sp'+str(sp)+'"']={}
        net=pardir+'/'+str(s)+'/'+str(c)+'/'+netfile
        size=int(netfile.split('_')[-1])
        f=open(net,'r')
        fixed_header=['motif']
        for line in f:
          if len(line.split())>1 and 'X' not in line.split()[0]:
            for mot in line.split():
              fixed_header.append(mot)
              full_motifs[size].add(mot)
          elif len(line.split())>1 and 'X' in line.split()[0]:  
            spnm='"sp'+line.split()[0][2:-1]+'"'
            par_dict[s][c][i][spnm][size]={}
            for j in range(1,len(fixed_header)):
              par_dict[s][c][i][spnm][size][fixed_header[j]]=line.split()[j]  
        f.close()
  return par_dict, full_motifs


# Creates a dictionary with position counts listed by species, connectance, network, species.
def read_role_file(roledir):
  role_dict={}
  full_positions=set()
  for s in range(50,110,10):
    role_dict[s]={}
    for c in os.listdir(roledir+str(s)):
      role_dict[s][c]={}
      for netfile in os.listdir(roledir+str(s)+'/'+str(c)+'/'):
        i=netfile.split('.roles')[0].split('_')[2]
        if '.roles_3' in netfile:
          role_dict[s][c][i]={}
          net=roledir+'/'+str(s)+'/'+str(c)+'/'+netfile
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
                sp='"sp'+line.split()[0].split('X')[1]
                role_dict[s][c][i][sp]={}
                for j in range(1,len(fixed_header)):
                  role_dict[s][c][i][sp][fixed_header[j]]=line.split()[j]

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
      for s in partic_dict:
        for c in partic_dict[s]:
          for net in partic_dict[s][c]:
            for sp in partic_dict[s][c][net]:
              motif_counts[motif].add(partic_dict[s][c][net][sp][size][motif])
      if motif_counts[motif]!=set(['0']):
        good_motifs[size].add(motif)
  return good_motifs


def sort_nonzero_roles(role_dict,full_positions):
  good_motifs=set()
  position_counts={}
  for motif in full_positions:
    position_counts[motif]=set()
    for s in role_dict:
      for c in role_dict[s]:
        for net in role_dict[s][c]:
          for sp in role_dict[s][c][net]:
            position_counts[motif].add(role_dict[s][c][net][sp][motif])
    if position_counts[motif]!=set(['0']):
      good_motifs.add(motif)
  return good_motifs


def collect_and_print_3sp(outfile,role_dict,good_positions,partic_dict,good_motifs,degdict,TLdict,persdict):
  f=open(outfile,'w')
  f.write('Size\tConnectance\tBasal_p\tNetwork\tSpecies\tPersistence\tin_Degree\tout_Degree\tSTL\tPATL\tm')
  f.write('\tm'.join(sorted(good_motifs[3]))+'\tp')
  f.write('\tp'.join(sorted(good_positions))+'\n')
  for s in range(50,110,10):
    for c in ['0.02','0.06','0.1','0.14','0.18']:
      for bp in [0.0,0.2,0.4,0.6,0.8,1.0]:
        for i in persdict[s][c][bp]:
          for sp in persdict[s][c][bp][i]: # Seems to be only non-basal species?
            f.write('\t'.join([str(s),c,str(bp),str(i),sp]))
            f.write('\t'+persdict[s][c][bp][i][sp])
            f.write('\t'+'\t'.join(degdict[s][c][bp][i][sp]))
            f.write('\t'+'\t'.join(TLdict[s][c][bp][i][sp]))
            for motif in sorted(good_motifs[3]):
              f.write('\t'+str(partic_dict[s][c][i][sp][3][motif]))
            for pos in sorted(good_positions):
              f.write('\t'+str(role_dict[s][c][i][sp][pos]))
            f.write('\n')
  f.close()

  return


def collect_and_print_4sp(outfile,partic_dict,good_motifs,degdict,TLdict,persdict):
  f=open(outfile,'w')
  f.write('Size\tConnectance\tBasal_p\tNetwork\tSpecies\tPersistence\tin_Degree\tout_Degree\tSTL\tPATL\tm')
  f.write('\tm'.join(sorted(good_motifs[4]))+'\n')
  for s in range(50,110,10):
    for c in ['0.02','0.06','0.1','0.14','0.18']:
      for bp in [0.0,0.2,0.4,0.6,0.8,1.0]:
        for i in persdict[s][c][bp]:
          for sp in persdict[s][c][bp][i]:
            f.write('\t'.join([str(s),c,str(bp),str(i),sp]))
            f.write('\t'+persdict[s][c][bp][i][sp])
            f.write('\t'+'\t'.join(degdict[s][c][bp][i][sp]))
            f.write('\t'+'\t'.join(TLdict[s][c][bp][i][sp]))
            for motif in sorted(good_motifs[4]):
              f.write('\t'+str(partic_dict[s][c][i][sp][4][motif]))
            f.write('\n')
  f.close()

  return


def main():
  # role_dict takes s, c, network, species, position: count
  # full_positions is a full set of positions
  role_dict,full_positions=read_role_file('../data/roles/empirical/')
  good_positions=sort_nonzero_roles(role_dict,full_positions)
  print 'roles read'
  # partic_dict takes s, c, network, species, size, motif: count
  # full_motifs is a dict for size: motif
  partic_dict,full_motifs=read_participation_file('../data/motif_participation/empirical/')
  good_motifs=sort_nonzero_positions(partic_dict,full_motifs)
  print 'participation read'
  # TLdict,degdict=read_TL_file('../data/TL/',role_dict) # Role dict to give a reference of i's
  # print 'TLs and degrees read'


  linear_TLdict,linear_degdict,linear_perdict=read_persistence_file('../data/all_empirical.tsv',role_dict)
  print 'persistence read'
  collect_and_print_3sp('../data/empirical_3sp_roles_participation.tsv',role_dict,good_positions,partic_dict,good_motifs,linear_degdict,linear_TLdict,linear_perdict)
  print '3sp role files complete'
      

if __name__ == '__main__':
  main()
