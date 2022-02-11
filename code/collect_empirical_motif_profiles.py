import os
import sys

# Creates a dictionary with motif counts listed by species, connectance, network, motif.
def read_role_file(roledir):
  role_dict={}
  for webfile in os.listdir(roledir):
    web=webfile.split('.')[0]
    role_dict[web]={}
    counts=[]
    if '.motifs_3' in webfile:
      net=roledir+webfile
      f=open(net,'r')
      for line in f:
        if len(line.split())>1 and line.split()[0]!='motif':
          motif=str(line.split()[0])
          count=int(line.split()[1])
          counts.append(count)
          if motif in ['6','12','36','38']:
            role_dict[web][motif]=count
      f.close()
    role_dict[web]['total']=sum(counts)
  return role_dict


def write_motif_file(role_dict,outfile):
  f=open(outfile,'w')
  f.write('Size\tConnectance\tNetwork\tm6\tm12\tm36\tm38\tTotal\n')
  for web in role_dict:
    f.write('\t'.join([str(Cdict[web][0]),str(Cdict[web][1]),web]))
    f.write('\t'+str(role_dict[web]['6']))
    f.write('\t'+str(role_dict[web]['12']))
    f.write('\t'+str(role_dict[web]['36']))
    f.write('\t'+str(role_dict[web]['38']))
    f.write('\t'+str(role_dict[web]['total'])+'\n')
  f.close()

  return

# Calculated S and C=L/(S*(S-1)) for empirical webs
Cdict={'kongsfjorden':(260,0.024),'loughhyne':(341,0.043),'reef_noSynodus':(242,0.056),
    'stmarks':(141,0.089),'weddell':(490,0.060),'ythanjacob':(88,0.051)}

def main():
  # role_dict takes s, c, network, motif: count
  # only motifs 6, 12, 36, 38 exist in the Bayesian networks
  motif_dict=read_role_file('../data/motif_frequencies/empirical/')

  write_motif_file(motif_dict,'../data/3sp_empirical_motif_profiles.tsv')
  print 'Motif roles read'


if __name__ == '__main__':
  main()