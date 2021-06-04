import os
import sys

# Creates a dictionary with motif counts listed by species, connectance, network, motif.
def read_role_file(roledir):
	role_dict={}
	for s in range(50,110,10):
		role_dict[s]={}
		for c in os.listdir(roledir+str(s)):
			role_dict[s][c]={}
			for netfile in os.listdir(roledir+str(s)+'/'+str(c)+'/'):
				i=netfile.split('.motifs')[0].split('_')[2]
				counts=[]
				if '.motifs_3' in netfile:
					role_dict[s][c][i]={}
					net=roledir+'/'+str(s)+'/'+str(c)+'/'+netfile
					f=open(net,'r')
					for line in f:
						if len(line.split())>1 and line.split()[0]!='motif':
							motif=str(line.split()[0])
							count=int(line.split()[1])
							counts.append(count)
							if motif in ['6','12','36','38']:
								role_dict[s][c][i][motif]=count
					f.close()
				role_dict[s][c][i]['total']=sum(counts)
	return role_dict


def write_motif_file(role_dict,outfile):
	f=open(outfile,'w')
	f.write('Size\tConnectance\tNetwork\tm6\tm12\tm36\tm38\tTotal\n')
	for s in range(50,110,10):
		for c in ['0.02','0.06','0.1','0.14','0.18']:
			for i in role_dict[s][c]:
				f.write('\t'.join([str(s),c,str(i)]))
				f.write('\t'+str(role_dict[s][c][i]['6']))
				f.write('\t'+str(role_dict[s][c][i]['12']))
				f.write('\t'+str(role_dict[s][c][i]['36']))
				f.write('\t'+str(role_dict[s][c][i]['38']))
				f.write('\t'+str(role_dict[s][c][i]['total'])+'\n')
	f.close()

	return


def main():
	# role_dict takes s, c, network, motif: count
	# only motifs 6, 12, 36, 38 exist in the Bayesian networks
	motif_dict=read_role_file('../data/motif_frequencies/acyclic/')

	write_motif_file(motif_dict,'../data/3sp_motif_profiles.tsv')
	print 'Motif roles read'

	motif_dict=read_role_file('../data/motif_frequencies/original/')

	write_motif_file(motif_dict,'../data/3sp_motif_profiles_cyclicwebs.tsv')


if __name__ == '__main__':
  main()
