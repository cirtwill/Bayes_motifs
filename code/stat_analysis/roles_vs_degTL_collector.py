import sys
import os
import decimal


def permfile_reader(permfile):
	temp={}
	f=open(permfile,'r')
	for line in f:
		if line.split()[0]=='"1"':
			output=float(line.split()[1])
	f.close()
	return output


def pval_calculator(model_F,perms):
	greater=[x for x in perms if x>=model_F]

	p=float(len(greater))/float(len(perms))

	return p


def recreate_aovtab(datadict,outfile):

	f=open(outfile,'w')
	f.write('S\tC\text_corr\tdf\tSS\tMS\tModel_F\tR2\tpval\n')
	for C in sorted(datadict):
		dats=datadict[C]
		f.write(str(dats['S'])+'\t'+str(dats['C'])+'\t'+str(dats['ext_corr'])+'\t')
		f.write(str(dats['df'])+'\t'+str(dats['SS'])+'\t'+str(dats['MS'])+'\t')
		f.write(str(dats['model_F'])+'\t'+str(dats['R2'])+'\t'+str(dats['pval'])+'\n')
	f.close()


def betadisper_reader(dispfile):
	f=open(dispfile,'r')
	for line in f:
		if line.split()[0]=='"Groups"':
			df=int(line.split()[1])
			SSE=float(line.split()[2])
			MSE=float(line.split()[3])
			F=float(line.split()[4])
			p=float(line.split()[5])
	f.close()

	result={'F':F,'p':p}
	return result


def tukey_reader(tukeyfile):
	results={}
	f=open(tukeyfile,'r')
	for line in f:
		if line.split()[0]!='"diff"':
			groups=line.split()[0][1:-1]
			diff=float(line.split()[1])
			lower=float(line.split()[2])
			upper=float(line.split()[3])
			p=float(line.split()[4])
			results[groups]=((diff,p))
	f.close()

	return results


def tukey_outfile(dicto,style):
	f=open('../../data/summaries/'+style+'_Tukey.tsv','w')
	f.write('Species\tC\tGroup1\tGroup2\tDiff\tP\n')

	for s in dicto:
		for c in dicto[s]:
			for group in sorted(dicto[s][c]['Tukey'].keys()):
				f.write(str(s)+'\t'+str(c)+'\t'+'\t'.join(group.split('-')))
				(diff,p)=dicto[s][c]['Tukey'][group]
				f.write('\t'+str(diff)+'\t'+str(p)+'\n')	
	f.close()
	print style, ' Tukey outfile done'


def betadisp_outfile(dicto,style):

	f=open('../../data/summaries/'+style+'_betadisper.tsv','w')
	f.write('Species\tC\tF\tP\n')
	for s in dicto:
		for c in dicto[s]:
			f.write(str(s)+'\t'+str(c)+'\t')
			val=dicto[s][c]['betadisper']
			f.write(str(val['F'])+'\t'+str(val['p'])+'\n')
	f.close()
	print style, ' betadisper outfile done'


def permanova_outfile(dicto,style):

	f=open('../../data/summaries/'+style+'_permanova.tsv','w')
	f.write('Species\tC\tP\n')
	for s in dicto:
		for c in dicto[s]:
			f.write(str(s)+'\t'+str(c)+'\t')
			f.write(str(dicto[s][c]['pval'])+'\n')
	f.close()
	print style, ' permanova outfile done'


def main():

	# Things I have: Pvals from observed and randomized anovas
	# betadisper anova table
	# betadisper Tukey test

	degdir='../../data/perms_deg/'
	TLdir='../../data/perms_STL/'

	TLdict={}
	degdict={}

	for s in sorted(os.listdir(degdir)):
		TL_outfile='../../data/summaries/'+s+'/roles3sp_vs_STL_permanova_'+s+'.tsv'
		deg_outfile='../../data/summaries/'+s+'/roles3sp_vs_deg_permanova_'+s+'.tsv'
		TLdict[s]={}
		degdict[s]={}
		for c in [x for x in sorted(os.listdir(degdir+s))]:
			TLlist=sorted(os.listdir(TLdir+s+'/'+c+'/perms/'))
			deglist=sorted(os.listdir(degdir+s+'/'+c+'/perms/'))

			TLdict[s][c]={'model_F':permfile_reader(TLdir+s+'/'+c+'/STL_vs_3sp_obs.tsv')}
			degdict[s][c]={'model_F':permfile_reader(degdir+s+'/'+c+'/deg_vs_3sp_obs.tsv')}

			all_TL_perms=[]
			all_deg_perms=[]

			for TLperm in TLlist:
				all_TL_perms.append(permfile_reader(TLdir+s+'/'+c+'/perms/'+TLperm))

			for degperm in deglist:
				all_deg_perms.append(permfile_reader(degdir+s+'/'+c+'/perms/'+degperm))


			if len(all_TL_perms)==999 and len(all_deg_perms)==999:
				TLdict[s][c]['pval']=pval_calculator(TLdict[s][c]['model_F'],all_TL_perms)
				degdict[s][c]['pval']=pval_calculator(degdict[s][c]['model_F'],all_deg_perms)
				print 'All present and correct with ',s,c
			else:
				print 'Missing ',str(999-len(all_TL_perms)),' perms for run: ',s,' ',c,str(999-len(all_deg_perms)),' perms for run: ',s,' ',c

			TLdict[s][c]['betadisper']=betadisper_reader(TLdir+s+'/'+c+'/STL_vs_3sp_disp_anova.tsv')
			TLdict[s][c]['Tukey']=tukey_reader(TLdir+s+'/'+c+'/STL_vs_3sp_disp_Tukey.tsv')
			degdict[s][c]['betadisper']=betadisper_reader(degdir+s+'/'+c+'/deg_vs_3sp_disp_anova.tsv')
			degdict[s][c]['Tukey']=tukey_reader(degdir+s+'/'+c+'/deg_vs_3sp_disp_Tukey.tsv')

	# Outfiles for Tukeys
	tukey_outfile(TLdict,'TL')
	tukey_outfile(degdict,'deg')

	# Outfiles for betadispers
	betadisp_outfile(TLdict,'TL')
	betadisp_outfile(degdict,'deg')

	# Outfiles for permanovas
	permanova_outfile(TLdict,'TL')
	permanova_outfile(degdict,'deg')


if __name__ == '__main__':
  main()
