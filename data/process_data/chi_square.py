# exec(open("chi_square.py").read())

from scipy import stats
import pandas as pd
import glob
import os
import scipy
from statsmodels.stats.multitest import multipletests


df = pd.read_csv("design.txt",sep="\t",header=None)
# print (df)
sample_list = df[1].unique().tolist()

g1 = "Treated_mRNA"
g2 = "Treated_protein"
control = "Control"



def parse_df(f):
	# print (f)
	# 
	df = pd.read_csv(f,sep="\t",index_col=0)
	label=f.split("_")[0]
	# if %Reads is 0, we won't get to know the total reads, then just set it as -1
	# print (df.sort_values('%Reads').head(10))
	df['%Reads'] = df['%Reads'].replace(0,-1)
	# print (df.sort_values('%Reads').head(10))
	# print (df.sort_values('%Reads').tail(10))
	df['%s_Total'%(label)] = df["#Reads"]/(df['%Reads']/100)
	
	df['%s_Total'%(label)] = df['%s_Total'%(label)].astype(int)
	df['%s_Reads'%(label)] = df["#Reads"]
	df = df[['%s_Total'%(label),'%s_Reads'%(label)]]
	try:
		df.loc[df['%s_Reads'%(label)]==-1]=[0,0]
	except:
		df.loc[df['%s_Reads'%(label)]==-1]=0
	return df


def parse_list(myList,label):
	# print (myList)
	df_list = [parse_df(glob.glob("%s*.allele.edit.tsv"%(f))[0]) for f in myList]
	# print (df_list)
	df = pd.concat(df_list,axis=1)
	
	df['%s_Total'%(label)]=df[['%s_Total'%(f) for f in myList]].sum(axis=1)
	df['%s_Reads'%(label)]=df[['%s_Reads'%(f) for f in myList]].sum(axis=1)
	df = df[['%s_Total'%(label),'%s_Reads'%(label)]]
	# print (df.head())
	return df

def fisher_exact(r,label):
	control = "Control"
	try:
		_,p,_,_ = stats.chi2_contingency([
		[r["%s_Reads"%(label)],r["%s_Total"%(label)]-r["%s_Reads"%(label)]],
		[r["%s_Reads"%(control)],r["%s_Total"%(control)]-r["%s_Reads"%(control)]]
								])
	except:
		# print (r)
		return 1
	# print (p)
	return p



def get_diff(r,label):
	control = "Control"
	treatment = float(r["%s_Reads"%(label)])/r["%s_Total"%(label)]
	control = float(r["%s_Reads"%(control)])/r["%s_Total"%(control)]
	return treatment - control

# chi2, p, dof, expected = sp.stats.chi2_contingency(x)
df_results = []
# print (df)
for s,d in df.groupby(1):
	print (s)
	df_list = []
	for s2,d2 in d.groupby(3):
		print (s2)
		# print (d2)
		# continue
		
		a = d2[d2[0]==g1]
		a.index = a[1].tolist()
		
		c = d2[d2[0]==control]
		c.index = c[1].tolist()

		
		b = d2[d2[0]==g2]
		b.index = b[1].tolist()

		
		df1 = parse_list(a[2].tolist(),g1)
		df2 = parse_list(b[2].tolist(),g2)
		df3 = parse_list(c[2].tolist(),control)
		current_df = pd.concat([df1,df2,df3],axis=1)
		

		current_df = current_df[current_df['Control_Total']!=0]

		current_df = current_df.sort_values("Control_Total")

		
		current_df['%s_vs_control_p_value'%(g1)] = current_df.apply(lambda r:fisher_exact(r,g1),axis=1)
		current_df['%s_vs_control_diff'%(g1)] = current_df.apply(lambda r:get_diff(r,g1),axis=1)
		current_df['%s_vs_control_p_value'%(g2)] = current_df.apply(lambda r:fisher_exact(r,g2),axis=1)
		current_df['%s_vs_control_diff'%(g2)] = current_df.apply(lambda r:get_diff(r,g2),axis=1)
		
		df_list.append(current_df)
	tmp = pd.concat(df_list)

	tmp['%s_vs_control_FDR'%(g1)] = multipletests(pvals=tmp['%s_vs_control_p_value'%(g1)].tolist(), alpha=0.05, method="fdr_bh")[1]
	tmp['%s_vs_control_FDR'%(g2)] = multipletests(pvals=tmp['%s_vs_control_p_value'%(g2)].tolist(), alpha=0.05, method="fdr_bh")[1]
	# print (tmp.columns)
	columns = ["Treated_mRNA_Total","Treated_mRNA_Reads","Treated_protein_Total","Treated_protein_Reads","Control_Total","Control_Reads","Treated_mRNA_vs_control_p_value","Treated_mRNA_vs_control_diff","Treated_mRNA_vs_control_FDR","Treated_protein_vs_control_p_value","Treated_protein_vs_control_diff","Treated_protein_vs_control_FDR"]
	tmp = tmp[columns]
	tmp['Control_%edit'] = tmp['Control_Reads']/tmp['Control_Total']
	tmp.columns = [s+"|"+x for x in tmp.columns]

	df_results.append(tmp)
	
df_final = pd.concat(df_results,axis=1)

FDR_columns = []
diff_columns = []

for c in df_final.columns:
	if not "mRNA" in c:
		continue
	if "FDR" in str(c):
		FDR_columns.append(c)
	if "diff" in str(c):
		diff_columns.append(c)

df_final['max_diff'] = df_final[diff_columns].max(axis=1)

def pass_cutoff(r):
	flag = False
	s_list = []
	for s in df[1].unique().tolist():
		FDR = "%s|Treated_mRNA_vs_control_FDR"%(s)
		diff = "%s|Treated_mRNA_vs_control_diff"%(s)
		if r[FDR] <0.05 and r[diff]>0.005:
			s_list.append(s)
			
	if len(s_list)==0:
		return "Not significant"
	else:
		return ",".join(s_list)
df_final['selected'] = df_final.apply(pass_cutoff,axis=1)

df_final = df_final.sort_values('max_diff',ascending=False)
# print (df_final.head())
df_final.to_csv("chi_square_test.csv")


	
