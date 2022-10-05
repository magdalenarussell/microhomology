import sys
import pandas as pd
import numpy as np
from tqdm import tqdm


def exact_match(seq5,res_seq3):
    min_seq = min(len(seq5), len(res_seq3))
    match = 0
    for i in range(1,min_seq):
        if seq5[-i]==res_seq3[-i]:
            match +=1
        else:
            break
    return match

def get_p_nuc_3prime(seq):
    complement = {'A':'T','T':'A','C':'G','G':'C'}
    return seq+complement[seq[-1]]+complement[seq[-2]]

def get_p_nuc_5prime(seq):
    complement = {'A':'T','T':'A','C':'G','G':'C'}
    return complement[seq[1]]+complement[seq[0]]+seq


def get_complement(seq):
    new_seq = ''
    complement = {'A':'T','T':'A','C':'G','G':'C'}
    for i in seq:
        new_seq+=complement[i]
    return new_seq


def exact_mh(s1, s2):
    fit = [0]
    for n in range(1,min(len(s1), len(s2))+1):
        same = np.sum([i==j for i,j in zip(s1[-n:],s2[:n])])
        if same == n:
            fit.append(same)
        else:
            fit.append(0)
        #print (s1[-n:], s2[:n])
    return np.max(fit)

def mismatch_mh(s1, s2):
    fit = [0]
    mismatches = [0]
    fit = np.sum([i==j for i,j in zip(s1,s2)])

    return fit

muscle = pd.read_csv('muscle_clustal_V_pastanchor_with_Pnuc.clw')
vgen = []
seqs = []
for i in range(muscle.shape[0]):
    cells = list(muscle.loc[i])[0].split(' ')
    for val in cells:
        if not val == '':
            if 'TRB' in val:
                vgen.append(val)
            else:
                seqs.append(val)
msa_vgenes = pd.DataFrame([vgen,seqs]).T
msa_vgenes.columns = ['name','value']

msa_vgenes['original_seq'] = [''.join(i.split('-'))[:-1] for i in msa_vgenes['value']]



muscle = pd.read_csv('muscle_clustal_J_withPnuc.clw')
vgen = []
seqs = []
for i in range(muscle.shape[0]):
    cells = list(muscle.loc[i])[0].split(' ')
    for val in cells:
        if not val == '':
            if 'TRB' in val:
                vgen.append(val)
            else:
                seqs.append(val)
msa_jgenes = pd.DataFrame([vgen,seqs]).T
msa_jgenes.columns = ['name','value']
msa_jgenes = msa_jgenes.iloc[:15]

msa_jgenes['original_seq'] = [''.join(i.split('-'))[:-1] for i in msa_jgenes['value']]


muscle = pd.read_csv('muscle_clustal_D_pastanchor.clw')
vgen = []
seqs = []
for i in range(muscle.shape[0]):
    cells = list(muscle.loc[i])[0].split(' ')
    for val in cells:
        if not val == '':
            if 'TRB' in val:
                vgen.append(val)
            else:
                seqs.append(val)
msa_dgenes = pd.DataFrame([vgen,seqs]).T
msa_dgenes.columns = ['name','value']
msa_dgenes = msa_dgenes.iloc[:3]
msa_dgenes['original_seq'] = [''.join(i.split('-'))[:-1] for i in msa_dgenes['value']]


pt_name = sys.argv[1]
folder = '/fh/fast/matsen_e/shared/tcr-gwas/emerson_stats/'
data = pd.read_csv(f'{folder}trim_stats_{pt_name}_w_pnucs_and_ties.tsv',sep='\t')


n_mh_exact = []
n_mh_mismatch = []
mismatch_numbers = []

for i in tqdm(range(data.shape[0])):
    thisj = data.iloc[i]['j_gene']
    thisd = data.iloc[i]['d_gene']
    if thisj in list(msa_jgenes['name']) and thisd in list(msa_dgenes['name']):
        thisj_seq = msa_jgenes[msa_jgenes['name']==thisj]['original_seq'].values[0]
        #print (thisj_seq)
        nb_p = max(data.iloc[i]['j_pnucs'],0)
        thisj_trim = data.iloc[i]['j_trim']+2-nb_p


        ### subseq is the trimmed off piece of the gene
        subseq = thisj_seq[:thisj_trim]


        thisd_seq = msa_dgenes[msa_dgenes['name']==thisd]['original_seq'].values[0]
        #print (thisd_seq)
        nb_p = max(data.iloc[i]['d1_pnucs'],0)
        thisd_trim = data.iloc[i]['d1_trim']+2-nb_p

        ### we keep only the gene until the trim
        thisd_seq = thisd_seq[:-thisd_trim]

        n_mh_exact.append(exact_match(thisd_seq, get_complement(subseq)))
        n_mh_mismatch.append(mismatch_mh(thisd_seq, get_complement(subseq)))
    else:
        n_mh_exact.append(0)
        n_mh_mismatch.append(0)


data['exact_mh_dj'] = n_mh_exact
data['mismatch_mh_dj'] = n_mh_mismatch





n_mh_exact = []
n_mh_mismatch = []
mismatch_numbers = []

for i in tqdm(range(data.shape[0])):
    thisv = data.iloc[i]['v_gene']
    thisd = data.iloc[i]['d_gene']
    if thisv in list(msa_vgenes['name']) and thisd in list(msa_dgenes['name']):
        thisv_seq = msa_vgenes[msa_vgenes['name']==thisv]['original_seq'].values[0]
        #print (thisj_seq)
        nb_p = max(data.iloc[i]['v_pnucs'],0)
        thisv_trim = data.iloc[i]['v_trim']+2-nb_p


        thisd_seq = msa_dgenes[msa_dgenes['name']==thisd]['original_seq'].values[0]
        #print (thisd_seq)
        nb_p = max(data.iloc[i]['d0_pnucs'],0)
        thisd_trim = data.iloc[i]['d0_trim']+2-nb_p

        #subseq = get_p_nuc_5prime(thisj_seq)
        subseq = thisd_seq[:thisd_trim] #d5' sequence that is trimmed off


        #thisd_seq = get_p_nuc_3prime(thisd_seq)
        thisv_seq = thisv_seq[:-thisv_trim]

        n_mh_exact.append(exact_match(thisv_seq, get_complement(subseq)))
        mh  = mismatch_mh(thisv_seq, get_complement(subseq))
        n_mh_mismatch.append(mh)
    else:
        n_mh_exact.append(0)
        n_mh_mismatch.append(0)


data['exact_mh_vd'] = n_mh_exact
data['mismatch_mh_vd'] = n_mh_mismatch


data.to_csv(f'processed_both/{pt_name}_withmh.csv')











