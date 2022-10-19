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


def get_complement(seq):
    new_seq = ''
    complement = {'A':'T','T':'A','C':'G','G':'C'}
    for i in seq:
        new_seq+=complement[i]
    return new_seq


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


dj_n_mh_exact = np.zeros(8)
dj_n_mh_mismatch = np.zeros(8)
dj_nb_junctions = 0


for thisj_trim in tqdm(range(20)):
    for thisd_trim in range(10):
        for thisj in msa_jgenes['name']:
            for thisd in msa_dgenes['name']:

                thisj_seq = msa_jgenes[msa_jgenes['name']==thisj]['original_seq'].values[0]
                thisd_seq = msa_dgenes[msa_dgenes['name']==thisd]['original_seq'].values[0]
                ### subseq is the trimmed off piece of the gene
                subseq = thisj_seq[:thisj_trim]


                ### we keep only the gene until the trim
                thisd_seq = thisd_seq[:-thisd_trim]

                em = int(exact_match(thisd_seq, get_complement(subseq)))
                em = min(em,7)
                dj_n_mh_exact[em]+=1
                mm = int(mismatch_mh(thisd_seq, get_complement(subseq)))
                mm = min(mm,7)
                dj_n_mh_mismatch[mm]+=1
                dj_nb_junctions+=1


vd_n_mh_exact = np.zeros(8)
vd_n_mh_mismatch = np.zeros(8)
vd_nb_junctions = 0

for thisv_trim in tqdm(range(20)):
    for thisd_trim in range(10):

        for thisv in msa_vgenes['name']:
            for thisd in msa_dgenes['name']:

                thisv_seq = msa_vgenes[msa_vgenes['name']==thisv]['original_seq'].values[0]
                thisd_seq = msa_dgenes[msa_dgenes['name']==thisd]['original_seq'].values[0]

                subseq = thisd_seq[:thisd_trim] #d5' sequence that is trimmed off

                thisv_seq = thisv_seq[:-thisv_trim]

                em = int(exact_match(thisv_seq, get_complement(subseq)))
                em = min(em,7)
                vd_n_mh_exact[em]+=1
                mm = int(mismatch_mh(thisv_seq, get_complement(subseq)))
                mm = min(mm,7)
                vd_n_mh_mismatch[mm]+=1
                vd_nb_junctions+=1


dj_n_mh_mismatch/=dj_nb_junctions
dj_n_mh_exact/=dj_nb_junctions


vd_n_mh_mismatch/=vd_nb_junctions
vd_n_mh_exact/=vd_nb_junctions


np.save('MH_all_DJ_junctions_mismatch.npy',dj_n_mh_mismatch)
np.save('MH_all_DJ_junctions_exact.npy',dj_n_mh_exact)
np.save('MH_all_VD_junctions_mismatch.npy',vd_n_mh_mismatch)
np.save('MH_all_VD_junctions_exact.npy',vd_n_mh_exact)




