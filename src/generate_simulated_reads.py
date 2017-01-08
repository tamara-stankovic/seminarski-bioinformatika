import pandas as pd
import numpy as np
import random


def cut_ref_genome(fasta, exons):
	tmp = exons.groupby('gene_id').exon_len.sum() - 100 + 1
	prob_samp = tmp.reset_index()
	prob_samp['prob'] = prob_samp.exon_len / prob_samp.exon_len.sum()
	prob_samp['cum_prob'] = prob_samp.prob.cumsum()

	ex_pairs = exons.apply(lambda row: (row['start'], row['end']), axis=1).values
	rez = ''
	index_list = []
	for x in ex_pairs:
	    rez = rez + fasta[x[0]:x[1]+1]
	    index_list.append([i for i in range(x[0], x[1] + 1)])
	index_list = np.concatenate(index_list)

	return (rez, index_list)
