import pandas as pd
import numpy as np

def load_data():
	gtf = pd.read_table('data/chr20.gtf', header=None)
	gtf.columns = ['chr', 'func', 'type', 'start', 'end', 'irr_1', 'irr_2', 'irr_3', 'add_info']

	exons = gtf[gtf['type'] == 'exon'].copy()
	exons['gene_id'] = exons['add_info'].apply(lambda x: x.split(';')[0])
	exons['exon_len'] = exons['end'] - exons['start'] + 1
	f = open ('data/chr20.fa')
	next (f)
	fasta = "".join([x.strip() for x in f.readlines()])

	return (exons,fasta)