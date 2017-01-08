import random
import numpy as np
import pysam

def picker(index_list, genome):
    n = len(genome) - 150
    seed = random.randint(0, n)
    
    return {'coordinate': index_list[seed], 'simulated_read': genome[seed:seed+151]}

def generateRandomSequence(n):
    l = ['T', 'C', 'G', 'A']
    return ''.join([random.choice(l) for x in range(n)])

def generateDel(sequence, n):
    if n > 50: raise Exception('len must be less than 50.')
    
    c = random.randint(0, 100)
    e = c + n
    return sequence[:c] + sequence[c + n:c+n+100-c]

def generateIns(seqence, n):
    if n > 50: raise Exception('len must be less than 50.')
        
    c = random.randint(0, 100)
    beg = seqence[:c]
    ins = generateRandomSequence(n)
    if(c + n < 100):
        return beg + ins + seqence[c + n:c+n+(100-c-n)]
    else:
        return beg + ins[:100-c-n]

def generateSample(n, indel_prob, del_prob, index_list, reference):
    rez_list = []
    for x in range(n):
        l = random.random()
        tmp = picker(index_list, reference)
        if l < indel_prob:
            l2 = random.random()
            if l2 < del_prob:
                read = {'read': generateDel(tmp['simulated_read'], random.randint(5, 10)),
                        'del': True,
                        'name': 'simulated_read_{:d}:{:d}'.format(x, tmp['coordinate']),
                        'coordinate': tmp['coordinate']
                        }
            else:
                read = {'read': generateIns(tmp['simulated_read'], random.randint(5, 10)),
                        'ins': True,
                        'name': 'simulated_read_{:d}:{:d}'.format(x, tmp['coordinate']),
                        'coordinate': tmp['coordinate']
                        }
                
        else:
            read = {'read': tmp['simulated_read'],
                    'name': 'simulated_read_{:d}:{:d}'.format(x, tmp['coordinate']),
                    'coordinate': tmp['coordinate']
                    }           
        rez_list.append(read)
    return rez_list
    
def generateFasta(seq_list, out_path):
    l = ['>' + x['name'] + '\n' + x['read'] + '\n' for x in seq_list]
    with open(out_path, 'w') as f:
        for y in l:
            f.write(y)

def false_align_precentage(bam_path, seq_len):
	sam = pysam.AlignmentFile(bam_path, 'r')

	c = 0	
	for x in sam:
	    t = str(x).split("\t")
	    name = t[0]
	    mapped_to = t[3]
	    real_coord = name.split(":")[-1]

	    if (int(real_coord) - int(mapped_to) != 0):
	    	c = c + 1
	return (c / seq_len) * 100
	    

        