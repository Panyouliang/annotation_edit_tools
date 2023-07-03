import sys
import numpy as np



plus_gene = {}
minus_gene = {}
with open(sys.argv[1],'r') as f:
    for line in f:
        if line.strip():
            line = line.strip().split('\t')
            chrs = line[0]
            if line[2] == 'mRNA':
                pos = [int(line[3]),int(line[4])]
                if line[6] == '+':
                    if chrs in plus_gene.keys():
                        plus_gene[chrs].append(pos)
                    else:
                        plus_gene[chrs] = [pos]
                if line[6] == '-':
                    if chrs in minus_gene.keys():
                        minus_gene[chrs].append(pos)
                    else:
                        minus_gene[chrs] = [pos]

plus_gene = {k: sorted(v) for k,v in plus_gene.items()}
minus_gene = {k: sorted(v) for k,v in minus_gene.items()}


def IGRs(dicts,lis):
    for k,v in dicts.items():
        if len(v) > 1:
            for i in range(len(v)-1):
                length = abs(int(v[i][1])-int(v[i+1][0]))

                lis.append(length)
    return lis

lis = []
lis = IGRs(plus_gene,lis)
lis = IGRs(minus_gene,lis)


lis = sorted(lis)
lis = np.array(lis)

p25 = np.percentile(lis,25)

media = np.median(lis)
mean = np.mean(lis)
maxi = np.max(lis)
mini = np.min(lis)

print('min: '+str(mini),'p25: '+str(p25),'mean: '+str(mean),'media: '+str(media),'max: '+str(maxi),sep='\t')
