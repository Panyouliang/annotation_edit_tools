# !/usr/bin/env python3
# -*- coding:utf-8 -*-

import sys,gzip,argparse
import numpy as np

#================================================================================
version = "v1.0"
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''
==================================================================
This is a script for evaluating the quality of UTRs infomations
of genome annotation.

Author: Panyouliang, panyouliang@genomics.cn
Version: v1.0
Date: 2023-06-13, yyyy-mm-dd
==================================================================''')

parser.add_argument('-v', '--version', action='version', version=version)
parser.add_argument('-gff', metavar='annotation file', type=str, required=True, help='Please input the annotation file(None isoform)')
parser.add_argument('-sp', metavar='Specie name', type=str, required=False, help='Please input the Specie name')
args = parser.parse_args()

#=================================================================================


def read_gff(F):

    gene_dict, UTR3, UTR5, coding = {}, 0, 0, {}

    if F[-3:] == '.gz':
        f = gzip.open(F,'rt')
    else:
        f = open(F,'r')
    utr5 = ['five_prime_UTR','five_prime_utr','UTR5']
    utr3 = ['three_prime_UTR','three_prime_utr','UTR3']

    for line in f:
        if line.strip():
            if line[0] != '#':
                line = line.strip().split('\t')
                if line[2] == 'mRNA':
                    gidx = line[8].find('ID')
                    gid = line[8][gidx:].split(';')[0].split('=')[1] +' '+line[6]
                    gene_dict[gid] = {}

                    coding[gid] = []

                else:
                    if line[2] in utr5:
                        if 'UTR5' in gene_dict[gid].keys():
                            gene_dict[gid]['UTR5'].append(abs(int(line[3])-int(line[4]))+1)
                        else:
                            gene_dict[gid]['UTR5'] = [abs(int(line[3])-int(line[4]))+1]

                    if line[2] in utr3:
                        if 'UTR3' in gene_dict[gid].keys():
                            gene_dict[gid]['UTR3'].append(abs(int(line[3])-int(line[4]))+1)
                        else:
                            gene_dict[gid]['UTR3'] = [abs(int(line[3])-int(line[4]))+1]
                    if line[2] == 'CDS':
                        pidx = line[8].find('Parent')
                        pid = line[8][pidx:].split(';')[0].split('=')[1] +' '+line[6]
                        length = abs(int(line[4])-int(line[3])) + 1
                        if pid in coding.keys():
                            coding[gid].append(length)

    coding_len = sum(sum(v) for v in coding.values())
    coding_len_per_gene = coding_len/len(coding.keys())


    for k,v in gene_dict.items():
        if len(v) > 0:
            for m,n in v.items():
                v[m] = sum(n)

    return gene_dict,coding_len_per_gene

def main():
    
    gff,sp = args.gff,args.sp

    gene_dict,coding_len_per_gene = read_gff(args.gff)

    gene_num,intact_UTR,UTR5,UTR3,UTR5_L,UTR3_L = 0,0,0,0,0,0
    UTR3_list,UTR5_list = [],[]


    for k,v in gene_dict.items():
        gene_num += 1
        if 'UTR3' in v.keys() and 'UTR5' in v.keys():
            intact_UTR += 1

        if 'UTR5' in v.keys():
            UTR5 += 1
            UTR5_L += v['UTR5']
            UTR5_list.append(v['UTR5'])

        if 'UTR3' in v.keys():
            UTR3 += 1
            UTR3_L += v['UTR3']
            UTR3_list.append(v['UTR3'])

    UTR3_list_np,UTR5_list_np = np.array(UTR3_list),np.array(UTR5_list)

    media_UTR3,media_UTR5 = np.median(UTR3_list_np),np.median(UTR5_list_np)

    if UTR5 > 0 or UTR3 > 0:
        #print('species','total genes','UTR5+UTR3 genes','UTR5+UTR3 proportion','UTR5 genes','UTR5 proportion','UTR3 genes','UTR3 proportion','coding len per gene','UTR5 len avg (bp)', 'UTR5 len median','UTR3 len avg (bp)', 'UTR3 len median','UTR5 total len (bp)','UTR3 total len (bp)',sep='\t')

        print(sp,gene_num,intact_UTR,format(intact_UTR/gene_num,'.2f'),UTR5,format(UTR5/gene_num,'.2f'),UTR3,format(UTR3/gene_num,'.2f'),int(coding_len_per_gene),int(sum(UTR5_list)/len(UTR5_list)),media_UTR5,int(sum(UTR3_list)/len(UTR3_list)),media_UTR3,UTR5_L,UTR3_L,sep='\t')
    else:
        print(sp,gene_num,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,sep='\t')


if __name__ == '__main__':
    main()





