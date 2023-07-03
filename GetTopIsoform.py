#! /usr/bin/env python3
# -*- coding:utf-8 -*-

import sys,gzip,argparse

#================================================================================
version = "v1.0"
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''
======================================================================
This is a script for extract top length isoform from annotation files
that including different isoforms infomation.

Author: Panyouliang, panyouliang@genomics.cn
Version: v1.0
Date: 2023-06-13, yyyy-mm-dd
======================================================================''')
parser.add_argument('-v', '--version', action='version', version=version)
parser.add_argument('-gff', metavar='annotation file', type=str, required=True, help='Please input the annotation file(including isoforms')
parser.add_argument('-sp', metavar='Specie name', type=str, required=True, help='Please input the Specie name')
args = parser.parse_args()
#================================================================================


def read(F):
    if F[-3:] == '.gz':
        f = gzip.open(F,'rt')
    else:
        f = open(F,'r')
    gff_file,all_gff_file,iso_len,gene,iso_num = {},{},{},{},0
    UTR_list = ['five_prime_UTR','three_prime_UTR','five_prime_utr','three_prime_utr','UTR3','UTR5']

    for line in f:
        if line.strip():
            if line[0] != '#':
                line = line.strip().split('\t')
                if line[2] == 'gene':

                    gidx = line[8].find('ID')
                    g_id = line[8][gidx:].split(';')[0].split('=')[1]+' '+line[6]
                    gff_file[g_id] = {}
                    all_gff_file[g_id] = {}
                    iso_len[g_id] = {}
                    gene[g_id] = line

                elif line[2] == 'mRNA': 

                    iso_num += 1
                    tidx = line[8].find('ID')
                    t_id = line[8][tidx:].split(';')[0].split('=')[1]

                    ggidx = line[8].find('Parent')
                    #ggidx = line[8].find('gene_id')
                    ggid = line[8][ggidx:].split(';')[0].split('=')[1]

                    if ggid not in all_gff_file.keys():
                        all_gff_file[ggid] = {t_id:[line]}
                        gff_file[ggid] = {t_id:[line]}
                        iso_len[ggid] = {t_id:0}

                        gene_lis = line[:]
                        gene_lis[2] = 'gene'
                        gene_lis[8] = 'ID='+ggid+';'
                        gene[ggid] = gene_lis
                    else:
                        all_gff_file[ggid][t_id] = [line]
                        gff_file[ggid][t_id] = [line]
                        iso_len[ggid][t_id] = 0

                    coding_length = 0

                elif line[2] == 'exon' or line[2] == 'CDS':
                #elif line[2] == 'CDS':
                    idx = line[8].find('Parent')
                    pid = line[8][idx:].split(';')[0].split('=')[1]
                    if pid in gff_file[ggid].keys():
                        gff_file[ggid][t_id].append(line)
                        all_gff_file[ggid][t_id].append(line)

                    if line[2] == 'CDS':
                        if pid in iso_len[ggid].keys():
                            coding_length += abs(int(line[4]) - int(line[3]))
                            iso_len[ggid][t_id] = coding_length

                elif line[2] in UTR_list:
                    uidx = line[8].find('Parent')
                    u_id = line[8][uidx:].split(';')[0].split('=')[1]

                    if u_id in all_gff_file[ggid].keys() or u_id == ggid:
                        all_gff_file[ggid][t_id].append(line)


    for k,v in all_gff_file.items():
        for m,n in v.items():
            if k.endswith('+'):
                model = sorted(n, key = (lambda x:int(x[3])),reverse=True)
                v[m] = model
            else:
                model = sorted(n[1:], key = (lambda x:int(x[3])),reverse=False)
                model.insert(0,n[0])
                v[m] = model


    return gff_file, all_gff_file, iso_len, gene


def main():

    gff, name = args.gff, args.sp

    gene_num,mRNA_num = 0,0

    gff_file,all_gff_file,iso_len,gene = read(gff)

    '''
    out1 = open(name+'.TopLength.CDS.gff','w')
    for k,v in iso_len.items():
        if len(v) > 0:
            tr_id=max(v,key=v.get)
            out1.write('\t'.join(x for x in gene[k])+'\n')
            for line in gff_file[k][tr_id]:
                out1.write('\t'.join(x for x in line)+'\n')
    '''

    out2 = open(name+'.TopLength.ALL.gff','w')
    for k,v in iso_len.items():
        if len(v) > 0:
            gene_num += 1
            mRNA_num += len(v)
            tr_id=max(v,key=v.get)
            out2.write('\t'.join(x for x in gene[k])+'\n')
            for line in all_gff_file[k][tr_id]:
                out2.write('\t'.join(x for x in line)+'\n')

    print('Gene numbers: '+str(gene_num), 'Isoforms: '+str(mRNA_num), sep='\t')

if __name__ == '__main__':
    main()
