#! /usr/bin/env python3
# -*- coding:utf-8 -*-


import sys,gzip,argparse
#================================================================================
version = "v1.0"
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''
======================================================================
This is a script for detect UTRs from annotation files that including CDS and exon tag.
Author: Panyouliang, panyouliang@genomics.cn
Version: v1.0
Date: 2023-06-13, yyyy-mm-dd
======================================================================''')

parser.add_argument('-v', '--version', action='version', version=version)
parser.add_argument('-gff', metavar='annotation file', type=str, required=True, help='Please input the annotation file')
args = parser.parse_args()
#================================================================================



def get_gff(gff):

    if gff[-3:] == '.gz':
        f = gzip.open(gff,'rt')
    else:
        f = open(gff,'r')

    cds_dict,exon_dict,all_gff,gene = {},{},{},{}
    gene_num = 0

    for line in f:
        line = line.strip().split('\t')
        
        if line[2] == 'gene':
            gid = line[8].split(';')[0].split('=')[1] +' '+line[6]
            all_gff[gid] = []
            cds_dict[gid],exon_dict[gid] = [],[]
            gene[gid] = [line]

        if line[2] == 'mRNA':
            gene_num += 1
            idx = line[8].find('Parent')
            pgid = line[8][idx:].split(';')[0].split('=')[1]
            tid = line[8].split(';')[0].split('=')[1]
            
            gene[gid].append(line)
            #all_gff[gid].append(line)

        if line[2] == 'CDS' or line[2] == 'exon':
            idx = line[8].find('Parent')
            pmid = line[8][idx:].split(';')[0].split('=')[1]

            if line[2] == 'CDS':
                cds_dict[gid].extend([int(line[4]),int(line[3])])
                all_gff[gid].append(line)

            if line[2] == 'exon':
                exon_dict[gid].extend([int(line[3]),int(line[4])])
                all_gff[gid].append(line)


    for k,v in cds_dict.items():
        v = v.sort()
    for k,v in exon_dict.items():
        v = v.sort()


    outlist = {}

    for k,v in cds_dict.items():

        cds_min,cds_max = min(v),max(v)
        min_list,max_list = [],[]

        min_list = [ x for x in exon_dict[k] if x < cds_min ]
        max_list = [ x for x in exon_dict[k] if x > cds_max ]

        max_list.insert(0,cds_max+1)

        min_list += [cds_min-1]

        outlist[k] = {'max':max_list, 'min':min_list}

    return outlist,all_gff,gene


def detect_utr(dicts):

    UTR5_dict,UTR3_dict = {},{}

    for k,v in dicts.items():

        UTR5_dict[k],UTR3_dict[k] = [],[]

        if k.endswith('-'):

            if len(v['max']) > 1:   # UTR5 detect
                for start,end in zip(v['max'][:-1:2], v['max'][1::2]):
                    UTR5_dict[k].append([start,end])

            if len(v['min']) > 1:   # UTR3 detect
                for start,end in zip(v['min'][:-1:2], v['min'][1::2]):
                    UTR3_dict[k].append([start,end])

        if k.endswith('+'):

            if len(v['max']) > 1:   # UTR3 detect
                for start,end in zip(v['max'][:-1:2], v['max'][1::2]):
                    UTR3_dict[k].append([start,end])

            if len(v['min']) > 1:   # UTR5 detect
                for start,end in zip(v['min'][:-1:2], v['min'][1::2]):
                    UTR5_dict[k].append([start,end])

    UTR5_dict = {k: v for k,v in UTR5_dict.items() if len(v) != 0}
    UTR3_dict = {k: v for k,v in UTR3_dict.items() if len(v) != 0}

    return UTR5_dict,UTR3_dict


def main():

    gff = args.gff

    outlist,all_gff,gene = get_gff(gff)

    UTR5,UTR3 = detect_utr(outlist)

    for k,v in all_gff.items():
        for line in gene[k]:
            print('\t'.join(x for x in line))

        tid = gene[k][1][8].split(';')[0].split('=')[1]

        if k in UTR5.keys():
            if k.endswith('+'):
                for piece in UTR5[k]:
                    print('\t'.join(x for x in gene[k][0][0:2]),'five_prime_utr',piece[0],piece[1],'\t'.join(x for x in gene[k][0][5:8]),'ID='+tid+';Parent='+tid+';',sep='\t')
            if k.endswith('-'):
                for piece in reversed(UTR5[k]):
                    print('\t'.join(x for x in gene[k][0][0:2]),'five_prime_utr',piece[0],piece[1],'\t'.join(x for x in gene[k][0][5:8]),'ID='+tid+';Parent='+tid+';',sep='\t')
        for line in v:
            print('\t'.join(x for x in line))

        if k in UTR3.keys():
            if k.endswith('+'):
                for piece in UTR3[k]:
                    print('\t'.join(x for x in gene[k][0][0:2]),'three_prime_utr',piece[0],piece[1],'\t'.join(x for x in gene[k][0][5:8]),'ID='+tid+';Parent='+tid+';',sep='\t')
            if k.endswith('-'):
                for piece in reversed(UTR3[k]):
                    print('\t'.join(x for x in gene[k][0][0:2]),'three_prime_utr',piece[0],piece[1],'\t'.join(x for x in gene[k][0][5:8]),'ID='+tid+';Parent='+tid+';',sep='\t')





if __name__ == '__main__':
    main()
