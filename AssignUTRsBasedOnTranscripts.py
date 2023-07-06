#! /usr/bin/env python3
# -*- coding:utf-8 -*-

import sys,gzip,argparse
import subprocess
import numpy as np

#================================================================================
version = "v1.0"
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''
================================================================
This is a script for evaluating the quality of genome annotation.

Author: Panyouliang, panyouliang@genomics.cn
Version: v1.0
Date: 2023-07-04, yyyy-mm-dd
================================================================''')
parser.add_argument('-v', '--version', action='version', version=version)
parser.add_argument('-core_gff', metavar='annotation file', type=str, required=True, help='Please input the annotation file(Select Top length isoforms before input)')
parser.add_argument('-Iso_gff',  metavar='genome file', type=str, required=True, help='Please input the stringtie assemble file of ISO-seq(gff format)')
parser.add_argument('-ngs_gff',  metavar='genome file', type=str, required=False, help='Please input the stringtie assemble file of RNA-seq(gff format)')
parser.add_argument('-sp', metavar='Specie name', type=str, required=False, default="husky", help='Please input the Specie name')
args = parser.parse_args()
#=================================================================================


def read_gff_CDS(F):
    if F[-3:] == '.gz':
        f = gzip.open(F,'rt')
    else:
        f = open(F,'r')

    gene_dict,gene2Iso,all_gene = {},{},{}
    gene_tag = {}
    gemoma_utr5,gemoma_utr5_lis,gemoma_utr3,gemoma_utr3_lis = {},{},{},{}
    
    UTR3 = ['three_prime_UTR','three_prime_utr','UTR3']
    UTR5 = ['five_prime_UTR','five_prime_utr','UTR5']

    for line in f.readlines():
        if line.strip():
            if line[0] != '#':
                line = line.strip().split('\t')
                if line[2] == 'gene':
                    geneID = line[8].split('=')[1].split(';')[0]+' '+line[6]
                    gene_tag[geneID] = line
                if line[2] == 'mRNA':
                    gidx,tidx = line[8].find('Parent'),line[8].find('ID')
                    gid,tid = line[8][gidx:].split('=')[1].split(';')[0] + ' ' + line[6],line[8][tidx:].split('=')[1].split(';')[0]
                    gene_dict[gid] = []
                    gene2Iso[gid] = tid
                    all_gene[gid] = {tid:[line]}

                if line[2] == 'CDS':
                    tidx = line[8].find('Parent')
                    tid = line[8][tidx:].split('=')[1].split(';')[0]
                    gene_dict[gid].extend([int(line[3]),int(line[4])])
                    all_gene[gid][tid].append(line)

                if line[2] in UTR5:
                    if geneID in gemoma_utr5.keys():
                        gemoma_utr5[geneID].append([int(line[3]),int(line[4])])
                        gemoma_utr5_lis[geneID].extend([int(line[3]),int(line[4])])
                    else:
                        gemoma_utr5[geneID] = [[int(line[3]),int(line[4])]]
                        gemoma_utr5_lis[geneID] = [int(line[3]),int(line[4])]

                if line[2] in UTR3:
                    if geneID in gemoma_utr3.keys():
                        gemoma_utr3[geneID].append([int(line[3]),int(line[4])])
                        gemoma_utr3_lis[geneID].extend([int(line[3]),int(line[4])])
                    else:
                        gemoma_utr3[geneID] = [[int(line[3]),int(line[4])]]
                        gemoma_utr3_lis[geneID] = [int(line[3]),int(line[4])]

    for gid,interval in gene_dict.items():
        interval = interval.sort()

    return gene_dict,gene2Iso,all_gene,gene_tag,gemoma_utr5,gemoma_utr3,gemoma_utr5_lis,gemoma_utr3_lis

def read_gtf_exon(F):
    if F[-3:] == '.gz':
        f = gzip.open(F,'rt')
    else:
        f = open(F,'r')

    trans_dict,transcript2gene = {},{}

    for line in f.readlines():
        if line.strip():
            if line[0] != '#':
                line = line.strip().split('\t')
                if line[2] == 'transcript':
                    gidx,tidx = line[8].find('geneID'),line[8].find('ID')
                    gid,tid = line[8][gidx:].split('=')[1].split(';')[0],line[8][tidx:].split('=')[1].split(';')[0]
                    #transcript2gene[tid] = gid
                    trans_dict[tid] = []

                if line[2] == 'exon':
                    ptidx = line[8].find('Parent')
                    ptid = line[8][ptidx:].split('=')[1].split(';')[0]

                    if ptid in trans_dict.keys():
                        trans_dict[ptid].extend([int(line[3]),int(line[4])])

    for tid,interval in trans_dict.items():
        interval = interval.sort()

    return trans_dict #,transcript2gene


def stat_overlap(gff1,gff2):

    cdsgID2exongID = {}

    result = subprocess.run(['perl', '/ldfssz1/ST_EARTH/P18Z10200N0107/panyouliang/my_work/00.Tools/02.Genome_Anno_evaluate/annotation_script/FindOverlapAtCDSlevel.exon.pl', gff1, gff2], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    output = result.stdout.decode()

    for line in output.splitlines():
        if line[0] != '#':
            line = line.strip().split('\t')
            if float(line[10]) >= 0.7:
                if line[0] in cdsgID2exongID.keys():
                    cdsgID2exongID[line[0]].append(line[1])
                else:
                    cdsgID2exongID[line[0]] = [line[1]]

    return cdsgID2exongID


def get_over(F):
    cdsgID2exongID = {}
    with open(F,'r') as f:
        for line in f:
            if line[0] != '#':
                line = line.strip().split('\t')
                if float(line[10]) >= 0.7 and float(line[11]) >= 0.3:
                    if line[0] in cdsgID2exongID.keys():
                        cdsgID2exongID[line[0]].append(line[1])
                    else:
                        cdsgID2exongID[line[0]] = [line[1]]

    return cdsgID2exongID




def cut_UTR(gene_dict,cdsgID2exongID,gene2Iso,trans_dict):

    UTR3,UTR5 = {},{}
    min_dict,max_dict = {},{}
    for gid,cds_interval in gene_dict.items():

        cds_min,cds_max = min(cds_interval),max(cds_interval)
        min_dict[gid],max_dict[gid] = {},{}
        max_len,min_len = {},{}
        if gene2Iso[gid] in cdsgID2exongID.keys():
            for exon_ptid in cdsgID2exongID[gene2Iso[gid]]:
                min_dict[gid][exon_ptid] = [x for x in trans_dict[exon_ptid] if x < cds_min ]
                max_dict[gid][exon_ptid] = [x for x in trans_dict[exon_ptid] if x > cds_max ]

                if len(max_dict[gid][exon_ptid]) > 0:
                    max_len[exon_ptid] = 0
                    if len(max_dict[gid][exon_ptid]) % 2 == 1:
                        max_dict[gid][exon_ptid].insert(0,cds_max+1)
                        for start,end in zip(max_dict[gid][exon_ptid][:-1:2], max_dict[gid][exon_ptid][1::2]):
                            max_len[exon_ptid] += abs(end-start)
                    else:
                        for start,end in zip(max_dict[gid][exon_ptid][:-1:2], max_dict[gid][exon_ptid][1::2]):
                            max_len[exon_ptid] += abs(end-start)

                if len(min_dict[gid][exon_ptid]) > 0:
                    min_len[exon_ptid] = 0
                    if len(min_dict[gid][exon_ptid]) % 2 == 1:
                        min_dict[gid][exon_ptid] += [cds_min-1]
                        for start,end in zip(min_dict[gid][exon_ptid][:-1:2], min_dict[gid][exon_ptid][1::2]):
                            min_len[exon_ptid] += abs(end-start)
                    else:
                        for start,end in zip(min_dict[gid][exon_ptid][:-1:2], min_dict[gid][exon_ptid][1::2]):
                            min_len[exon_ptid] += abs(end-start)

            if len(max_len) > 0:
                maxk = max(max_len, key=max_len.get)
                if gid.endswith('+'):
                    UTR3[gid] = max_dict[gid][maxk]
                if gid.endswith('-'):
                    UTR5[gid] = max_dict[gid][maxk]

            if len(min_len) > 0:
                mink = max(min_len, key=min_len.get)
                if gid.endswith('+'):
                    UTR5[gid] = min_dict[gid][mink]
                if gid.endswith('-'):
                    UTR3[gid] = min_dict[gid][mink]


    UTR3_dict = {}
    for gid,exon_list in UTR3.items():
        UTR3_dict[gid] = []
        for start,end in zip(exon_list[:-1:2], exon_list[1::2]):
            UTR3_dict[gid].append([start,end])

    UTR5_dict = {}
    for gid,exon_list in UTR5.items():
        UTR5_dict[gid] = []
        for start,end in zip(exon_list[:-1:2], exon_list[1::2]):
            UTR5_dict[gid].append([start,end])


    return UTR3_dict,UTR5_dict,UTR5,UTR3

def update_gene_coordinate(gene_dict,UTR5,UTR3):
    gene_pos = {}
    for gid,interval in gene_dict.items():
        if gid in UTR3.keys() and gid in UTR5.keys():
            all_pos = interval + UTR3[gid] + UTR5[gid]
        if gid in UTR3.keys() and gid not in UTR5.keys():
            all_pos = interval + UTR3[gid]
        if gid in UTR5.keys() and gid not in UTR3.keys():
            all_pos = interval + UTR5[gid]
        if gid not in UTR5.keys() and gid not in UTR3.keys():
            all_pos = interval

        gene_pos[gid] = [min(all_pos),max(all_pos)]

    return gene_pos





def out_print(all_gene,gene_pos,UTR5_dict,UTR3_dict,gene_tag):
    for gid,m_value in all_gene.items():
        utr5,utr3,exon = 0,0,0
        for k,v in m_value.items():
            tid = v[0][8].split('=')[1].split(';')[0]
            gidx = v[0][8].find('Parent')
            geneID = v[0][8][gidx:].split('=')[1].split(';')[0]
            print('\t'.join(x for x in v[0][0:2]),'gene',gene_pos[gid][0],gene_pos[gid][1],'\t'.join(x for x in v[0][5:8]),'ID='+geneID+';',sep='\t')
            print('\t'.join(x for x in v[0][0:3]),gene_pos[gid][0],gene_pos[gid][1],'\t'.join(x for x in v[0][5:]),sep='\t')


        if gid in UTR5_dict.keys():
            for piece in UTR5_dict[gid]:
                utr5 += 1
                exon += 1
                print('\t'.join(x for x in gene_tag[gid][0:2]),'five_prime_UTR',piece[0],piece[1],'\t'.join(x for x in gene_tag[gid][5:8]),'Parent='+tid+';',sep='\t')
                print('\t'.join(x for x in gene_tag[gid][0:2]),'exon',piece[0],piece[1],'\t'.join(x for x in gene_tag[gid][5:8]),'Parent='+tid+';',sep='\t')


        for tid,cds in m_value.items():
            for line in cds[1:]:
                exon += 1
                print('\t'.join(x for x in line))
                print('\t'.join(x for x in gene_tag[gid][0:2]),'exon',line[3],line[4],'\t'.join(x for x in gene_tag[gid][5:8]),'Parent='+tid+';',sep='\t')

        if gid in UTR3_dict.keys():
            for piece in UTR3_dict[gid]:
                utr3 += 1
                exon += 1
                print('\t'.join(x for x in gene_tag[gid][0:2]),'exon',piece[0],piece[1],'\t'.join(x for x in gene_tag[gid][5:8]),'Parent='+tid+';',sep='\t')
                print('\t'.join(x for x in gene_tag[gid][0:2]),'three_prime_UTR',piece[0],piece[1],'\t'.join(x for x in gene_tag[gid][5:8]),'Parent='+tid+';',sep='\t')




def main():

    gff1,gff2,gff3,sp = args.core_gff,args.Iso_gff,args.ngs_gff,args.sp

    gene_dict,gene2trans,all_gene,gene_tag,gemoma_utr5,gemoma_utr3,gemoma_utr5_lis,gemoma_utr3_lis = read_gff_CDS(gff1)
    
    trans_dict_Iso = read_gtf_exon(gff2)

    cdsgID2exongID_Iso = stat_overlap(gff1,gff2)
    #cdsgID2exongID_Iso = get_over('overlap.txt')

    UTR3_iso_dict,UTR5_iso_dict,UTR5_iso,UTR3_iso = cut_UTR(gene_dict,cdsgID2exongID_Iso,gene2trans,trans_dict_Iso)
    gene_pos_iso = update_gene_coordinate(gene_dict,UTR5_iso,UTR3_iso)

    if gff3:
        cdsgID2exongID_ngs = stat_overlap(gff1,gff3)
        trans_dict_ngs = read_gtf_exon(gff3)
        #cdsgID2exongID_ngs = get_over('overlap.ngs.txt')
        UTR3_ngs_dict,UTR5_ngs_dict,UTR5_ngs,UTR3_ngs = cut_UTR(gene_dict,cdsgID2exongID_ngs,gene2trans,trans_dict_ngs)
        gene_pos_ngs = update_gene_coordinate(gene_dict,UTR5_ngs,UTR3_ngs)

        gene_pos = {}
        for gid,cds_v in all_gene.items():
            utr5,utr3,exon = 0,0,0
            if gid in UTR5_iso_dict.keys() or gid in UTR5_ngs_dict.keys():
                if gid in UTR5_iso_dict.keys():
                    pos_utr5 = UTR5_iso[gid]
                else:
                    pos_utr5 = UTR5_ngs[gid]
            else:
                if gid in gemoma_utr5_lis.keys():
                    pos_utr5 = gemoma_utr5_lis[gid]
                else:
                    pos_utr5 = gene_dict[gid]

            if gid in UTR3_iso_dict.keys() or gid in UTR3_ngs_dict.keys():
                if gid in UTR3_iso_dict.keys():
                    pos_utr3 = UTR3_iso[gid]
                else:
                    pos_utr3 = UTR3_ngs[gid]
            else:
                if gid in gemoma_utr3_lis.keys():
                    pos_utr3 = gemoma_utr3_lis[gid]
                else:
                    pos_utr3 = gene_dict[gid]


            interval = pos_utr5 + pos_utr3
            gene_pos[gid] = [min(interval),max(interval)]


            for tid,coding in cds_v.items():
                ggidx = coding[0][8].find('Parent')
                geneID = coding[0][8][ggidx:].split('=')[1].split(';')[0]
                print('\t'.join(x for x in coding[0][0:2]),'gene',gene_pos[gid][0],gene_pos[gid][1],'\t'.join(x for x in coding[0][5:8]),'ID='+geneID+';',sep='\t')
                print('\t'.join(x for x in coding[0][0:3]),gene_pos[gid][0],gene_pos[gid][1],'\t'.join(x for x in coding[0][5:]),sep='\t')

            if gid in UTR5_iso_dict.keys() or gid in UTR5_ngs_dict.keys():
                if gid in UTR5_iso_dict.keys():
                    for piece in UTR5_iso_dict[gid]:
                        print('\t'.join(x for x in gene_tag[gid][0:2]),'UTR5',piece[0],piece[1],'\t'.join(x for x in gene_tag[gid][5:8]),'Parent='+tid+';',sep='\t')
                        print('\t'.join(x for x in gene_tag[gid][0:2]),'exon',piece[0],piece[1],'\t'.join(x for x in gene_tag[gid][5:8]),'Parent='+tid+';',sep='\t')
                else:
                    for piece in UTR5_ngs_dict[gid]:
                        print('\t'.join(x for x in gene_tag[gid][0:2]),'UTR5',piece[0],piece[1],'\t'.join(x for x in gene_tag[gid][5:8]),'Parent='+tid+';',sep='\t')
                        print('\t'.join(x for x in gene_tag[gid][0:2]),'exon',piece[0],piece[1],'\t'.join(x for x in gene_tag[gid][5:8]),'Parent='+tid+';',sep='\t')

            if gid not in UTR5_iso_dict.keys() and gid not in UTR5_ngs_dict.keys():
                if gid in gemoma_utr5.keys():
                    for piece in gemoma_utr5[gid]:
                        utr5 += 1
                        print('\t'.join(x for x in gene_tag[gid][0:2]),'UTR5',piece[0],piece[1],'\t'.join(x for x in gene_tag[gid][5:8]),'Parent='+tid+';',sep='\t')
                        print('\t'.join(x for x in gene_tag[gid][0:2]),'exon',piece[0],piece[1],'\t'.join(x for x in gene_tag[gid][5:8]),'Parent='+tid+';',sep='\t')

            for tid,cds in cds_v.items():
                for line in cds[1:]:
                    exon += 1
                    print('\t'.join(x for x in line))
                    print('\t'.join(x for x in gene_tag[gid][0:2]),'exon',line[3],line[4],'\t'.join(x for x in gene_tag[gid][5:8]),'Parent='+tid+';',sep='\t')


            if gid in UTR3_iso_dict.keys() or gid in UTR3_ngs_dict.keys():
                if gid in UTR3_iso_dict.keys():
                    for piece in UTR3_iso_dict[gid]:
                        print('\t'.join(x for x in gene_tag[gid][0:2]),'UTR3',piece[0],piece[1],'\t'.join(x for x in gene_tag[gid][5:8]),'Parent='+tid+';',sep='\t')
                        print('\t'.join(x for x in gene_tag[gid][0:2]),'exon',piece[0],piece[1],'\t'.join(x for x in gene_tag[gid][5:8]),'Parent='+tid+';',sep='\t')

                else:
                    for piece in UTR3_ngs_dict[gid]:
                        print('\t'.join(x for x in gene_tag[gid][0:2]),'UTR3',piece[0],piece[1],'\t'.join(x for x in gene_tag[gid][5:8]),'Parent='+tid+';',sep='\t')
                        print('\t'.join(x for x in gene_tag[gid][0:2]),'exon',piece[0],piece[1],'\t'.join(x for x in gene_tag[gid][5:8]),'Parent='+tid+';',sep='\t')

            if gid not in UTR3_iso_dict.keys() and gid not in UTR3_ngs_dict.keys():
                if gid in gemoma_utr3.keys():
                    for piece in gemoma_utr3[gid]:
                        utr3 += 1
                        print('\t'.join(x for x in gene_tag[gid][0:2]),'exon',piece[0],piece[1],'\t'.join(x for x in gene_tag[gid][5:8]),'Parent='+tid+';',sep='\t')
                        print('\t'.join(x for x in gene_tag[gid][0:2]),'UTR3',piece[0],piece[1],'\t'.join(x for x in gene_tag[gid][5:8]),'Parent='+tid+';',sep='\t')
                        

    else:
        out_print(all_gene,gene_pos_iso,UTR5_iso_dict,UTR3_iso_dict,gene_tag)



if __name__ == '__main__':
    main()

