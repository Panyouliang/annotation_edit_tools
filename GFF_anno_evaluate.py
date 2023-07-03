#! /usr/bin/env python3
# -*- coding:utf-8 -*-


import sys,gzip,argparse
import numpy as np

#================================================================================
version = "v1.0"
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''
============================================================
This is a script for evaluating the quality of genome annotation.

Author: Panyouliang, panyouliang@genomics.cn
Version: v1.0
Date: 2023-06-13, yyyy-mm-dd
============================================================''')
parser.add_argument('-v', '--version', action='version', version=version)
parser.add_argument('-gff', metavar='annotation file', type=str, required=False, help='Please input the annotation file(None isoform)')
parser.add_argument('-ref',  metavar='genome file', type=str, default="None", required=False, help='Please input the genome file (.fasta)')
parser.add_argument('-sp', metavar='Specie name', type=str, required=True, help='Please input the Specie name')
args = parser.parse_args()
#=================================================================================



def get_gff(gff):

    cds_dict, exon_dict, orf_dict, tran_dict, intron_list, gene_num = {},{},{},{},[],0

    if gff[-3:] == '.gz':
        f = gzip.open(gff,'rt')
    else:
        f = open(gff,'r')

    for line in f:
        if line.strip():
            if line[0] != '#':
                line = line.strip().split('\t')
                if line[2] == 'mRNA':
                    gene_num += 1
                    flag = 0
                    g_id = line[8].split(';')[0].split('=')[1]
                    tran_id = g_id +' '+ line[6]
                    cds_dict[g_id] = []
                    exon_dict[g_id] = []
                    orf_dict[g_id] = []
                    tran_list = []

                else:
                    if line[2] == 'CDS':
                        cds_dict[g_id].append(abs(int(line[4])-int(line[3]))+1)
                        orf_dict[g_id] += [int(line[3]),int(line[4])]
                        exon_dict[g_id].append([line[3],line[4]])
                        tran_pos = [line[0],line[3],line[4]]
                        tran_list.append(tran_pos)
                        tran_dict[tran_id] = tran_list
    gff_dict = {}
    for k,v in tran_dict.items():
        str_list = sorted(v, key = (lambda x:int(x[1])))
        gff_dict[k] = str_list


    return cds_dict, exon_dict, orf_dict, gff_dict, intron_list, gene_num


def get_genome(F1):
    if F1[-3:] == '.gz':
        f = gzip.open(F1,'rt')
    else:
        f = open(F1,'r')
    genome_dict,seqlist,chr = {},[],''

    for line in f:
        line = line.strip()
        if line[0] == '>':
            genome_dict[chr] = ''.join(x for x in seqlist)
            chr = line.split()[0][1:]
            seqlist = []
        else:
            seqlist.append(line)

    genome_dict[chr] = ''.join(x for x in seqlist)
    return genome_dict


def get_mRNA(gff_dict,genome_dict):
    sequence = {}
    for k,v in gff_dict.items():
        if k:
            gid = k
            sequence[gid] = ''
            fragment = ''
            for line in v:
                fragment = genome_dict[line[0]][int(line[1])-1:int(line[2])]
                sequence[gid] += fragment
    return sequence

def rev(seq):
    base = {
            'A':'T','T':'A','G':'C','C':'G','N':'N','n':'n','a':'t','t':'a','c':'g','g':'c'
            }
    seq_list = list(reversed(seq))
    seq_rev = [base[k] for k in seq_list]
    seq_list = ''.join(seq_rev)
    return seq_list

def extract_mRNA(fragment):
    mRNA = {}
    for k,v in fragment.items():
        if k.endswith('-'):
            string = rev(v)
            mRNA[k] = string.upper()
        elif k.endswith('+'):
            mRNA[k] = v.upper()

    return mRNA


code = {
    'GCA':'A','GCC':'A','GCG':'A','GCT':'A',
    'TGC':'C','TGT':'C',
    'GAC':'D','GAT':'D',
    'GAA':'E','GAG':'E',
    'TTC':'F','TTT':'F',
    'GGA':'G','GGC':'G','GGG':'G','GGT':'G',
    'CAC':'H','CAT':'H',
    'ATA':'I','ATC':'I','ATT':'I',
    'AAA':'K','AAG':'K',
    'CTA':'L','CTC':'L','CTG':'L','CTT':'L','TTA':'L','TTG':'L',
    'ATG':'M',
    'AAC':'N','AAT':'N',
    'CCA':'P','CCC':'P','CCG':'P','CCT':'P',
    'CAA':'Q','CAG':'Q',
    'CGA':'R','CGC':'R','CGG':'R','CGT':'R','AGA':'R','AGG':'R',
    'TCA':'S','TCC':'S','TCG':'S','TCT':'S','AGC':'S','AGT':'S',
    'ACA':'T','ACC':'T','ACG':'T','ACT':'T',
    'GTA':'V','GTC':'V','GTG':'V','GTT':'V',
    'TGG':'W',
    'TAC':'Y','TAT':'Y',
    'TAA':'*','TAG':'*','TGA':'*'
    }

stop_codon = ['TAA','TAG','TGA']


def main():

    gff,ref,name = args.gff,args.ref,args.sp

    if gff:

        cds_dict, exon_dict, orf_dict, gff_dict, intron_list, gene_num = get_gff(gff)

        # total genes
        total_gene = gene_num

        # single exon genes
        single_exon_gene = 0
        for k,v in exon_dict.items():
            if len(v) == 1:
                single_exon_gene += 1

        # exons per gene
        cds_list = sum(cds_dict.values(), [])
        exon_per_gene = len(cds_list) / gene_num

        # exon length (mean)
        exon_length_mean = sum(cds_list) / len(cds_list)

        # mRNA length (mean)
        mRNA_length_mean = sum(cds_list) / gene_num

        # ORF length (mean)
        orf = []
        for k,v in orf_dict.items():
            if len(v) > 0:
                mini = min(v)
                maxi = max(v)
                orf.append(abs(maxi-mini)+1)
        gene_length_mean = sum(orf) / len(orf)

        # intron length (mean)
        intron_list = []
        for v in exon_dict.values():
            if len(v) > 1:
                for i in range(len(v)-1):
                    intron_list.append(abs(int(v[i+1][0])-int(v[i][1])))
        intron_length = sum(intron_list) / len(intron_list)


    intact_orf,pseudogene = 0,0

    if ref != 'None':
        fragment = get_mRNA(gff_dict, get_genome(ref))
        all_mRNA = extract_mRNA(fragment)

        for k,v in all_mRNA.items():
            protein = ''
            for i in range(0, len(v), 3):
                codon = v[i:i+3]
                if codon in code:
                    protein += code[codon]
                else:
                    protein += 'X'

            judge_list = list(filter(None,protein.split('*')))
            if len(judge_list) > 1:
                pseudogene += 1

        intactORF_list = []
        for k,v in all_mRNA.items():
            if v[-3:] in stop_codon:
                intactORF_list.append(k)


        for k,v in gff_dict.items():
            if k not in intactORF_list:
                if k.endswith('-'):
                    if int(v[0][1]) > 3:
                        v[0][1] = int(v[0][1]) - 3
                else:
                    v[-1][2] = int(v[-1][2]) + 3


        orf_RNA = extract_mRNA(get_mRNA(gff_dict, get_genome(ref)))

        for k,v in orf_RNA.items():
            if v[0:3] == 'ATG':
                if v[-3:] in stop_codon:
                    intact_orf += 1

    #print('species','total genes','intact-ORF genes', 'intact-ORF proportion', 'pre-terminated genes', 'pre-terminated proportion', 'single exon genes', 'single exon proportion', 'gene length avg (bp)','coding region length avg (bp)','exon length avg (bp)','exons per gene','intron length avg (bp)',sep='\t')

    print(name,total_gene, intact_orf, format(intact_orf/total_gene,'.2f'), pseudogene, format(pseudogene/total_gene,'.4f'), single_exon_gene, format(single_exon_gene/total_gene,'.2f'), int(gene_length_mean), int(mRNA_length_mean), int(exon_length_mean), int(exon_per_gene), int(intron_length), sep='\t')


if __name__ == '__main__':
    main()

