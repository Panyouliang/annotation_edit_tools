#! /usr/bin/env python3
# -*- coding:utf-8 -*-


import sys


def read_gff(F):
    all_gene,gene_d,all_mRNA,coding_rangeLen,cds_num_d = {},{},{},{},{}
    with open(F,'r') as f:
        for line in f:
            if line.strip():
                if line[0] != '#':
                    line = line.strip().split('\t')
                    if line[2] == 'gene':
                        geneID = line[8].split('=')[1].split(';')[0]
                        all_gene[geneID] = {}
                        all_mRNA[geneID] = {}
                        coding_rangeLen[geneID] = []
                        cds_num_d[geneID] = []
                        gene_d[geneID] = line
                    else:
                        if line[2] == 'mRNA':
                            pidx = line[8].find('Parent')
                            pid = line[8][pidx:].split('=')[1].split(';')[0]
                            tid = line[8].split('=')[1].split(';')[0]
                            all_gene[pid][tid] = [line]
                            all_mRNA[pid][tid] = [0,0]
                            cds_num,cds_len = 0,0
                        else:
                            tpidx = line[8].find('Parent')
                            tpid = line[8][tpidx:].split('=')[1].split(';')[0]
                            all_gene[pid][tpid].append(line)

                            if line[2] == 'CDS':
                                cds_num += 1
                                cpidx = line[8].find('Parent')
                                cpid = line[8][cpidx:].split('=')[1].split(';')[0]
                                cds_len += abs(int(line[4]) - int(line[3])) + 1
                                coding_rangeLen[pid].append(cds_len)
                                cds_num_d[pid].append(cds_num)
                                all_mRNA[pid][cpid] = [cds_num,cds_len]

    return all_gene, gene_d, all_mRNA,coding_rangeLen,cds_num_d


def main():
    gff = sys.argv[1]

    all_gene,gene_d,all_mRNA,coding_rangeLen,cds_num_d = read_gff(gff)

    for geneID,mRNAInfo in all_mRNA.items():

        if max(cds_num_d[geneID]) == 1:
            if max(coding_rangeLen[geneID]) >= 300:
                print('\t'.join(x for x in gene_d[geneID]))
        else:
            if max(coding_rangeLen[geneID]) >= 150:
                print('\t'.join(x for x in gene_d[geneID]))

        for mrnaID,cdsLen in mRNAInfo.items():
            if cdsLen[0] == 1 and cdsLen[1] >= 300:
                for line in all_gene[geneID][mrnaID]:
                    print('\t'.join(x for x in line))

            elif cdsLen[0] >= 2 and cdsLen[1] >= 150:
                for line in all_gene[geneID][mrnaID]:
                    print('\t'.join(x for x in line))

if __name__ == '__main__':
    main()







