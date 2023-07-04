#! /usr/bin/env python3
# -*- coding:utf-8 -*-

import sys,gzip,argparse
#================================================================================
version = "v1.0"
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''
======================================================================
This is a script for transform gtf to gff file.
Author: Panyouliang, panyouliang@genomics.cn
Version: v1.0
Date: 2023-07-03, yyyy-mm-dd
======================================================================''')
parser.add_argument('-v', '--version', action='version', version=version)
parser.add_argument('-filelist', metavar='The path list', type=str, required=True, help='Please input the path file')
#parser.add_argument('-sp', metavar='Specie name', type=str, required=True, help='Please input the Specie name')
args = parser.parse_args()
#================================================================================


def read_gff(F,UTR5,UTR3):
    all_cds = {}
    trans = {}
    with open(F,'r') as f:
        for line in f:
            line = line.strip().split()
            if line[2] == 'mRNA':
                ids = line[8].split(';')[0].split('=')[1]+' '+line[6]
                trans[ids] = line
                all_cds[ids] = []
                if ids not in UTR5.keys():
                    UTR5[ids] = []
                if ids not in UTR3.keys():
                    UTR3[ids] = []

            if line[2] == 'UTR5':
                pid = line[8].split(';')[0].split('=')[1]+' '+line[6]
                pos = [line[3],line[4]]
                UTR5[ids].append(pos)

            if line[2] == 'UTR3':
                pid = line[8].split(';')[0].split('=')[1]+' '+line[6]
                pos = [line[3],line[4]]
                UTR3[ids].append(pos)

            if line[2] == 'CDS':
                pid = line[8].split(';')[0].split('=')[1]+' '+line[6]
                all_cds[pid].append(line)

    return UTR5,UTR3,all_cds,trans


def merge_utr_locus(dicts):
    new_utr_list = {}
    for tid,utr_list in dicts.items():
        utr_list = sorted(utr_list, key=lambda x:x[0], reverse=True)
        merge_position = []
        for interval in utr_list:
            if not merge_position:
                merge_position.append(interval)
            else:
                merged = False
                for utr_range in merge_position:
                    if (interval[0] >= utr_range[0] and interval[0] <= utr_range[1]) or (interval[1] >= utr_range[0] and interval[1] <= utr_range[1]):
                        utr_range[0] = min(utr_range[0],interval[0])
                        utr_range[1] = max(utr_range[1],interval[1])
                        merged = True
                        break
                if not merged:
                    merge_position.append(interval)

        new_utr_list[tid] = merge_position

    return new_utr_list


def main():
    f = args.filelist
    UTR5,UTR3 = {},{}
    for input_f in open(f,'r'):
        F = input_f.strip()
        UTR5,UTR3,all_cds,trans = read_gff(F,UTR5,UTR3)

    UTR5,UTR3 = merge_utr_locus(UTR5),merge_utr_locus(UTR3)

    for tid,cds_list in all_cds.items():

        print('\t'.join(x for x in trans[tid]))

        for interval in UTR5[tid]:
            types5 = trans[tid][:]
            types5[2],types5[3],types5[4],types5[8] = 'UTR5',interval[0],interval[1],'Parent='+tid.split()[0]
            print('\t'.join(x for x in types5))

        for line in cds_list:
            print('\t'.join(x for x in line))

        for interval in UTR3[tid]:
            types3 = trans[tid][:]
            types3[2],types3[3],types3[4],types3[8] = 'UTR3',interval[0],interval[1],'Parent='+tid.split()[0]
            print('\t'.join(x for x in types3))


if __name__ == '__main__':
    main()
