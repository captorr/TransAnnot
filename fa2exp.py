"""
Author: Zhang Chengsheng, @2020.07.30
"""

import os,sys
import re


def faExtract(fa):
    keywords = 'full_length_coverage'
    subKeyword = 'num_subreads'
    result = {}
    with open(fa,'r') as f:
        for i in f.readlines():
            if i.startswith('>'):
                id = i.strip().split(' ')[0].lstrip('>')
                FLC,SRC = 1,0
                for j in i.strip().split(';'):
                    if keywords in j:
                        res = re.findall("{}=(.*)".format(keywords),j)
                        FLC = int(res[0]) if res else FLC
                    if subKeyword in j:
                        res = re.findall("{}=(.*)".format(subKeyword),j)
                        SRC = int(res[0]) if res else SRC
                result[id] = [FLC,SRC]
    return result


def fa2exp(fa,expFile):
    exp = faExtract(fa)
    totalEXP = sum([exp[i][0] for i in exp])
    with open(fa,'r') as f, open(expFile,'w') as o:
        o.write('#ID\tFLC\tTPM\n')
        for i in f.readlines():
            if i.startswith('>'):
                id = i.strip().split(' ')[0].lstrip('>')
                if id in exp:
                    FLC = exp[id][0]
                    TPM = float(FLC)/totalEXP*1000000
                    o.write('{}\t{}\t{}\n'.format(id,FLC,TPM))


def option(argv):
    from argparse import ArgumentParser as AP
    usages = "python3 {} -f fasta -o output".format(argv[0])
    p = AP(usage=usages)
    p.add_argument("-f", dest="fa", metavar="[fa]", help="fasta file", required=True)
    p.add_argument("-o", dest="output", metavar="[outputFile]", help="output file",required=True)

    if len(argv) == 1:
        p.print_help()
        exit(1)
    return p.parse_args(argv[1:])


if __name__ == '__main__':
    args = option(sys.argv)
    fa2exp(args.fa,args.output)
