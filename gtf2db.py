"""
Author: Zhang Chengsheng, @2020.03.06
"""

import re,sys
import pickle


def dict_make(gtf,source='ensembl',gtf_out=''):
    global db_dict
    db_dict = {}
    with open(gtf,'r') as f:
        for line in f.readlines():
            if line.startswith('#'):
                continue
            context_parse(line,source)
    db_patch()
    if not gtf_out:
        return db_dict
    else:
        pickle_make(db_dict,gtf_out)


def pickle_make(dict_in,file_out):
    with open(file_out, 'wb') as o:
        pickle.dump(dict_in, o)


def load_pickle(file_in):
    return pickle.load(open(file_in, 'rb'))


def context_parse(line, source='ensembl'):
    global db_dict
    header = ['chr', 'db', 'record', 'start', 'end', 'tmp', 'strand', 'tmp', 'info']
    keywords = ['gene_id','transcript_id','gene_name','transcript_name','exon_number','gene_biotype']
    geneKeyword = 'gene_name'
    transcriptKeyword = 'transcript_id' if source in ['ncbi'] else 'transcript_name'
    context = line.strip('\n').split('\t')
    record = context[2]
    gene_biotype = 'unknown'
    chr = context[0]
    start, end = int(context[3]), int(context[4])
    start, end = min(start, end), max(start, end)
    strand = context[6]  # +/-
    #strings = context[8]
    if chr not in db_dict:
        db_dict[chr] = {}
    if record == 'gene':
        gene_id = re_find_keyword(line,geneKeyword)
        if gene_id not in db_dict[chr]:
            db_dict[chr][gene_id] = [gene_biotype,strand,[start,end],{}]
        else:
            if not db_dict[chr][gene_id][2]:
                db_dict[chr][gene_id][2] = [start,end]
            if not db_dict[chr][gene_id][1]:
                db_dict[chr][gene_id][1] = strand
    if record == 'transcript':
        gene_id = re_find_keyword(line, geneKeyword)
        transcript_id = re_find_keyword(line, transcriptKeyword)
        if gene_id not in db_dict[chr]:
            db_dict[chr][gene_id] = [gene_biotype,strand,[],{}]
        if transcript_id not in db_dict[chr][gene_id][3]:
            db_dict[chr][gene_id][3][transcript_id] = [strand,[start,end],{}]
        else:
            if not db_dict[chr][gene_id][3][transcript_id][1]:
                db_dict[chr][gene_id][3][transcript_id][1] = [start,end]
    if record == 'exon':
        gene_id = re_find_keyword(line, geneKeyword)
        transcript_id = re_find_keyword(line, transcriptKeyword)
        if gene_id not in db_dict[chr]:
            db_dict[chr][gene_id] = [gene_biotype,strand,[],{}]
        if transcript_id not in db_dict[chr][gene_id][3]:
            db_dict[chr][gene_id][3][transcript_id] = [strand,[],{}]
        if not list(db_dict[chr][gene_id][3][transcript_id][2]):
            exon_number = 1
        else:
            exon_number = max(list(db_dict[chr][gene_id][3][transcript_id][2]))+1
        if exon_number not in db_dict[chr][gene_id][3][transcript_id][2]:
            db_dict[chr][gene_id][3][transcript_id][2][exon_number] = [start,end]


def re_find_keyword_wtf_gencode(strings,keyword='exon_number'):
    return re.findall(keyword + ' (.*?);', strings)[0]


def re_find_keyword(strings,keyword='gene_name'):
    return re.findall(keyword + ' "(.*?)"', strings)[0]


def db_patch():
    def _get_transcript_region(dict):
        p = []
        for exon in dict:
            p.append(dict[exon][0])
            p.append(dict[exon][1])
        return [min(p),max(p)]

    def _get_gene_region(dict):
        p = []
        for transcript in dict:
            p.append(dict[transcript][1][0])
            p.append(dict[transcript][1][1])
        return [min(p),max(p)]

    global db_dict
    for chr in db_dict:
        for gene in db_dict[chr]:
            gene_exon_idx = {}  # 4
            gene_exon = []
            gene_ss = []  # 5
            gene_exon_num = []  # 6
            gene_transcipt_exon_idx = []  # 7
            for transcript in db_dict[chr][gene][3]:
                if not db_dict[chr][gene][3][transcript][1]:
                    db_dict[chr][gene][3][transcript][1] = _get_transcript_region(db_dict[chr][gene][3][transcript][2])
                for exon in db_dict[chr][gene][3][transcript][2]:
                    e1,e2 = db_dict[chr][gene][3][transcript][2][exon]
                    if db_dict[chr][gene][3][transcript][2][exon] not in gene_exon:
                        gene_exon.append(db_dict[chr][gene][3][transcript][2][exon])
                    if e1 not in gene_ss:
                        gene_ss.append(e1)
                    if e2 not in gene_ss:
                        gene_ss.append(e2)
                exon_num = len(db_dict[chr][gene][3][transcript][2])
                gene_exon_num.append(exon_num)
            gene_exon = sorted(gene_exon,key=lambda x:(x[0],x[1]))
            for _idx,i in enumerate(gene_exon):
                gene_exon_idx[_idx] = i
                for transcript in db_dict[chr][gene][3]:
                    transcript_exon_idx = []
                    for exon in db_dict[chr][gene][3][transcript][2]:
                        transcript_exon_idx.append(gene_exon.index(db_dict[chr][gene][3][transcript][2][exon]))
                    db_dict[chr][gene][3][transcript].append(transcript_exon_idx)
                    gene_transcipt_exon_idx.append(transcript_exon_idx)
            if not db_dict[chr][gene][2]:
                db_dict[chr][gene][2] = _get_gene_region(db_dict[chr][gene][3])
            db_dict[chr][gene].append(gene_exon_idx)
            db_dict[chr][gene].append(gene_ss)
            db_dict[chr][gene].append(gene_exon_num)
            db_dict[chr][gene].append(gene_transcipt_exon_idx)


def optinos(argv):
    from argparse import  ArgumentParser as AP
    usages = "python3 {} -i gtf -o db".format(argv[0])
    p = AP(usage=usages)
    p.add_argument("-i", dest="gtf", metavar="gtf_file", help="Annotation file in format of GTF",required=True)
    p.add_argument("-o", dest="db", metavar="file_out", help="output file path",required=True)
    p.add_argument("-s", dest="source", metavar="source", help="source of GTF file,[ensembl,ncbi,gencode], default:ensembl", type=str, default='ensembl',choices=['ensembl','ncbi','gencode'])
    if len(argv) == 1:
        p.print_help()
        exit(1)
    return p.parse_args(argv[1:])


if __name__ == '__main__':
    args = optinos(sys.argv)
    dict_make(args.gtf,source=args.source,gtf_out=args.db)
