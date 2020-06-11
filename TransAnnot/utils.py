"""

"""
__version__ = '0.0.1'
__email__ = 'zhangchengsheng1@qq.com'
__codex__ = 'https://github.com/captorr/TransAnnot'

import os,pickle
import Long2Short
import Sam2Bed


def load_pickle(file_in):
    return pickle.load(open(file_in, 'rb'))


def pickle_make(dict_in,file_out):
    with open(file_out, 'wb') as o:
        pickle.dump(dict_in, o)


def bed_read(bed):
    isoforms_bed = {}
    with open(bed,'r') as f:
        for idx,line in enumerate(f.readlines()):
            if not line.strip():
                continue
            id = line.split('\t')[3]
            if not idx:
                id_temp = id
                bed_string = line
                continue
            if id == id_temp:
                bed_string += line
            else:
                if id_temp in isoforms_bed:
                    print('error: id repeat in bed files: {}'.format(id_temp))
                isoforms_bed[id_temp] = bed_string
                id_temp = id
                bed_string = line
        else:
            isoforms_bed[id_temp] = bed_string
    return isoforms_bed


def check_overlap(ref, gene, length=0):
    flag = 0
    overlap_length = 0
    if ref[0] <= gene[0]:
        if ref[1] <= gene[0]:
            flag = 0  # 不挨着ref在左
        elif ref[1] >= gene[1]:
            flag = 2  # ref 包含gene
            overlap_length = abs(gene[1] - gene[0] + 1)
        else:
            flag = 1  # overlap
            overlap_length = abs(ref[1] - gene[0] + 1)
    else:
        if ref[0] >= gene[1]:
            flag = 0  # 不挨着ref在右
        elif ref[1] <= gene[1]:
            flag = 2  # gene包含ref
            overlap_length = abs(ref[1] - ref[0] + 1)
        else:
            flag = 1  # overlap
            overlap_length = abs(gene[1] - ref[0] + 1)
    if length:
        return overlap_length
    return flag


def ss_compare(ref,target):
    MAX_DRAFT = 10
    target_start, target_end = [i[0] for i in target], [i[1] for i in target]
    target_range = [min(target_start), max(target_end)]
    if not ref:
        return 0,target,target_range
    ref_start,ref_end = [i[0] for i in ref],[i[1] for i in ref]
    ref_range = [min(ref_start),max(ref_end)]
    tk,rk=0,0
    match = []
    for t_idx in range(len(target_start)):
        start_match,end_match = 0,0
        if not check_overlap(ref_range,target[t_idx]):
            continue
        for j in range(len(ref_start[rk:])):
            r_idx = rk+j
            if abs(target_start[t_idx]-ref_start[r_idx]) < MAX_DRAFT:
                start_match = 1
            if abs(target_end[t_idx]-ref_end[r_idx]) < MAX_DRAFT:
                end_match = 1
            if start_match or end_match:
                rk += j
                break
        if not start_match+end_match:
            match = []
            break
        match.append(start_match+end_match)
    same = 0
    if len(match) > 1:
        if 1 not in match[1:-1]:
            same = 1
            if abs(ref_range[-1]-ref_range[0]) > abs(target_range[-1]-target_range[0]):
                target = ref
                target_range = ref_range
    return same,target,target_range


def runSplit(c):
    faSplit = os.path.join(c['OUTPUT_DIR'],c['SAMPLE_UNIQUE_NAME']+c['FA_SPLIT_SUFFIX'])
    Long2Short.main(c['FASTA'],faSplit,int(c['READ_LENGTH']),int(c['READ_OVERLAP']),int(c['MIN_READ_LENGTH']))


def runSam2Bed(sam,bed,Type,process,readLength,readOverlap):
    Sam2Bed.main(sam,bed,type=Type,process=int(process),READ_LENGTH=int(readLength),OVERLAP_LENGTH=int(readOverlap))


def runMinimap2(c):
    print('run minimap2')
    sam = os.path.join(c['OUTPUT_DIR'],c['SAMPLE_UNIQUE_NAME']+c['SAM_SUFFIX_MINIMAP2'])
    samSort = os.path.join(c['OUTPUT_DIR'],c['SAMPLE_UNIQUE_NAME']+c['SAM_SORT_SUFFIX_MINIMAP2'])
    cmd1 = '{} -ax splice -t {} -uf --secondary=no -C5 {} {} > {}'.format(c['MINIMAP2'],c['PROCESS'],c['GENOME_FA'],c['FASTA'],sam)
    cmd2 = '{} sort -n -o {} {}'.format(c['SAMTOOLS'],samSort,sam)
    os.system(cmd1)
    os.system(cmd2)
    os.remove(sam)
    bed = os.path.join(c['OUTPUT_DIR'],c['SAMPLE_UNIQUE_NAME']+c['BED_SUFFIX_MINIMAP2'])
    runSam2Bed(samSort,bed,'long',c['PROCESS'],c['READ_LENGTH'],c['READ_OVERLAP'])
    return bed


def runHisat2(c):
    print('run hisat2')
    faSplit = os.path.join(c['OUTPUT_DIR'],c['SAMPLE_UNIQUE_NAME']+c['FA_SPLIT_SUFFIX'])
    sam = os.path.join(c['OUTPUT_DIR'],c['SAMPLE_UNIQUE_NAME']+c['SAM_SUFFIX_HISAT2'])
    samSort = os.path.join(c['OUTPUT_DIR'],c['SAMPLE_UNIQUE_NAME']+c['SAM_SORT_SUFFIX_HISAT2'])
    cmd1 = '{} -p {} -f --score-min L,0,-0.8 -x {} -U {} -S {}'.format(c['HISAT2'],c['PROCESS'],c['HISAT2_INDEX'],faSplit,sam)
    cmd2 = '{} sort -n -o {} {}'.format(c['SAMTOOLS'],samSort,sam)
    os.system(cmd1)
    os.system(cmd2)
    os.remove(sam)
    bed = os.path.join(c['OUTPUT_DIR'], c['SAMPLE_UNIQUE_NAME'] + c['BED_SUFFIX_HISAT2'])
    runSam2Bed(samSort, bed, 'short', c['PROCESS'], c['READ_LENGTH'], c['READ_OVERLAP'])
    return bed


def runGmap(c):
    print('run gmap')
    sam = os.path.join(c['OUTPUT_DIR'], c['SAMPLE_UNIQUE_NAME'] + c['SAM_SUFFIX_GMAP'])
    samSort = os.path.join(c['OUTPUT_DIR'], c['SAMPLE_UNIQUE_NAME'] + c['SAM_SORT_SUFFIX_GMAP'])
    cmd1 = ''
    cmd2 = '{} sort -n -o {} {}'.format(c['SAMTOOLS'], samSort, sam)
    #os.system(cmd1)
    #os.system(cmd2)
    #os.remove(sam)
    bed = os.path.join(c['OUTPUT_DIR'], c['SAMPLE_UNIQUE_NAME'] + c['BED_SUFFIX_GMAP'])
    #runSam2Bed(samSort, bed, 'long', c['PROCESS'], c['READ_LENGTH'], c['READ_OVERLAP'])
    return bed

