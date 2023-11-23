"""
Author: Zhang Chengsheng, @2020.04.22
"""

import os,sys,random,pickle


def refseq_extract(site_chr,site_site,range,fa,strand=1,flank=0):
    faidx = fa+'.fai'
    if not os.path.exists(faidx):
        print('ERROR: {} not found'.format(faidx))
        exit(1)
    prefix = ''  # chr
    site_chr = prefix + site_chr.lstrip('chr')
    faidx_dict = {}
    with open(faidx,'r') as o:
        for i in o.readlines():
            faidx_dict[i.strip().split('\t')[0]] = i.strip().split('\t')[1:]
    fbuffer = open(fa,'r')
    if flank:
        start,end = int(site_site)-flank,int(site_site)+flank
    else:
        start,end = range
    offset = int(faidx_dict[site_chr][1])
    line = int(faidx_dict[site_chr][2])
    size = int(faidx_dict[site_chr][3])
    location = offset + int(start / line) + start - 1
    length = (int(end / line) - int(start / line)) * (size - line) + end - start + 1
    fbuffer.seek(location, 0)
    sequence = fbuffer.read(length).replace('\n','')
    if strand in [-1,'-1','-']:
        sequence = trans(sequence,transform=True,reverse=True)
    return sequence


def trans(seq, transform=True,reverse=False):
    dict1 = {"A":"T",
             "T":"A",
             "G":"C",
             "C":"G",
             "a":"t",
             "t":"a",
             "c":"g",
             "g":"c",
             "n":"n",
             "N":"N"}
    new_seq = ""
    if transform:
        for i in seq:
            new_seq += dict1[i]
    else:
        new_seq = seq
    if reverse:
        new_seq = new_seq[::-1]
    return new_seq


def pickle_make(dict_in,file_out):
    with open(file_out, 'wb') as o:
        pickle.dump(dict_in, o)


def load_pickle(file_in):
    return pickle.load(open(file_in, 'rb'))


def check_overlap(ref, gene):
    flag = 0
    if ref[0] <= gene[0]:
        if ref[1] <= gene[0]:
            flag = 0  # 不挨着ref在左
        elif ref[1] >= gene[1]:
            flag = 2  # ref 包含gene
        else:
            flag = 1  # overlap
    else:
        if ref[0] >= gene[1]:
            flag = 0  # 不挨着ref在右
        elif ref[1] <= gene[1]:
            flag = 2  # gene包含ref
        else:
            flag = 1  # overlap
    return flag


def calc_transcript_length(dict):
    length = 0
    for i in dict:
        length += abs(dict[i][-1] - dict[i][0]) + 1
    return length


def random_weight(weight):
    total = sum(weight.values())
    p = random.randint(1,total)
    t = 0
    for i in weight:
        t += weight[i]
        if t >= p:
            return i


class fa_simulate:
    def __init__(self,database,genome):
        self.db = database
        self.gemone = genome
        self.lengthRange = [2000,10000]
        self.minExonNum = 4
        self.dbWeight = {chrom:len(self.db[chrom]) for chrom in self.db}

    def getFullLength(self,chrom=0):
        database = self.db
        chrom = random_weight(self.dbWeight) if not chrom else chrom
        while chrom in ['MT','Y']:
            chrom = random.choice(list(database))
        gene = random.choice(list(database[chrom]))
        gene_type = database[chrom][gene][0]
        exon_loc_idx = database[chrom][gene][4]  # 所有外显子按位置排序的index
        exon_ts_idx = database[chrom][gene][7]  # 所有转录本各自的外显子index internal_incomplete_known_transcript
        exon_count = database[chrom][gene][6]  # 各转录本外显子个数
        splice_sites = database[chrom][gene][5]  # 全剪切点列表
        transcript = random.choice(list(database[chrom][gene][3]))
        strand = database[chrom][gene][3][transcript][0]
        [t_start, t_end] = database[chrom][gene][3][transcript][1]
        exons = database[chrom][gene][3][transcript][2]
        exon_idx = database[chrom][gene][3][transcript][3]
        length = calc_transcript_length(exons)
        exon_num = len(exons)
        return exons, chrom, gene_type, gene, transcript, strand, length, exon_idx, exon_loc_idx, exon_ts_idx, exon_count, splice_sites,exon_num

    def simOneFullLength(self,chrom=0,minExonNum=0):
        minExonNum = minExonNum if minExonNum else self.minExonNum
        circle_length = 50
        while circle_length:
            circle_length -= 1
            exons, chrom, gene_type, gene, transcript, strand, length, exon_idx, exon_loc_idx, exon_ts_idx, exon_count, splice_sites,exon_num = self.getFullLength(chrom=chrom)
            if not self.exonLengthCheck(exons):
                continue
            if len(exons) >= minExonNum and self.lengthRange[0] <= length <= self.lengthRange[1]:
                return exons, chrom, gene_type, gene, transcript, strand, length, exon_idx, exon_loc_idx, exon_ts_idx, exon_count, splice_sites, exon_num
        else:
            return 0

    def exonLengthCheck(self,exons,minLength=15):
        for i in exons:
            if abs(exons[i][-1]-exons[i][0]) < minLength:
                return 0
        return 1

    def endCutOff(self,exon,strand,side='5',length_cut=0):
        """随机两端裁剪exon,side=['5','3','53']，length_cut: 1:随机两端减去碱基，0:无操作，-1:随机两端增加长度"""
        newExon = exon.copy()
        cut5, cut3 = 0, 0
        cutExon5,cutExon3 = 0,0
        if length_cut != -1:
            cutMax = int(len(newExon)/3)
            cutExon5 = random.randint(1,max(cutMax,1)) if '5' in side else 0
            cutExon3 = random.randint(1,max(cutMax,1)) if '3' in side else 0
            for i in range(cutExon3):
                s = newExon.pop(max(newExon))
                cut3 += abs(s[-1] - s[0])
            for i in range(cutExon5):
                s = newExon.pop(min(newExon))
                cut5 += abs(s[-1] - s[0])
            if length_cut == 1:
                MIN_EXON_LENGTH = 100
                MIN_AFTER_CUT_LENGTH = 70
                if strand in [1,'1','+']:
                    s = newExon[min(newExon)]
                    if s[-1] - s[0] > MIN_EXON_LENGTH:
                        newS = random.randint(s[0]+1,s[-1]-MIN_AFTER_CUT_LENGTH)
                        cut5 += abs(newS - s[0])
                        newExon[min(newExon)][0] = newS
                    s = newExon[max(newExon)]
                    if s[-1] - s[0] > MIN_EXON_LENGTH:
                        newS = random.randint(s[0]+MIN_AFTER_CUT_LENGTH+1, s[-1]-1)
                        cut3 += abs(s[-1] - newS)
                        newExon[max(newExon)][-1] = newS
                else:
                    s = newExon[min(newExon)]
                    if s[-1] - s[0] > MIN_EXON_LENGTH:
                        newS = random.randint(s[0]+MIN_AFTER_CUT_LENGTH+1, s[-1]-1)
                        cut5 += abs(s[-1] - newS)
                        newExon[min(newExon)][-1] = newS
                    s = newExon[max(newExon)]
                    if s[-1] - s[0] > MIN_EXON_LENGTH:
                        newS = random.randint(s[0]+1,s[-1]-MIN_AFTER_CUT_LENGTH)
                        cut3 += abs(newS - s[0])
                        newExon[max(newExon)][0] = newS
            cut3 = -cut3
            cut5 = -cut5
        else:
            add_length = random.randint(300,1000)
            if '5' in side:
                if strand in [1,'1','+']:
                    newExon[min(newExon)][0] = newExon[min(newExon)][0]-add_length
                else:
                    newExon[min(newExon)][-1] = newExon[min(newExon)][-1]+add_length
                cut5 += add_length
            if '3' in side:
                if strand in [1,'1','+']:
                    newExon[max(newExon)][-1] = newExon[max(newExon)][-1]+add_length
                else:
                    newExon[max(newExon)][0] = newExon[max(newExon)][0]-add_length
                cut3 += add_length
        return newExon,cut5,cut3,cutExon5,cutExon3

    def simpleNICmod(self,exon):
        if len(exon) < 4:
            return 0
        subtypesWeight = {"exon_deletion":1,
                          "intron_extension":1,}
        newExon = exon.copy()
        for i in range(random.randint(1,2)):
            st = random_weight(subtypesWeight)
            rm = random.choice(list(newExon)[1:-1])
            if st == "exon_deletion":
                newExon.pop(rm)
            else:
                cutFlag = 0
                temp = newExon[rm]
                for i in newExon:
                    if i == rm:
                        cutFlag = 1
                        continue
                    if cutFlag:
                        newExon[i] = [min(newExon[i]+temp),max(newExon[i]+temp)]
                        newExon.pop(rm)
                        break
                else:
                    print('???')
                    return 0
        return newExon

    def simpleNNCmod(self,exons,splice_sites):
        if not exons:
            return 0
        MIN_DISTANCE_OF_NOVEL_SS = 100
        MAX_DISTANCE_OF_NOVEL_SS = 2000
        ss_list = []
        for i in exons:
            ss_list.append(exons[i][0])
            ss_list.append(exons[i][-1])
        l,r = min(ss_list),max(ss_list)
        ss_list.remove(l)
        ss_list.remove(r)
        cCount = 0
        while 1:
            cCount += 1
            if cCount > 30:
                return 0
            novelSS = random.randint(min(ss_list),max(ss_list))
            if novelSS in splice_sites:
                continue
            for idx, i in enumerate(ss_list):
                ss_temp = i
                distance_temp = abs(novelSS - ss_temp)
                if not idx:
                    distance = distance_temp
                    ss = ss_temp
                    continue
                if distance_temp < distance:
                    distance = distance_temp
                    ss = ss_temp
            if MIN_DISTANCE_OF_NOVEL_SS < distance < MAX_DISTANCE_OF_NOVEL_SS:
                break
        for e in exons:
            if ss in exons[e]:
                exons[e] = [novelSS if i == ss else i for i in exons[e]]
        return exons

    def simpleFUSIONmod(self,exons1,exons2):
        newExons = {}
        count = 0
        for i in exons1:
            count += 1
            newExons[i] = exons1[i]
        for i in exons2:
            count += 1
            newExons[i] = exons2[i]
        return newExons

    def simISM(self,num=5):
        subtypes = {'5ism':'5_ISM_EXON',
                   '3ism':'3_ISM_EXON',
                   '53ism':'INTERNAL_ISM_EXON',
                   '5flank':'5_FLANK',
                   '3flank':'3_FLANK',
                   '53flank':'INTERNAL_FLANK'}
        subtypesWeight = {'5ism': 20,
                          '3ism': 20,
                          '53ism': 20,
                          '5flank': 4,
                          '3flank': 4,
                          '53flank': 4}
        fa,stat,bed = '','',''
        count = 0
        while count < num:
            id = 'ISM_{}'.format(count)
            print("\r{}  ".format(id), end='', flush=1, ),
            info = self.simOneFullLength()
            if not info:
                continue
            exons, chrom, gene_type, gene, transcript, strand, length, exon_idx, exon_loc_idx, exon_ts_idx, exon_count, splice_sites, exon_num = info
            ct = random_weight(subtypesWeight)
            newExons, cut5, cut3, cutExon5, cutExon3 = self.endCutOff(exons, strand, side=ct, length_cut=-1 if 'flank' in ct else 1)
            seq = self.get_sequence(newExons, chrom, strand)
            subtype = subtypes[ct]
            header = self.header_format(id, gene, chrom, strand, calc_transcript_length(exons), len(exons), 'ISM', subtype)
            stat += self.stat_format(id,'ISM',subtype,gene,transcript,chrom,strand,calc_transcript_length(exons),len(exons),cutExon5,cutExon3,cut5,cut3)
            bed += self.bed_format(exons,chrom,strand,id,fusion=0)
            fa += '{}{}\n'.format(header,seq)
            count += 1
        return fa,stat,bed

    def simFSM(self,num=5):
        fa, stat, bed = '', '', ''
        count = 0
        while count < num:
            id = 'FSM_{}'.format(count)
            print("\r{}  ".format(id), end='', flush=1, ),
            info = self.simOneFullLength()
            if not info:
                continue
            exons, chrom, gene_type, gene, transcript, strand, length, exon_idx, exon_loc_idx, exon_ts_idx, exon_count, splice_sites, exon_num = info
            seq = self.get_sequence(exons, chrom, strand)
            header = self.header_format(id, gene, chrom, strand, calc_transcript_length(exons), exon_num, 'FSM', 'FSM')
            stat += self.stat_format(id, 'FSM', 'FSM', gene, transcript, chrom, strand, calc_transcript_length(exons), exon_num, 0, 0, 0, 0)
            bed += self.bed_format(exons, chrom, strand, id, fusion=0)
            fa += '{}{}\n'.format(header, seq)
            count += 1
        return fa,stat,bed

    def simNIC(self,num=5):
        fa, stat, bed = '', '', ''
        subtypesWeight = {'5ism': 10,
                          '3ism': 10,
                          '53ism': 4}
        count = 0
        while count < num:
            id = 'NIC_{}'.format(count)
            print("\r{}  ".format(id), end='', flush=1, ),
            info = self.simOneFullLength()
            if not info:
                continue
            exons, chrom, gene_type, gene, transcript, strand, length, exon_idx, exon_loc_idx, exon_ts_idx, exon_count, splice_sites, exon_num = info
            if calc_transcript_length(exons) > 20000:
                continue
            ct = random_weight(subtypesWeight)
            exons, cut5, cut3, cutExon5, cutExon3 = self.endCutOff(exons, strand, side=ct,length_cut=1)
            exons = self.simpleNICmod(exons)
            if not exons:
                continue
            seq = self.get_sequence(exons, chrom, strand)
            header = self.header_format(id, gene, chrom, strand, calc_transcript_length(exons), len(exons), 'NIC', 'NIC')
            stat += self.stat_format(id, 'NIC', 'NIC', gene, transcript, chrom, strand, calc_transcript_length(exons), len(exons), cutExon5,cutExon3,cut5,cut3)
            bed += self.bed_format(exons, chrom, strand, id, fusion=0)
            fa += '{}{}\n'.format(header, seq)
            count += 1
        return fa,stat,bed

    def simNNC(self,num=5):
        fa, stat, bed = '', '', ''
        subtypesWeight = {'5ism': 10,
                          '3ism': 10,
                          '53ism': 4}
        count = 0
        while count < num:
            id = 'NNC_{}'.format(count)
            print("\r{}  ".format(id), end='', flush=1, ),
            info = self.simOneFullLength()
            if not info:
                continue
            exons, chrom, gene_type, gene, transcript, strand, length, exon_idx, exon_loc_idx, exon_ts_idx, exon_count, splice_sites, exon_num = info
            ct = random_weight(subtypesWeight)
            exons, cut5, cut3, cutExon5, cutExon3 = self.endCutOff(exons, strand, side=ct, length_cut=1)
            for i in range(random.randint(1,2)):
                exons = self.simpleNNCmod(exons,splice_sites)
            if not exons:
                continue
            seq = self.get_sequence(exons, chrom, strand)
            header = self.header_format(id, gene, chrom, strand, calc_transcript_length(exons), len(exons), 'NNC', 'NNC')
            stat += self.stat_format(id, 'NNC', 'NNC', gene, transcript, chrom, strand, calc_transcript_length(exons), len(exons), cutExon5,cutExon3,cut5,cut3)
            bed += self.bed_format(exons, chrom, strand, id, fusion=0)
            fa += '{}{}\n'.format(header, seq)
            count += 1
        return fa,stat,bed

    def simFUSION(self,num=5):
        def idMerge(id1,id2):
            return '{}|{}'.format(id1,id2)

        fa, stat, bed = '', '', ''
        subtypesWeight = {'5ism': 10,
                          '3ism': 10,
                          '53ism': 10}
        count = 0
        while count < num:
            id = 'FUSION_{}'.format(count)
            print("\r{}  ".format(id), end='', flush=1, ),
            info1 = self.simOneFullLength()
            info2 = self.simOneFullLength()
            if not info1 or not info2:
                continue
            exons1, chrom1, gene_type1, gene1, transcript1, strand1, length1, exon_idx1, exon_loc_idx1, exon_ts_idx1, exon_count1, splice_sites1, exon_num1 = info1
            exons2, chrom2, gene_type2, gene2, transcript2, strand2, length2, exon_idx2, exon_loc_idx2, exon_ts_idx2, exon_count2, splice_sites2, exon_num2 = info2
            exons1, cut51, cut31, cutExon51, cutExon31 = self.endCutOff(exons1, strand1, side=random_weight(subtypesWeight), length_cut=1)
            exons2, cut52, cut32, cutExon52, cutExon32 = self.endCutOff(exons2, strand2, side=random_weight(subtypesWeight), length_cut=1)
            #newExon = self.simpleFUSIONmod(exons1,exons2)
            seq1 = self.get_sequence(exons1, chrom1, strand1)
            seq2 = self.get_sequence(exons2, chrom2, strand2)
            seq = seq1+seq2
            header = self.header_format(id, idMerge(gene1,gene2), idMerge(chrom1,chrom2), idMerge(strand1,strand2), idMerge(calc_transcript_length(exons1),calc_transcript_length(exons2)), idMerge(len(exons1),len(exons2)), 'FUSION', 'FUSION')
            stat += self.stat_format(id, 'FUSION', 'FUSION', idMerge(gene1,gene2), idMerge(transcript1,transcript2), idMerge(chrom1,chrom2), idMerge(strand1,strand2), idMerge(length1,length2), idMerge(len(exons1),len(exons2)), 'NA','NA','NA','NA')
            bed += self.bed_format([exons1,exons2], [chrom1,chrom2], [strand1,strand2], id, fusion=1)
            fa += '{}{}\n'.format(header, seq)
            count += 1
        return fa,stat,bed

    def get_sequence(self,exons,chrom,strand):
        seq = ''
        for i in exons:
            seq += refseq_extract(str(chrom), 0, exons[i], fa=self.gemone,strand=strand)
        return seq

    def bed_format(self,exons,chrom,strand,id,fusion=0):
        bed_txt = ''
        if fusion:
            exons,exons2 = exons
            chrom,chrom2 = chrom
            strand,strand2 = strand
        for i in exons:
            bed_txt += '{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrom,exons[i][0],exons[i][1],id,'*',strand)
        if fusion:
            for i in exons2:
                bed_txt += '{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrom2, exons2[i][0], exons2[i][1], id, '*', strand2)
        return bed_txt

    def header_format(self,id,gene,chrom,strand,length,exon_num,classification,subtype):
        return '>{} Gene={} chrom={} stand={} length={} exon_num={} type={} subtype={}\n'.format(id,gene,chrom,strand,length,exon_num,classification,subtype)

    def stat_format(self,id,classification,subtype,gene,transcript,chrom,strand,length,exon_num,cutExon5,cutExon3,cut5,cut3):
        return '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(id,classification,subtype,gene,transcript,chrom,strand,length,exon_num,cutExon5,cutExon3,cut5,cut3)

    def oWrite(self,context,fileopen):
        fa,stat,bed = context
        f,s,b = fileopen
        f.write(fa)
        s.write(stat)
        b.write(bed)

    def run(self,output,num=5):
        with open(output+'.fa','w') as f, open(output+'.stat','w') as s, open(output+'.bed','w') as b:
            s.write('ID\tType\tSubtype\tGene\tTranscript\tChrom\tStrand\tLength\tExonNum\t5`exonDeletion\t3`exonDeletion\t5`lengthCut\t3`lengthCut\n')
            self.oWrite(self.simFSM(num=num),[f,s,b])
            self.oWrite(self.simISM(num=num),[f,s,b])
            self.oWrite(self.simNIC(num=num),[f,s,b])
            self.oWrite(self.simNNC(num=num),[f,s,b])
            self.oWrite(self.simFUSION(num=num),[f,s,b])


def options(argv):
    from argparse import ArgumentParser as AP
    usages = "python3 {} -d db -o fa -n 100".format(argv[0])
    p = AP(usage=usages)
    p.add_argument("-g", dest="genome", metavar="genome fa", help="reference genome fa", type=str, default=100)
    p.add_argument("-d", dest="db", metavar="db", help="gtfDB file",required=True)
    p.add_argument("-o", dest="fa", metavar="fa_out", help="Fasta file output",required=True)
    p.add_argument("-n", dest="num", metavar="[int]", help="number of transcript simulation, default: 100", type=int, default=100)
    if len(argv) == 1:
        p.print_help()
        exit(1)
    return p.parse_args(argv[1:])


if __name__ == '__main__':
    args = options(sys.argv)
    db = load_pickle(args.db)
    A = fa_simulate(db,args.genome)
    A.run(args.fa,args.num)
