# TransAnnot 使用手册

[English](https://github.com/captorr/TransAnnot/blob/master/README.md) | 中文文档 | [gitee](https://gitee.com/captor1/TransAnnot)

TransAnnot 是一款用于三代全长转录本测序结果注释的工具。

## 软件安装

	git clone https://github.com/captorr/TransAnnot.git

## 环境依赖

 * [HISAT2](https://github.com/DaehwanKimLab/hisat2)/[MINIMAP2](https://github.com/lh3/minimap2)/[GMAP](http://research-pub.gene.com/gmap/)中至少一款软件
 * samtools
 * python3

## FAST RUN

	python TransAnnot.py -f [fasta] -g [genome fasta] -o [output directory] -a [annot gtf] -p [process] --use_minimap2 1 --use_hisat2 [hisat2 index]

## 软件运行

1. 软件所需参数可由以下3种方式传入，序号大的参数会覆盖低序参数：
	1. 软件目录下的基础配置文件TransAnnot.Config
	2. 参数 `-c [config]` 导入的临时配置文件
	3. 命令行输入参数

2. 首次运行先调整基础配置文件[TransAnnot.Config](https://github.com/captorr/TransAnnot/blob/master/TransAnnot.Config), 主要设置以后运行时不需要频繁改动的参数：
	* 相关软件的调用路径
	* 默认使用的映射软件、映射软件的INDEX路径
	* 参考基因组(FASTA)、参考基因组注释文件（GTF）、默认PROCESS数

3. 设置好基础配置文件后，运行某一样本时仅需提供**输入fasta文件**和**输出文件夹路径**即可。可由`-c config`设置，也可由命令行输入`-f [fatsa] -o [output]`

4. 软件支持reads\transcript\gene水平的表达量统计，需要通过`--tpm`参数或`-c config`输入基于**reads ID**的表达量列表

5. 软件支持多样本的联合分析\表达量统计，可通过单独运行[TransAnnotMerge.py]()实现，具体操作说明见后续。

## 示例

这里提供一个[流程示例](https://github.com/captorr/TransAnnot/tree/master/example)

## 运行结果

如未改动基础配置文件[TransAnnot.Config](https://github.com/captorr/TransAnnot/blob/master/TransAnnot.Config)中的后缀参数，输出结果将有：

* `{sample_id}.annot.bed` 所有reads的bed格式结果

* `{sample_id}.annot.stat`  所有reads的注释统计结果

* `{sample_id}.annot.db.pickle`  用于后续可视化分析

* `{sample_id}.annot.cluster.gene`  有注释reads的gene水平聚类结果

* `{sample_id}.annot.cluster.transcript`  有注释reads的transcript水平聚类结果

* `{sample_id}.annot.cluster.reads`  有注释reads的聚类信息

* `{sample_id}.annot.junction` 剪切点信息

* `{sample_id}.annot.multiAnno` 多重注释信息

### {sample_id}.annot.stat 各列信息：

* `ID`： reads ID
* `Classification`： classificatioin of reads
* `Subtype`: subtype of reads
* `Gene`: gene annotation or region in genome[chr1:100000-100500]
* `Transcript`: transcript annotation
* `Chrom`: chromosome
* `Strand`: strand
* `Seq_length`: reads length
* `Seq_exon`: exon number of reads
* `Ref_length`: length of transcript annotation
* `Ref_exon_num`: exon number of transcript annotatioin
* `diff_to_gene_start`: 5` site difference of reads and annotation gene in reference genome
* `diff_to_gene_end`: 3` site difference of reads and annotation gene in reference genome
* `diff_to_transcript_start`: 5` site difference of reads and annotation transcript in reference genome
* `diff_to_transcript_end`: 3` site difference of reads and annotation transcript in reference genome
* `exon_miss_to_transcript_start`: number of exon missed in 5` site between reads and transcript annotation
* `exon_miss_to_transcript_end`: number of exon missed in 3` site between reads and transcript annotation

## TransAnnot.Config

配置文件分为**基础配置文件**和**临时配置文件**两种。

基础配置文件：需放置在软件目录下，由软件自动检查，无需参数导入。该配置文件仅允许修改参数，不允许删除或增加行。该配置文件中的参数，优先级最低，通常用来设置不需要频繁修改的参数。不确定的参数可以不填。

临时配置文件：可存放于任意位置，由`-c`参数导入。该配置文件可任意增删修改，参数优先级高于基础配置文件，如二者有冲突则使用临时配置文件中的参数。通常用来单独配置随样本变化的参数。不随样本或项目变动的参数可删除以保持视觉上的清爽。

配置文件相关参数释义：

* `FASTA`: `[path]`，输入文件，fa格式
* `OUTPUT_DIR`: `[path]`，输出文件夹
* `GENOME_FA`: `[path]`，参考基因组fa文件
* `GTF_ANNOTATION`:`[path]`，注释gtf文件
* `PROCESS`: `[int]`,并行进程数
* `SAMPLE_UNIQUE_NAME`:`[string]`,样本名，后续如果有多样本合并时，该名称需唯一
* `SAMTOOLS`:`[path]`，samtools的运行路径
* `USE_HISAT2`: `[int]`，是否使用Hisat2，0代表不用，1代表使用
* `HISAT2`: `[path]`,Hisat2的运行路径
* `HISAT2_INDEX`: `[path]`，Hisat2 index的路径，由`hisat2-build`生成
* `USE_MINIMAP2`: `[int]`，是否使用Minimap2，0代表不用，1代表使用
* `MINIMAP2`: `[path]`，Minimap2的运行路径
* `USE_GMAP`: `[int]`，是否使用GMAP，0代表不用，1代表使用
* `GMAP`: `[path]`，GMAP的运行路径
* `GMAP_INDEX`: `[path]`,GMAP index的路径，由`gmap_build`生成
* `TPM_LIST`: `[path]`，Isoform的表达量文件，可选参数，如有则结果将包含定量信息，空白则不包含
* `READ_LENGTH`: `[int]`，二代映射fa的read length，默认100
* `READ_OVERLAP`: `[int]`，二代映射fa的read overlap，默认80
* `MIN_READ_LENGTH`: `[int]`，二代映射fa的最短长度，默认30
* `REPORT_IN_RUNNING`: `[int]`，是否在运行中实时打印进程报告，0表示静默，1表示实时输出


## TransAnnotMerge

用于多样本的联合分析

举个例子，当样本1中发现了geneA的一个新转录本，样本2中同样发现了geneA的一个新转录本，然而他们在各自的结果文件中都叫做`geneA_NNC_1`，此时如何判断两个`geneA_NNC_1`是否真的具有相同结构、是否真为同一种新转录本？

另一个例子，我想得到多个样本中所有基因\转录本的表达量矩阵，然而各样本中都存在独特的、其他样本中没有的基因\转录本，该怎么办？

此时可交由`TranAnnotMerge`处理。`TranAnnotMerge`将所有样本中的Isoform重新编排聚类，对新转录本、基因统一分配新的ID，输出新的聚类结果，当提供表达量数据时还能输出全样本的表达量矩阵

### TransAnnotMerge运行命令

`python TranAnnotMerge.py -c MergeConfig -o outputdir -m [TPM/FLC/None]`

包含3个参数：

* `-c`： Merge Config，由4列内容组成，分别为`{sample_id}`，`{sample_id}.annot.stat`文件路径，`{sample_id}.annot.bed`文件路径，`{sample_id}.annot.db.pickle`文件路径。形如

| #sample | stat | bed | db |
| ------- | ---- | --- | -- |
| ------- | ---- | --- | -- |

* `-o`: 融合结果的输出目录
* `-m`：选定表达量矩阵中的数值表示(需要前序分析中导入了表达量文件)，FLC指 full length count，不加此参数则无表达量矩阵

### TransAnnotMerge流程演示

1. 从fa中提取Isoform表达量文件(exp)：`python fa2exp.py -f [fa] -o [exp]`
2. 将各个样本的exp文件通过**配置文件**中的**TPM_LIST**参数或用`-tpm`导入，然后分别运行TransAnnot
3. 汇总各样本结果写入`MergeConfig`文件，运行TranAnnotMerge: `python TranAnnotMerge.py -c MergeConfig -o outputdir -m TPM`

### TransAnnotMerge结果文件

* `{sample_id}.reads.exp`： 各样本read水平结果文件
* `{sample_id}.transcript.exp`： 各样本transcript水平结果文件
* `{sample_id}.gene.exp`： 各样本gene水平结果文件
* `gene.exp`： 全样本gene水平表达量矩阵
* `transcript.exp`： 全样本transcript水平表达量矩阵
* `merge.db.pickle`： 全样本分析数据，用于后续可视化分析


## TransAnnot的转录本分类方式

基于SQANTI的分类做了一些改动，目前有以下几类：

### 转录本水平分类

* `FSM`: 已知转录本，full splice site match
* `ISM`: 不完整的已知转录本，incomplete splice site match
* `NIC`: 剪切点已知的新转录本，novel in catalog
* `NNC`: 含有新剪切点的新转录本，novel not in catalog
* `GENIC`: 基因内的，genic
* `INTERGENIC`: 基因间的，intergenic
* `FUSION`: 融合，fusion
* `UNKNOWN`： unknown

### 外显子水平分类

* `KE`: known exon
* `LEKE`: 左边缘外显子, left end known exon
* `REKE`: 右边缘外显子, right end known exon
* `NEKSLE`: 左边缘外显子融合, novel exon with known splice site in left end exon and has the unique region overlap with at least two known exons
* `NEKSRE`: 右边缘外显子融合, novel exon known splice site in right end exon and has the unique region overlap with at least two known exons
* `IE`: 内含子保留, intron retention: two known splice sites from the same transcript's sequential exon
* `NEDT`:两端剪切点都已知的新外显子,  novel exon with two known splice sites from different transcript
* `NELS`: 左端新剪切点, novel exon with novel left splice site
* `NERS`: 右端新剪切点, novel exon with novel right splice site
* `LEE`: 左端不对齐延伸外显子, left exon_extension： the novel splice site in the left end of the exon which is longer than any exons overlap with it
* `REE`: 右端不对齐延伸外显子, right exon_extension: the novel splice site in the right end of the exon which is longer than any exons overlap with it
* `NEDS`: 双端新剪切点, novel exon:double novel splice sites overlap with at least one known exon
* `NEIG`: 双端新剪切点基因内无注释的新外显子, novel exon inner-gene：novel exon inside the gene and without any overlap with known exon
* `NEOG`: 双端新剪切点基因外无注释的新外显子, novel exon inter-gene：novel exon outside the gene
* `NELE`: 左边缘外显子新剪切点, novel exon with novel splice site in the far left exon
* `NERE`: 右边缘外显子新剪切点, novel exon with novel splice site in the far right exon
* `MDNS`: 单外显子转录本双端不对齐, monoexon with double novel splice sites

## 其他开发工具

TransAnnot开发过程中做的实用工具，主要服务于TransAnnot本体，因此仅做了些基本功能

### LRSIM

LRSIM (Long Reads SIMulation), 用于模拟生成全长转录本,可批量生成FSM、ISM、NIC、NNC、FUSION五类模拟转录本数据。仅生成纯净序列，不做ONT/PB随机扰动修饰。

### TransAnnotViewer

基于TransAnnot结果的可视化工具，windows桌面端, 可对结果绘图
![transannotviewer](https://github.com/captorr/TransAnnot/raw/master/static/transannotviewer.png)
![transannotviewer](https://github.com/captorr/TransAnnot/raw/master/static/ARH.png)

### TransAnnotReport

基于TransAnnot结果的统计工具，生成一份PDF报告，功能简陋。

## Citation

TransAnnot的主要功能已集成到[TAGET](https://github.com/gx-health/TAGET)中，如此工具有助于您的研究，请引用以下文章

Xia, Y., Jin, Z., Zhang, C. et al. TAGET: a toolkit for analyzing full-length transcripts from long-read sequencing. Nat Commun 14, 5935 (2023). https://doi.org/10.1038/s41467-023-41649-0