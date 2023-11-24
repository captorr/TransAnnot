## 流程测试

[English](https://github.com/captorr/TransAnnot/blob/master/example/readme.md) | 中文文档

### step1 生成模拟数据

#### 生成GTF_DB
	
	# 大约耗时1分钟
	cd /your_path/TransAnnot
	python gtf2db.py -i [gtf_file] -o [gtf_db] -s ensembl

#### 生成模拟数据

生成两套模拟数据
	
	# 每组大约耗时1分钟，可并行
	cd /your_path/TransAnnot/LRSIM
	python LRSIM.py -g [reference genome fa] -d [gtf_db] -o [sim_fa1] -n 50 -r
	python LRSIM.py -g [reference genome fa] -d [gtf_db] -o [sim_fa2] -n 50 -r

#### 提取READ_NUM/FLC

同样提取两个数据(不需要表达量的可以省略此步)

	# 耗时越几秒钟	
	cd /your_path/TransAnnot
	python fa2exp.py -f [sim_fa1] -o [sim_fa1_exp]
	python fa2exp.py -f [sim_fa2] -o [sim_fa1_exp]


### step2 运行TransAnnot

#### 配置基础config

配置TransAnnot.Config文件,填入依赖软件环境变量或运行路径，填入参考基因组、注释GTF等文件路径，注释gtf默认使用ensembl版本

#### 命令行运行

基础配置文件中写好依赖软件的前提下，可运行以下命令处理两个数据（不需要表达量的可以省略`--tpm`参数）
	
	#每组耗时越4分钟，可并行
	cd /your_path/TransAnnot
	python TransAnnot.py -f [sim_fa1] -o [output_dir1] -n sim1 -p 8 --tpm [sim_fa1_exp]
	python TransAnnot.py -f [sim_fa2] -o [output_dir2] -n sim2 -p 8 --tpm [sim_fa2_exp]

#### 直接运行

不填写基础配置文件的情况下，也可由以下命令直接运行（不需要表达量的可以省略`--tpm`参数）

	cd /your_path/TransAnnot
	python TransAnnot.py -f [sim_fa1] -o [output_dir1] -g [reference_genome] -a [annotation_gtf] -p 8 --use_minimap2 1 --use_hisat2 [hisat2_index] -n sim1 --tpm [sim_fa1_exp]
	python TransAnnot.py -f [sim_fa2] -o [output_dir2] -g [reference_genome] -a [annotation_gtf] -p 8 --use_minimap2 1 --use_hisat2 [hisat2_index] -n sim2 --tpm [sim_fa1_exp]

主程序至此已完成

### step3 运行TransAnnotMerge

如有合并多样本的需求，可执行此步完成。单样本可跳过此阶段。

#### 填写Merge.config

按文件格式填写[merge.config]，见示例文件[merge.config](https://github.com/captorr/TransAnnot/blob/master/example/merge.config)

#### 运行TransAnnotMerge

不需要表达量的可以省略`-m`参数

	cd /your_path/TransAnnot
	python TransAnnotMerge.py -c [merge.config] -o [output_path] -m TPM

### step4 TransAnnotViewer可视化

单样本或TransAnnotMerge的输出文件夹中均有 **[*.db.pickle]** 文件，此文件用于接下来的转录本的可视化，这一步需要在windows等具备窗口视图的系统中运行。

	cd /your_path/TransAnnot/TransAnnotViewer
	python TransAnnotViewer.py -i [db.pickle]