## EXAMPLE

### step1 generate simulate fasta

#### make GTF_DB
	
	# about 1 min
	cd /your_path/TransAnnot
	python gtf2db.py -i [gtf_file] -o [gtf_db] -s ensembl

#### generate fasta

Here we generate two sets fasta to test the merge pipline, you can generate only one if no required.
	
	# about 1min per set
	cd /your_path/TransAnnot/LRSIM
	python LRSIM.py -g [reference genome fa] -d [gtf_db] -o [sim_fa1] -n 50 -r
	python LRSIM.py -g [reference genome fa] -d [gtf_db] -o [sim_fa2] -n 50 -r

#### extract READ_NUM/FLC

You can skip this step if you don't need Expression matrix.

	# about several sec	
	cd /your_path/TransAnnot
	python fa2exp.py -f [sim_fa1] -o [sim_fa1_exp]
	python fa2exp.py -f [sim_fa2] -o [sim_fa1_exp]


### step2 running TransAnnot

#### edit BASE CONFIG

Edit BASE CONFIG file [TransAnnot.Config](https://github.com/captorr/TransAnnot/blob/master/demo/TransAnnot.Config), you can follow the example file in this directory.

#### running by BASE CONFIG

Run TransAnnot by these command with the edited BASE CONFIG. 

If you don't need Expression matrix, don't enter the parameter `--tpm`
	
	# about 4 min per set
	cd /your_path/TransAnnot
	python TransAnnot.py -f [sim_fa1] -o [output_dir1] -n sim1 -p 8 --tpm [sim_fa1_exp]
	python TransAnnot.py -f [sim_fa2] -o [output_dir2] -n sim2 -p 8 --tpm [sim_fa2_exp]

#### fast run

You can also run TransAnnot without BASE CONFIG by these command.
Same as before, don't enter the parameter `--tpm` if you don't need Expression matrix

	cd /your_path/TransAnnot
	python TransAnnot.py -f [sim_fa1] -o [output_dir1] -g [reference_genome] -a [annotation_gtf] -p 8 --use_minimap2 1 --use_hisat2 [hisat2_index] -n sim1 --tpm [sim_fa1_exp]
	python TransAnnot.py -f [sim_fa2] -o [output_dir2] -g [reference_genome] -a [annotation_gtf] -p 8 --use_minimap2 1 --use_hisat2 [hisat2_index] -n sim2 --tpm [sim_fa1_exp]

The work of TransAnnot is done here.

### step3 running TransAnnotMerge

This step is for multi-samples, you can skip this step if you have only one sample

#### edit Merge.config

edit [merge.config] file in format，see [merge.config](https://github.com/captorr/TransAnnot/blob/master/demo/merge.config)

#### run TransAnnotMerge

you can skip `-m` if you don't need Expression matrix

	cd /your_path/TransAnnot
	python TransAnnotMerge.py -c [merge.config] -o [output_path] -m TPM

### step4 TransAnnotViewer

Both of output directory of single sample and TransAnnotMerge have **[*.db.pickle]** file，which was used for visualization of TransAnnotViewer.
This step required windows desktop.

	cd /your_path/TransAnnot/TransAnnotViewer
	python TransAnnotViewer.py -i [db.pickle]