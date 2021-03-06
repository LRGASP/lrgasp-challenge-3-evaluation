# Challenge 3 Evaluation

The evaluation of submitted transcript models will be done using [SQANTI3](https://github.com/ConesaLab/SQANTI3) descriptors and number of LRGASP-agreed evaluation metrics. The general procedure to assess these evaluation metrics on your data is to run **sqanti3_lrgasp.challenge3.py**, an adapted version of the original code that will generate automatically an HTML report with the results of the evaluation for challenge 3. Please, download/clone the entire [sqanti3_evaluation](https://github.com/LRGASP/lrgasp-challenge-3-evaluation.git) directory, including the **utilities** subfolder.

Challenge 3 of LRGASP was thought to evaluate the performance of different pipelines to recontruct a transcriptome without using previous annotation nor a reference genome (except for the *long and genome* and *free-style* categories, for which it is allowed to use a genome). However, for its evaluation, we will use a *de novo* ONT-based genome available [here](https://www.synapse.org/#!Synapse:syn25683366) to obtain information about the submitted transcriptomes.

## Setting up the environment

In order to install all the dependencies needed by **sqanti3_lrgasp.challenge3.py**, please use the [YML](https://github.com/LRGASP/lrgasp-challenge-3-evaluation/blob/main/sqanti3_lrgasp.yml) file to build a conda environment. 

```
conda env create -f sqanti3_lrgasp.yml
source activate sqanti3_lrgasp
```

SQANTI3 also takes advantage of some scripts of [cDNA_Cupcake](https://github.com/Magdoll/cDNA_Cupcake/wiki#install). Please install it with the *sqanti3_lrgasp* environment activated and add `cDNA_Cupcake`and `cDNA_Cupcake/sequence` to your `PYTHONPATH`.

```
(sqanti3_lrgasp)$ git clone https://github.com/Magdoll/cDNA_Cupcake.git
(sqanti3_lrgasp)$ cd cDNA_Cupcake
(sqanti3_lrgasp)$ python setup.py build
(sqanti3_lrgasp)$ python setup.py install

(sqanti3_lrgasp)$ export PYTHONPATH=$PYTHONPATH:<path_to>/cDNA_Cupcake/sequence/
(sqanti3_lrgasp)$ export PYTHONPATH=$PYTHONPATH:<path_to>/cDNA_Cupcake/

```

Remember to activate the *sqanti3_lrgasp* environment and setting up the `PYTHONPATH` every time you are running the evaluation.

## Run SQANTI3

When running [SQANTI3](https://github.com/ConesaLab/SQANTI3), your transcript-models are usually compared against some kind of genome annotation, but in this case, that's not a possibility since we are replicating the situation in which there isn't any type of information beyond sequencing data. For Challenge 3 the only annotation file that will be provided to SQANTI3 will be the spike-ins information. [Here you can find that GTF](https://www.synapse.org/#!Synapse:syn25683630). It contains the annotation of Lexogen's SIRVs-Set 4.

LRGASP will be using Illumina junction coverage to evaluate your transcript models. We therefore recommend you run **sqanti3_lrgasp.challenge3.py** enabling this analysis. To do so:

-   **SJ coverage**:  As SJ information is dependent on the sample being analyzed, it is necessary to run previously STAR to map the Illumina reads against the genome and identify possible splice-junctions using the `--twopassMode`. Then, the resulting _*SJ.oyut.tab_ file can be input to **sqanti3_lrgasp.challenge3.py** with the parameter `-c`. This is an example of how we normally run STAR for this SJ detection:

1. Create a genome index with STAR without providing the reference annotation. We don't want to bias the SJ-detection towards the annotated splice sites.

```
STAR --runThreadN <num_threads> --runMode genomeGenerate --genomeDir <star_index> --genomeFastaFiles <reference_genome_FASTA> --outTmpDir <index_dir_tmp> 
```

2. Map short-reads using `--twopassMode`

```
STAR --runThreadN <num_threads> --genomeDir <star_index> --readFilesIn <read1> <read2> --outFileNamePrefix <output_prefix> --alignSJoverhangMin 8  --alignSJDBoverhangMin 1 --outFilterType BySJout --outSAMunmapped Within --outFilterMultimapNmax 20 --outFilterMismatchNoverLmax 0.04 --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --sjdbScore 1 --genomeLoad NoSharedMemory --outSAMtype BAM SortedByCoordinate --twopassMode Basic
```

It is also neccessary to provide metadata files in JSON format. To complete a submission it is mandatory to create experiment and entry JSON files and those are the ones that `sqanti3_lrgasp.challenge1.py` is expecting. [Here](https://lrgasp.github.io/lrgasp-submissions/docs/metadata.html) you can find more information about metadata files and some templates for them.

If you don't have those JSON files ready yet, now `--experiment_json` and `--entry_json` are optional arguments. If they are not defined, inside `utilities/` folder there are a couple of *dummy* files that will do the job. Please, make sure that `species` field of the **experiment_dummy.json** is right. Manatee samples used a different set of spike-ins compared to mouse samples. In consequence, if the `species` are not set properly, it will result into a wrong evaluation of SIRVs.


### Example

This is an example of how to run the **sqanti3_lrgasp.challenge3.py** script:

```
python sqanti3_lrgasp.challenge3.py manatee_submitted.fasta lrgasp_sirv4.gtf lrgasp_manatee_sirv1.fasta \
	--experiment_json manatee_experiment.json --entry_json manatee_entry.json \
	-c my_test.SJ.out.tab -d /my_output/directory -o manatee_submission_test
```

This will generate in `/my_output/directory` the two main SQANTI3 outputs (`*_classification.txt`and `*_junctions.txt`) and a HTML file that will be called, in this case, `manatee_submission_test_Evaluation_report.html`.
