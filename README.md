#QueryFuse
-QueryFuse is software package for gene fusion detection using RNA-Seq data on a list of genes of interested. By using the small gene list, QueryFuse can give report in a fast, accurate and sensitive way, which can also help users focus on the genes they mostly care about. The software use alignment output from general RNA-Seq aligners as the input to avoid unnecessary realignment and save time. Then it uses local aligner to cluster discordant paired end reads and locate the fusion boundary. The software integrates filters and ranking methods to reduce false positive and prioritize important events. Fusion events that have only spanning reads1 will be reported in the second category. Alignment graph with reference template and all the related supporting reads will be generated automatically.
-1Spanning reads: We define spanning reads as pair-end-aligned reads that have the fusion boundaries in the gap between the paired ends.


##Publications


##RNA-Seq Datasets
-Publicly available datasets are all described in the publication. The TL-03 sample location? (need to ask bejoern)
-A simulated test set (with 100 splitting reads and 26 spanning reads for each fusion) is in the simulation-test-data folder. (accepted_hits.bam, unmapped.bam and 1 simulated fusion list) Also, parameter_file.txt for this dataset is provided.(available in: http://smonti.bumc.bu.edu/~montilab/zoho/QueryFuse/simulated_data_example/)



##Setup
-To run QueryFuse, you need to have anaconda installed in your Linux system.
-Get the QueryFuse code
-Download and untar/gunzip the QueryFuse code package to a directory. No further installation is needed for the code itself, but anaconda environment is needed. "dos2unix" (or similar commands) are recommended, in order to confirm all the files are in the correct format for your computer. Now, QueryFuse can only run under Linux/Unix, which supports shell script, because QueryFuse is written in shell, Python, Perl and R.
-Setup running environment using anaconda
-In each QueryFuse version, there will be a file called "anaconda_env_setting_readme.txt". Please follow it to set up the running environment. There are external tools required for running QueryFuse and they will all be setup by anaconda.



##Reference Dataset (available in: http://smonti.bumc.bu.edu/~montilab/zoho/QueryFuse/reference_files/)
-1.	Genome.fa and build the bowtie index. For bowtie website. [Note: the fa file must use "chr*" as chromosome rather than use just 1-22, X, Y and M.] The hg19.fa I used can be also downloaded in the reference folder. (You need to build bowtie index yourself based on the fa by using "bowtie2-build".)
-2.	Gene annotation. They can be built in bio-mart. Columns to include and following this order. For example, the one for hg19 (hg19_whole_gene_list.bed in the reference folder). In order to have better extension, I generated a 5k bp extended version on both sides of each gene. (hg19_whole_gene_list_5k_expanded.bed in the reference folder). Because the chrM is shortly, need to pay attention that not to extent the boundary over chrM's length. However, generally, chrM is not considered in fusion detection and should be filtered out.)



##Input data format
-QueryFuse now takes only paired end aligned bam file from aligners as input. To convert sam to bam, samtools command "samtools view" can be used.



##How to run
-Running QueryFuse is as simple as running a single command line by just typing "python QF_path/QF_multi_query_warpper.py parameter_file". The parameter file should include all the parameters listed below. (But you can also type in the parameter in the command line without changing the parameter file.)
-[Details of all parameters: please see the QF_para_file_example.txt in each QueryFuse version.]
*Note: the output directory should be different from the directory containing the input bam file.*



##Output
###Final output
-For a specific sample, all the most important outputs are in the "results" folder and each query gene will have its own subfolder. The "whole_fusion_sum_filtered.txt" file is the key one with score and ranked result after filtering. The "whole_fusion_sum_all.txt" file contains all the detect fusions before filtering and has scores of features but no ranks. If you want to have a look at the detail graph of ranked fusion events, you can find it in the "fusion_supporting_graph_final" folder by the fusion details.
-There are three other folders in a run: "logs" folder with all running logs; "bams" folder with shared preprocess files for all queries; "intermedias" folder with query specific intermediate files for each query. 


###Identification
-PARTNER_GENE_NAME: the name of the partner gene, which can be ambiguous.
-PARTNER_GENE_ID: the ENSEMBL ID of the partner gene, which is unique.
-QUERY_GENE_NAME: the name of the query gene.
-QUERY_GENE_ID: the ENSEMBL ID of the query gene.

###Evidence 
-SPLIT_NUM: number of splitting reads supporting this fusion event. If it is fusion event with multiple alignments, the sum of all splitting reads is divided by the number of multiple locations.
-SPAN_NUM: number of spanning reads supporting this fusion event.
-SUPPORT_SUM_NUM: Total supporting read number of this fusion event. It is the sum of SPLIT_NUM and SPAN_NUM.
-SPLIT_PVAL: the value to indicate how well the splitting reads spread around the breakpoint. (The bigger, the better.) It is the p-value from the KS-test.
-SHIFT_RANGE: the length of shared sequence at the breakpoints of the pair of partners. Shifting the breakpoint in this region will not affect the fusion sequence. Biologically, this is locations that an enzyme can bind on and the double-strand DNA break or splicing event can happen at.
-DINUCLEOTIDE_ENTROPY: entropy of 16 possible dinucleotides around the breakpoint.
-MULTI-ALIGN_NUM: number of possible alignment locations for this fusion reference.
-RANK_SCORE: The sum of ranks of all previous features.


###Annotation
-CHR_PARTNER_GENE: the chromosome of the partner gene.
-BREAKPOINT_PARTNER_GENE: the breakpoint location of the partner gene.
-DIRECTION_PARTNER_GENE: the chromosomal connection direction of the partner gene. ("F" means forward direction, which means it is 3' end on plus strand or 5' end on the minus strand and the detected arm is on the left (side with smaller location number than the breakpoint) of a chromosome. "R" means reverse direction, which means it is 5' end on the plus strand or 3' end on the minus strand and the detected arm is on the right (side with bigger location number than the breakpoint) of a chromosome.
-CHR_QUERY_GENE: the chromosome of the query gene
-BREAKPOINT_QUERY_GENE: the breakpoint location of the query gene.
-DIRECTION_QUERY_GENE: the chromosomal connection direction of the query gene.
