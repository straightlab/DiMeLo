# DiMeLo
# Scripts for analysis of DiMeLo data from in vitro and in situ DiMeLo-seq

#04.04.23

#This folder contains scripts for DiMeLo-seq analyses split into invitro and insitu analyses, and v1 (until Nature Methods 2022 publication) and v2 (after Nature Methods 2022 publication)

##
# In vitro analyses:
##

#The following folders and files correspond to in vitro analyses

# Workflow_for _in_vitro_DiMeLo_analyses_KS040323.sh

###The following code is for analyses of in vitro DiMeLo sequencing data output by ONT in the form of fast5 files. The final output is in the form of bam file for each barcode, with reads filtered by length.
###The code below corresponds to the following workflow

#Input (from nanopore) - fast5 (filename.fast5)
#Fast basecalling (Guppy) - This step is only for demultiplexing the fast5's for faster processing times
#fast5 --> fastq
#Demultiplexing (fastq --> demultiplexed fastq's in different folders, list of readid's per barcode)
#splitting fast5's (using readid's)
#Methylation sensitive basecalling and mapping to template (Megalodon, uses Guppy)
##Corrects sequence based on template
##Maps sequence onto template
##Calls methylation at A's and C's
#fast5 --> fastq and bam files
#bam files have both sequence and methylation information that we will read out using Pysam in a Jupyter notebook

# DiMeLo_v2_invitro_template_040323_KS.ipynb

#Functions and code for plotting MEATseq m6A probability stats and read clustering
#jup notebook template for analyzing 1x601 Chromatin DiMeLo-seq and sensitivity measurements, chromatinization in extract etc.
#Only 500 reads per barcode, this can be changed but typically 500 reads is sufficient for these analyses.
#Input Megaladon output - mod_mappings.bam (filter by length to keep file size small)
#Template - 1x601.fa
#Functions:
#parse_bam6_fast - takes in bam file, reads through each read entry and outputs prob score at each A and coverage at each position along the given template. Other outputs include smoothed/sliding window averaged prob score, binary output for m6A call at each A position, number of modified A's and number of A's per read for other statistics, readids.
#extract_Ml_Aa_fast - takes in Mm and Ml from Megaladon output for each read. converts to array with just m6A calls with the correct index (removes mC prob scores)
#m6A_seqmap_prob_fast - takes in positions of A's and T's along template and assigns prob scores from Ml of megaladon output

# DiMeLo_v1_invitro/

#This folder contains files for in vitro analyses corresponding to Figures 2 and ext. figure 1,2 of Altemose et al 2022 Nature Methods DiMeLo-seq paper

# DiMeLo_v2_invitro_trials/

#This folder contains Guppy and Megalodon commands to make bam files for in vitro analyses (for DiMeLo-seq performed between 10/2022 and 04/2023)

# DiMeLo_v2_invitro_jupyter/

#This folder contains jupyter notebooks to perform aggregate and individual read analyses of methylation for in vitro analyses (for DiMeLo-seq performed between 10/2022 and 04/2023)

# Template_fasta/

#Template DNA sequences used for in vitro DiMeLo-seq analyses

#pGEM_3z_601.fa - fasta file of template for 1x601 DNA sequence (full plasmid sequence)
#pGEM_3z_1X601_730bp.fa - fasta file of template for 1x601 DNA sequence
#ASP696_18x601.fa - fasta file of template for 18x601 DNA sequence (full plasmid sequence)

######
##
# In situ analyses:
##

#The following folders and files correspond to in vitro analyses

# DiMeLo_allmerge_align_Chm13_ChrXY_041122.sh

#Strategy for basecalling and aligning to Hybrid genome with ChrX and ChrY.
#From DiMeLo-seq paper:
#Basecalling for centromere enriched samples was performed twice both times using Guppy (5.0.7). The first basecalling used the “super accuracy” basecalling model (dna_r9.4.1_450bps_sup.cfg), followed by alignment to the CHM13+HG002X+hg38Y reference genome using Winnowmap (v2.03). These alignments were then filtered for only primary alignments and mapq score greater than 10 using samtools view -F 2308 -q10. A second round of basecalling was then performed again using Guppy (5.0.7) but now with the rerio all-context basecalling model (res_dna_r941_min_modbases-all-context_v001.cfg) with --bam_out and --bam_methylation_threshold 0.0. Modified basecalls were then merged by read id with winnowmap alignments to generate bam files with high confidence alignments combined with modification calls for downstream processing. For CENP-A-directed experiments four independent biological replicates were used, and for controls (IgG-directed, free-floating pA-Hia5, and untreated), two independent biological replicates were used. For all samples the first replicate was sequenced on two separate flow cells and all sequencing runs were merged for the final analysis.

# DiMeLo2_winnowmap_to_merge_cleanbam_bw.sh

#This file contains scripts for running Winnowmap run on Guppy Superaccuracy output fastq, then combining this bam with the modified base calls in bam output by Guppy All context. The combined bam is then cleaned (to remove secondary alignments etc.) and then sorted and indexed for further processing 

#Making merged bam files in barcode*_mod/pass/ to sort before merging with alignment bam_out
#This is code to be run after guppy superacc and guppy allcontext basecalling have been performed.

# DiMeLo_v1_insitu/

#This folder contains files for in situ analyses corresponding to Figure 6 and ext. figure 10 of Altemose et al 2022 Nature Methods DiMeLo-seq paper

# DiMeLo_v2_insitu_trials/

#This folder contains Guppy and Winnowmap commands to make bam files for in situ analyses performed on HG002 cells (for DiMeLo-seq performed between 06/2022 and 08/2022 along with Pragya Sidhwani (Notes for these trials are available with Pragya)). Alpha in the title corresponds to samplers after performing alpha-HORRES
#Briefly, barcodes in these trials correspond to the following:
#1 - CENP-A DiMeLo-seq
#2 - H3K9Me3 DiMeLo-seq
#3 - SETDB DiMeLo-seq
#4 - Taser DiMeLo-seq
#5 - CTCF DiMeLo-seq
#6 - Cohesin DiMeLo-seq
#7 - CENP-C DiMeLo-seq
#8 - CENP-T DiMeLo-seq
#9 - IgG DiMeLo-seq control
#10 - Fiberseq control

# basecalling_and_lmnb1_from_NickAltemose/
#files from Nick Altemose as uploaded along with DiMeLo-seq 2022 Nature Methods paper
#aggregate analyses used for optimization of in situ DiMeLo-seq protocol
#aggregate analyses for Lamin (lmnb1)

# ctcf_and_h3k9me3_DiMeLo_from_AnnieMaslan/
#files from Annie Maslan as uploaded along with DiMeLo-seq 2022 Nature Methods paper
#aggregate and single molecule analyses at CTCF sites
#aggregate and single molecule analyses at H3K9Me3 sites
