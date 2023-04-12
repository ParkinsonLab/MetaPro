# __MetaPro:  Meta-analysis pipeline for Transcriptomics + Genomics__
---
The MetaPro Meta Transcriptomics/Genomics Pipeline is a software tool that will perform Transcriptomics and Genomic analysis
on Paired and Single Reads of FASTQ data.
This readme is for getting started with setting up and using MetaPro.

We have prepared a tutorial demonstrating bioinformatics concepts through the use of MetaPro.  It can be found here:
https://github.com/ParkinsonLab/MetaPro_tutorial

# __How to get the pipeline__
---
This package is meant to work in conjunction with Docker/Singularity.  All of the prerequisite tools, including the pipeline code
is delivered via the Docker Hub. https://hub.docker.com/r/parkinsonlab/metapro/
Alternatively, individual parts of the pipeline are avaiable from this Github repository.

Therefore, to use this pipeline, Docker (https://www.docker.com/) or Singularity (https://www.sylabs.io/guides/2.6/user-guide/) is needed.
This also means there is nothing to install (tools and code) besides Docker/Docker CE/Singularity


# __How to install MetaPro__
---
There is no installation for the pipeline. MetaPro is a script, bundled into a docker image. 
However, in order to use MetaPro, there are a few steps to perform:
## Database files
We recommend using the Databases and libraries we have created for 
We have created a convenient downloader/extractor program for you to use, bundled with the software

    singularity exec <singularity image.sif> python3 /pipeline/lib_downloader.py <your destination> 

This script will download, and unpack the libraries we host to your desired location.  The details of the libraries can be found further down this document, under the configurations section.

## MetaGeneMark license
MetaPro uses MetaGeneMark in its workflow.  To use this software, the user must obtain a free license for the 64-bit version here:

http://exon.gatech.edu/license_download.cgi

The resulting license file (.gm_key) must be placed in the home directory of the docker/singularity environment

    /home/<your username>

# __How to use MetaPro__
---
Since the pipeline is a script, the commands require the user to have a basic understanding of how to invoke a python script. 
This pipeline comes with a config.ini file.  The user is meant to modify the file to point to the location of local files and Databases.
Our config file is written with Python's ConfigParser, using the basic interpretation.  
The docker invocation command should look like:

    docker run python3 /pipeline/MetaPro.py -c <your config.ini> -1 <your forward read.fastq> -2 <your reverse read.fastq> -o <your output directory>
If you are running MetaPro with singletons, the command would change to:
    
    docker run python3 /pipeline/MetaPro.py -c <your config.ini> -s <your singletons.fastq> -o <your output directory>

if you are invoking the pipeline using Singularity:

    singularity exec <your singularity image.sif> python3 /pipeline/MetaPro.py -c <your config file> -1 <your forward read.fastq> -2 <your reverse read.fastq> -o <your output directory>

If you are running MetaPro with singletons, the command would change to:

    singularity exec <your singularity image.sif> python3 /pipeline/MetaPro.py -c <your config file> -s <your     singletons.fastq> -o <your output directory>

MetaPro is built with runtime arguments and options:
    
    --help: call a menu displaying the runtime functions
    -o | --output: sets the output path where all files will be generated
    -c | --config: sets the path to the configuration file
    -1 | --pair1: expects the absolute path for the forward reads of a paired-end dataset
    -2 | --pair2: expects the absolute path for the reverse reads of a paired-end dataset
    -s | --single: expects the absolute path for single-end data (singletons).  Note: If the path to singletons are set, the data is assumed to be single-end, even if -1 and -2 are set. 
    -con | --contig: (for tutorial use only) This argument expects contig segments
    --nhost | --no-host: option to set MetaPro to bypass host filtering
    --tutorial: calls MetaPro in tutorial mode.  This mode lets users run one step at a time. It requires specific keywords to invoke, which are explained below 

Tutorial mode lets users call MetaPro to run only 1 specific stage at a time.
To use this mode, the user must set the paths to the input data and call MetaPro with the tutorial arguemnt, and a keyword.
The input data is entered through the arguments -1, -2, -s, and -con
The keywords are:
- quality: to run quality filtering
- host: to run host filtering
- vector: to remove vectors
- rRNA: to remove rRNA and tRNA, leaving only putative mRNA
- repop: to restore duplicates pulled from the quality step (requires the quality step to be run first)
- contigs: to assemble contigs from reads
- GA: performs gene annotation (with all 3 tools)
- TA: performs taxonomic annotation (requires the GA step to be run first)
- EC: performs enzyme classification (requires the GA step to be run first)
- output: collects all the data from GA, TA, and EC to make tables. (requires GA, TA, and EC to be run first)

# __How to schedule on a job scheduler__
---
MetaPro is a python script.  To schedule MetaPro on a job scheduler, the user should call on singularity to call on python, to call on the MetaPro.py script.
An example would be:

    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --cpus-per-task=40
    #SBATCH --time=24:00:00
    #SBATCH --job-name mpro_501_cqy
    #SBATCH --output=/scratch/j/jparkin/billyc59/kneaddata_run/run_logs/mpro_501_cqy_out.txt
    #SBATCH --error=/scratch/j/jparkin/billyc59/kneaddata_run/run_logs/mpro_501_cqy_err.txt
    image='/home/j/jparkin/billyc59/metapro_develop.sif'
    read1='/scratch/j/jparkin/billyc59/kneaddata_run/NOD501CecQY_r1_new.fastq'
    read2='/scratch/j/jparkin/billyc59/kneaddata_run/NOD501CecQY_r3_new.fastq'
    config='/scratch/j/jparkin/billyc59/kneaddata_run/config_mouse.ini'
    output='/scratch/j/jparkin/billyc59/kneaddata_run/mpro/501_cqy'
    singularity exec -B /home -B /project -B /scratch $image python3 /pipeline/MetaPro.py -c $config -1 $read1 -2 $read2 --verbose_mode leave -o $output

# __Configuration file__
---
MetaPro relies on a configuration file to set program operation directives, paths, and various filter limits.
To get the most out of MetaPro, we recommend becoming familiar with the following:

### Databases
---
The libraries MetaPro uses are critically important for the user to use the pipeline properly.

We intentionally did not include any database files in the Docker/Singularity image due to size constraints.
We also want to highlight the fact that users can use their own databases for MetaPro based on their own needs.

However, many of these databases require indexing, so we are curating ready-to-go versions these databases on our webserver:

    https://compsysbio.org/metapro_libs/
    
There is also a library downloader script embedded in MetaPro, located at:
    
    python3 /pipeline/lib_downloader.py <your destination folder> <optional keyword to download selective elements>
    
This downloader script is built automate the downloading and extraction of MetaPro's databases from our webserver.  However, it also comes with the ability to download and extract individual portions. The optional keywords are:

    all: downloads + extacts everything, default for leaving the option blank is <all>
    detect: detect2's libraries
    ec_pathway: EC-to-pathway map <needed for reporting>
    priam: PRIAM libraries
    rfam: rfam libraries
    wevote: wevote names and nodes
    accession2taxid: accession -> taxid map for TA
    centrifuge: centrifuge libraries
    chocophlan: This option downloads + extracts all 3 copies of our taxa-grouped ChocoPhlAn
    choco_genus: downloads + extracts only the genus-grouped ChocoPhlAn
    choco_family: downloads + extracts only the family-grouped ChocoPhlAn
    choco_order: downloads + extracts only the order-grouped ChocoPhlAn
    kraken2: kraken2 libraries. k2_standard_20230314
    nr: NR from October 10, 2021 
    path_to_superpath: pathway-to-superpathway map
    swissprot: june 10, 2020
    trimmomatic: for version 0.36
    univec: july 05, 2018
    
    
Below is a description of each of the database fields we used.  

    Database_path

This field isn't a part of the parameters that the pipeline accepts.  It's a shortcut argument that makes filling the path to each database easier.

    Univec_Core

The Univec_Core Database is used in the Vector Contaminents removal stage.  A copy can be found at: https://www.ncbi.nlm.nih.gov/tools/vecscreen/univec/

    Adapter

The Adapter Database is used by Trimmomatic to trim Adapter segments out of the sequence files
A copy can be found inside the Trimmomatic tool installer, located at: http://www.usadellab.org/cms/?page=trimmomatic
This pipeline was built and tested using the TruSeq3-PE-2.fa Adapter Database

    Host

The Host Database is used to filter out Host Contaminents from your sequence file.  This stage is agnostic to which database you use.  We currently use coding sequences (CDS). You will need to change this with the CDS database of whichever animal was used in your experiment.
We get our CDS databases from the NCBI, eg: ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human

    Rfam

The Rfam Database is used by Infernal, the rRNA filter.
A copy can be found here: http://rfam.xfam.org/

    DNA_DB

The DNA DB is used on both BWA and BLAT for gene annotation.  We use a taxa-grouped copy of the ChocoPhlAn database found in HuMAnN3.
A copy can be found at: http://huttenhower.sph.harvard.edu/humann2_data/chocophlan/chocophlan.tar.gz
However, downloading a copy from the Huttenhower lab will not contain our modifications

However, we also allow the use of custom databases.  Custom databases will still need to be indexed by BWA.
Note: If a database is larger than 5GB, it can still be used, but it will need to be split.  Each file of the split database must not exceed 5GB.  pBLAT cannot handle a file larger than 5GB. 

    Prot_DB

The Prot_DB is the protein db.  We use the non-redundant database from NCBI.  It will need to be indexed by DIAMOND before usage. (see DIAMOND for more details: https://github.com/bbuchfink/diamond)
It can be found here: ftp://ftp.ncbi.nlm.nih.gov/blast/db/

    accession2taxid

This database links each accession to a taxid number.  It's used as part of a custom program in the pipeline.
It can be found at: ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/

    nodes
    names

These files are used in various parts in the pipeline.  They are NCBI's list of known taxa <nodes> and their associated scientific names <names>
It can be found at: ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip

    kraken2_db

The Kraken2 Database is used by Kraken2 for taxonomic annotation.  The files can be found at: https://benlangmead.github.io/aws-indexes/k2
The creators of this database update it frequently.  We recommend getting a copy from the source.  However, in case the database is not available, we host a copy 
on our webserver.

    Centrifuge_db

The Centrifuge Database is used for the Centrifuge tool, which is a part of the enzyme annotation phase.  More information can be found here: https://ccb.jhu.edu/software/centrifuge/manual.shtml
We use the Nucleotide Database, after it has been indexed.  Details can be found at the link to Centrifuge

    SWISS_PROT

SWISS_PROT is the SWISS Prot database, now called UniProt.  We use it in DIAMOND during the enzyme annotation stage.  A copy can be found here: https://www.uniprot.org/downloads
Note that this also needs to be indexed by DIAMOND prior to use.

    SWISS_PROT_map

This is a file generated by a setup script we have included with the pipeline.  To obtain this file, copy the file "Download_Annotated_EC_mapping.py" from /pipeline/setup_scripts to a 
desired location, and run it.  The script will generate the tab-separated-value map file needed for the pipeline

    PriamDB

The database used by PRIAM.  To obtain this, the user needs to download a distribution of PRIAM.  The PriamDB path is looking for the location where the following files are found:
- PROFILES folder
- annotation_rules.xml
- genome_rules.xml
PRIAM can be downloaded here: http://priam.prabi.fr/

    DetectDB

Detect is an enzyme annotation tool. 

    WEVOTEDB

WeVote is a taxonomy consensus tool that determines taxonomy, given a collection of possible results.

    EC_pathway

This file translates ECs to pathways.  It is needed for the CytoScape data generation

    path_to_superpath

This file is also for the CytoScape data generation

    MetaGeneMark_model: /pipeline_tools/mgm/MetaGeneMark_v1.mod

The MetaGeneMark model used by MetaGeneMark to identify genes within contigs.  This is included in the MetaPro docker/singularity image.  No need to change this path.

    taxid_class_map:
    
    

## Tools
---
MetaPro uses this section of the configuration to control what version of which tool is called.
All of the default tools are installed onto the docker/singularity container.

    Python: python3 
    Java: java -jar (For running jarfiles)
    cdhit_dup: version 4.6.8-2017-1208
    adapterremoval: version 2.1.7 
    vsearch: version 2,7.1
    BWA: version 0.7.17
    SAMTOOLS: verison 1.8
    BLAT: pBLAT version 2.0
    DIAMOND: version 0.9.19 
    BlastP: version 2.7.1+
    Barrnap: version 0.8
    Needle: A tool within emboss (using version 6.6.0)
    Infernal: version 1.1.2-linux
    Kaiju: version 1.6.2
    Centrifuge: latest version
    PRIAM: We host this one ourselves, due to there being intermittent downtimes on the site
    DETECT: Using a modified copy of DETECTv2, built for parallel execution (Detect 2.2.9)
    BLAST_dir: prerequisite for PRIAM 
    WEVOTE: latest version
    Spades: version 3.15.2
    MetaGeneMark: There is only 1 version


If the user wants to change the tool version, they are free to, but the pipeline is not gauranteed to work if there are significant formatting changes to the input, or output

# Settings
---
The 3rd section of MetaPro's configuration file controls the various runtime resource settings.

    AdapterRemoval_minlength: sets the minimum basepair length for AdapterRemoval (default 30)
    Show_unclassified: tells MetaPro's RPKM report to include unclassified reads as a column (default: No)
    bypass_log_name: lets debuggers choose the file that will be the bypass log (default: bypass_log.txt)
    debug_stop_flag: lets debuggers stop MetaPro after a specific stage (default: none)
    num_threads: lets users set the number of threads that MetaPro can use (default: max system limit)
    taxa_exist_cutoff: decides what the percentage cutoff of representative taxa found in the GA-pre-scan.  (default: 0.1)
    DNA_DB_mode: indicates to MetaPro which database is being used (options: custom, default: chocophlan). If set to custom, MetaPro will bypass GA-pre-scan, but will perform a database check for the custom entry.
    filter_stringency: controls what MetaPro will do to resolve corner cases in paired-end read filtering. Options: high (accept reads only when forward and reverse end pass filter. |  low (reject reads only when forward and reverse match to a filter database). default: high
  
  
The mem threshold settings control the amount of remaining RAM that must exist for additional instances of each tool to be launched (available RAM is measured during the initial launch)
    
    BWA_mem_threshold: default 50(%)
    BLAT_mem_threshold: default 50(%)
    DIAMOND_mem_threshold: default 50(%)
    BWA_pp_mem_threshold: default 50(%)
    BLAT_pp_mem_threshold: default 50(%)
    DIAMOND_pp_mem_threshold: default 50(%)
    DETECT_mem_threshold: default 50(%)
    Infernal_mem_threshold: default 50(%)
    Barrnap_mem_threshold: default 50(%)
    TA_mem_threshold: default 50(%)
    
Job limits: controls the amount of concurrent jobs allowed to launch at any one time.  This setting should be changed to fit the user's single-node thread-count

    BWA_job_limit: default (80% of max CPUs on system)    
    BLAT_job_limit: default (80% of max CPUs on system)
    DIAMOND_job_limit: default (80% of max CPUs on system)
    BWA_pp_job_limit: default (80% of max CPUs on system)
    BLAT_pp_job_limit: default (80% of max CPUs on system)
    DIAMOND_pp_job_limit: default (80% of max CPUs on system)
    DETECT_job_limit: default (80% of max CPUs on system)
    Infernal_job_limit: default (80% of max CPUs on system)
    Barrnap_job_limit: default (80% of max CPUs on system)
    TA_job_limit: default (80% of max CPUs on system)
    
The job delay settings control the amount of time the master controller waits until a new process is started.  This is used in conjuction with the memory limit settings to avoid out-of-memory issues, and process-kill issues as the full RAM usage isn't realized until some of the programs get underway.

    BWA_job_delay: default 5 (seconds)
    BLAT_job_delay: default 5 (seconds)
    DIAMOND_job_delay: default 5 (seconds)
    BWA_pp_job_delay: default 5 (seconds)
    BLAT_pp_job_delay: default 5 (seconds)
    DIAMOND_pp_job_delay: default 5 (seconds)
    TA_job_delay: default 5 (seconds)
    
Keep settings: controls whether the interim data is kept in MetaPro. 

    keep_all | default: yes
    keep_quality | default: no
    keep_host | default: no
    keep_vector | default: no
    keep_rRNA | default: no
    keep_repop | default: no
    keep_assemble_contigs | default: no
    keep_GA_BWA | default: no
    keep_GA_BLAT | default: no
    keep_GA_DIAMOND | default: no
    keep_GA_final | default: no
    keep_TA | default: no
    keep_EC | default: no
    keep_outputs | default: no
    
## Labels
---
Labels are how MetaPro names each stage.  In cases where users want to test different settings in GA/TA/EC, or to try a battery of trials using the same data, we now let users name their stages in customized ways:

    quality_filter | default: quality_filter
    host_filter | default: host_filter
    vector_filter | default: vector_filter
    rRNA_filter | default: rRNA_filter
    repop | default: duplicate_repopulation
    assemble_contigs | default: assemble_contigs
    GA_pre_scan | default: GA_pre_scan
    GA_split | default: GA_split
    GA_BWA | default: GA_BWA
    GA_BLAT | default: GA_BLAT
    GA_DIAMOND | default: GA_DMD
    GA_final_merge | default: GA_final_merge
    ta | default: taxonomic_annotation
    ec | default: enzyme_annotation
    outputs | default: outputw
    


## Important note for MetaGeneMark
---
MetaGeneMark is free to use, but requires a license.  If your MetaGeneMark license is not valid, MetaPro will not proceed past the Contig Assembly stage.
To acquire a valid MetaGeneMark license, please visit: http://exon.gatech.edu/GeneMark/license_download.cgi
Select MetaGeneMark, along with the linux 64bit version, and fill in the required fields.  Then download the 64bit key.  Extract the key to obtain the .gm_key file.
The expected location of the license file (.gm_key) is in the home folder of your singularity instance.
Place your license in:
> /home/<your user name>

# Important Features
---
This pipeline was built with a few key features that make it easier to use than other pipelines.
## Verbose-mode
This pipeline will produce many interim files along the way.  In an effort to maintain transparency, and allow for debugging and inspection of stages, the interim files have been saved, and compressed at each stage of the pipeline.  If the user wishes to avoid compression (saving roughly 1/100th of the runtime), the compress flag argument should be set to "verbose".  Otherwise, the default is to compress the interim files.

## Auto-resume
The pipeline is capable of skipping any and all stages it has run before.  If the user wishes to run a subset of stages, the user has to remove the stage folders, and any compressed files of the accompanying stage, and the pipeline will re-run the missing stages.  Note:  The removed stages do not have to be contiguous, but it is recommended that they are, to ensure the accuracy of the pipeline results.  Eg: a pipeline runs A -> B -> C -> D -> E.  If A, and D are removed, The pipeline will re-run A and D, leaving B, C, and E alone.  However, if D depends on C, then the changes from A -> B -> C will not propagate downward.  This feature also has the benefit of being able to resume running if the pipeline run was operating on a job-scheduler-controlled system, and inadequate time was allocated to the job.  The pipeline will simply pick up where it left off until the job is complete. 

## Auto-death
In an effort to save computational resources, the pipeline will shut itself down if the dependencies of a stage are not met.  This is to ensure that if a pipeline run is part of a batch-processing scheduler, the pipeline will not continue to waste resources on an erroneous job.  

# Increasing performance
---
## Operating mode
MetaPro runs in a Singularity instance.  As of writing (Sept 28, 2018), Singularity does not support multi-machine parallelism.  This pipeline does not utilize MPI, but instead strives to use all the cores made available by the singularity machine through the Python Multiprocessing module.  To increase the performance of the pipeline, more cores should be given to the host machine, and increasing the number of cores the pipeline is allowed to use.


## Guide to outputs
---
MetaPro produces many outputs for the user to use:
### read count
The read count table (read_count.tsv) is the full accounting of the sequence reads as they progress through MetaPro
The columns are:

    Total reads: The number of unique raw sequence reads in the sample
    High quality: The number of unique sequence reads after the low-quality filtering step
    % high quality: High quality, as a percentage of total reads
    host reads: the number of host sequence reads found by MetaPro
    % host reads: host reads, as a percentage of total reads
    vector reads: the number of vector reads found by MetaPro
    % vector reads: vector reads, as a percentage of total reads
    rRNA + tRNA reads: the number of rRNA and tRNA reads removed by MetaPro
    % rRNA + tRNA reads: rRNA + tRNA reads, as a percentage of total reads
    Putative reads: the number of mRNA reads that will be annotated by MetaPro
    % Putative reads: putative reads, as a percentage of total reads
    annotated mRNA reads: the number of putative reads that were assigned to a gene or protein
    % annotated reads: annotated mRNA reads, as a percentage of putative reads
    Unique transcripts: the number of unique genes and proteins found by MetaPro
    High-Quality enzymes: The number of enzymes found by MetaPro (using the high-quality EC filter settings)
    Low-Quality enzymes: the number of enzymes found by MetaPro (usiing the low-quality EC filter settings)

### Contig stats
---
This file displays the N50 and L50 read lengths of the contigs assembled by rnaSPADes.  

### Gene map
---
The Gene map shows all of the genes and proteins that have been identified by either BWA, BLAT, or DIAMOND.
The map shows every read that has been identified to those genes/proteins.

The columns in the Gene map are: 

    Gene ID / Protein accession
    Length of the gene / protein
    Number of reads annotated to the gene/protein
    The read IDs that annotated to the gene/protein

### RPKM table
---
The RPKM table shows the abundance of each gene, expressed as RPKM of the reads annotated to the gene/protein.
The columns are:

    The ID of the gene/protein
    The length of the gene/protein
    The number of reads associated with each gene/protein
    The ECs associated with each gene/protein (delineated with a "|" to show multiple ECs per gene/protein)
    The RPKM of the gene/protein
    The rest of the columns are taxa, representing at least 1% of the sample (on default), and the RPKM associated by taxa.

### Enzyme superpathway heatmap
---
This table is not meant to be used externally.  It is for generating the EC heatmap
This table contains no new information from RPKM, and is a transformed table to flatten the data where there are multiple ECs per gene/protein
Each row corresponds to each enzyme superpathway, while each column will be one of the major organisms (at least 1%) in the sample

### Taxa Classifications
---
This table shows the taxa classification for every putative read that MetaPro processed.
The columns are:
    
    Classified/Unclassified.  C is classified, U is unclassified.  
    The read ID
    The full taxonomic tree of the read

### Taxa Summary
---
This table is a tally of all unique taxa trees found in the sample
The columns are:

    taxa: the full taxonomic tree 
    count: the number of times a sequence read annotated to this specific taxa tree

