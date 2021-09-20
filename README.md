# MetaPro:  Meta-analysis pipeline for Transcriptomics + Genomics
---
The MetaPro Meta Transcriptomics/Genomics Pipeline is a software tool that will perform Transcriptomics and Genomic analysis
on Paired and Single Reads of FASTQ data.

# How to install
---
This package is meant to work in conjunction with Docker/Singularity.  All of the prerequisite tools, including the pipeline code
is delivered via the Docker Hub. https://hub.docker.com/r/parkinsonlab/metapro/
Alternatively, individual parts of the pipeline are avaiable from this Github repository.

Therefore, to use this pipeline, Docker (https://www.docker.com/) or Singularity (https://www.sylabs.io/guides/2.6/user-guide/) is needed.
This also means there is nothing to install (tools and code) besides Docker/Docker CE/Singularity


# How to use
---
This pipeline comes with a config.ini file.  The user is meant to change, configure, and contort the file to point to the location of local files and Databases.
Our config file is written with Python's ConfigParser, using the basic interpretation.  
The docker invocation command should look like:

> docker run python3 /pipeline/MetaPro.py -c <your config> -1 <your forward read.fastq> -2 <your reverse read.fastq> -o <your output directory>

if you are invoking the pipeline using Singularity:
> singularity exec <your singularity image.sif> python3 /pipeline/MetaPro.py -c <your config file> -1 <your forward read.fastq> -2 <your reverse read.fastq> -o <your output directory>


The following is an outline of wwhat each of the sections mean:

# Parameters
---
* Output_Folder
Output_Folder is where the user would indicate to the pipeline where you want the output files to be dumped

* Threads
Threads is the number of threads that the pipeline is allowed to use.  The pipeline is dependent on threads and parallelization to operate efficiently

# Configuration File
---
The configuration file is an essential portion of MetaPro's operation.  It governs the internal settings, and tool paths 
that MetaPro will use in the analysis of the sample data.

# Sequences
---
* Single
Single is for single reads.  Only fill this in if your sequence is a single file

* Pair 1 and Pair 2
Pair 1 is for Forward Reads.  Pair 2 is for Reverse Reads.  Some portions of the code rely on the quality of the Forward Read file as a leading indicator of quality for the data filtering processes.  It is imperative that Pair 1 be used for your Forward Reads only, else you run the risk of a bad analysis. 

## Databases
---
We intentionally did not include any database files in this distribution so as to allow the user more flexibility in how they want to use the pipeline.  Databases also become obselete quickly, and the image size would be enormous.  

However, many of these databases require indexing, so we are curating ready-to-go versions these databases on our website:
https://compsysbio.org/metapro_libs/

We have created a convenient downloader/extractor program for you to use, bundled with the software.
> singularity exec <singularity image> python3 /pipeline/lib_downloader.py <your destination> 


Below is a description of each of the fields we used.  
* database_path
This field isn't a part of the parameters that the pipeline accepts.  It's a shortcut argument that makes filling the path to each database easier.

* Univec_Core
The Univec_Core Database is used in the Vector Contaminents removal stage.  A copy can be found at: https://www.ncbi.nlm.nih.gov/tools/vecscreen/univec/
* Adapter
The Adapter Database is used by Trimmomatic to trim Adapter segments out of the sequence files
A copy can be found inside the Trimmomatic tool installer, located at: http://www.usadellab.org/cms/?page=trimmomatic
This pipeline was built and tested using the TruSeq3-PE-2.fa Adapter Database
* Host
The Host Database is used to filter out Host Contaminents from your seuqence file.  You will need to change this with the CDS database of whichever animal was used in your experiment.
We get our CDS databases from the NCBI, eg: ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human
* Rfam
The Rfam Database is used by Infernal, the rRNA filter.
A copy can be found here: http://rfam.xfam.org/
* DNA_DB
The DNA DB is what we use to annotate the sequence data against.  We use the ChocoPhlAn database.
A copy can be found at: http://huttenhower.sph.harvard.edu/humann2_data/chocophlan/chocophlan.tar.gz
* DNA_DB_Split
The ChocoPhlAn database too large for pBLAT to process.  We split it up and process the chunks simultaneously.  The pipeline will split it and dump the chunks at this location
* Prot_DB
The Prot_DB is the protein db.  We use the non-redundant database from NCBI.  It will need to be indexed by DIAMOND before usage. (see DIAMOND for more details: https://github.com/bbuchfink/diamond)
It can be found here: ftp://ftp.ncbi.nlm.nih.gov/blast/db/
* accession2taxid
This database links each accession to a taxid number.  It's used as part of a custom program in the pipeline.
It can be found at: ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/
* nodes
This file is used in various parts in the pipeline.
It can be found at: ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip
* names
It can be found at the same location as nodes (above)
* Kaiju_db
The Kaiju Database is used by Kaiju for taxonomic annotation.  Kaiju requires that the database be indexed before usage.  
Please see: https://github.com/bioinformatics-centre/kaiju/blob/master/README.md for more information.  Once the indexing is complete, the path to the database needs to be provided in this location
* Centrifuge_db
The Centrifuge Database is used for the Centrifuge tool, which is a part of the enzyme annotation phase.  More information can be found here: https://ccb.jhu.edu/software/centrifuge/manual.shtml
We use the Nucleotide Database, after it has been indexed.  Details can be found at the link to Centrifuge
* SWISS_PROT
SWISS_PROT is the SWISS Prot database, now called UniProt.  We use it in DIAMOND during the enzyme annotation stage.  A copy can be found here: https://www.uniprot.org/downloads
Note that this also needs to be indexed by DIAMOND prior to use.
* SWISS_PROT_map
This is a special file generated by a setup script we have included with the pipeline.  To obtain this file, copy the file "Download_Annotated_EC_mapping.py" from /pipeline/setup_scripts to a 
desired location, and run it.  The script will generate the tab-separated-value map file needed for the pipeline
* PriamDB
The database used by PRIAM.  To obtain this, the user needs to download a distribution of PRIAM.  The PriamDB path is looking for the location where the following files are found:
- PROFILES folder
- annotation_rules.xml
- genome_rules.xml
PRIAM can be downloaded here: http://priam.prabi.fr/
* DetectDB
Detect is an enzyme annotation tool. 

* WEVOTEDB
WeVote is a taxonomy consensus tool that determines taxonomy, given a collection of possible results.

* EC_pathway
This file translates ECs to pathways.  It is needed for the CytoScape data generation

* path_to_superpath
This file is also for the CytoScape data generation

* MetaGeneMark_model
The MetaGeneMark model used by MetaGeneMark to identify genes within contigs.

## Settings
---
These settings are used by MetaPro to orchestrate the concurrency of the tools in action, as well as control some quality-of-life functions within MetaPro
* AdapterRemoval_minlength
This is the minimum length for an adapter

*Target_Rank
Sets the 




## Important note for MetaGeneMark
---
MetaGeneMark is free to use, but requires a license.  If your MetaGeneMark license is not valid, MetaPro will not proceed past the Contig Assembly stage.
To acquire a valid MetaGeneMark license, please visit: http://exon.gatech.edu/GeneMark/license_download.cgi
Select MetaGeneMark, along with the linux 64bit version, and fill in the required fields.  Then download the 64bit key.  Extract the key to obtain the .gm_key file.
The expected location of the license file (.gm_key) is in the home folder of your singularity instance.


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

## Verbose-mode
---
The "keep" and "quiet" settings to verbose_mode will use additional time to compress (keep) or delete (quiet) the interim files produced by the pipeline.  If the performance of a single run is the priority, the "verbose" option should be used to avoid this overhead. 


## Guide to outputs
---
MetaPro produces many outputs for the user to use:
- Gene map
The Gene map shows all of the genes and proteins that have been identified by either BWA, BLAT, or DIAMOND.
The map shows every read that has been identified to those genes/proteins.
The columns in the Gene map are: 
-- Gene ID / Protein accession
-- Length of the gene / protein
-- Number of reads annotated to the gene/protein
-- The read IDs that annotated to the gene/protein

- RPKM table
The RPKM table shows the abundance of each gene, expressed as RPKM of the reads annotated to the gene/protein.
The columns are:
-- The ID of the gene/protein
-- The length of the gene/protein
-- The number of reads associated with each gene/protein
-- The ECs associated with each gene/protein (delineated with a "|" to show multiple ECs per gene/protein)
-- The RPKM of the gene/protein
-- The rest of the columns are taxa, representing at least 1% of the sample (on default), and the RPKM associated by taxa.

- EC heatmap RPKM
This table is not meant to be used externally.  It is for generating the EC heatmap
This table contains no new information from RPKM, and is a transformed table to flatten the data where there are multiple ECs per gene/protein

- Taxa Classifications
This table shows the taxa classification for every putative read that MetaPro processed.
The columns are:
-- Classified/Unclassified.  C is classified, U is unclassified.  
-- The read ID
-- The full taxonomic tree of the read

