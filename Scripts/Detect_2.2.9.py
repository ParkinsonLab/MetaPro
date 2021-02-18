#!/usr/bin/env python

import os
import sqlite3
import subprocess
from argparse import ArgumentParser
from collections import defaultdict, OrderedDict
from operator import itemgetter
from string import whitespace
from sys import stdout
from Bio import SeqIO
from datetime import datetime as dt
import sys
import queue as q
import pandas as pd
import multiprocessing as mp
import psutil as psu
import time
import resource as rs

def mem_checker(threshold):
    #threshold is a percentage
    mem = psu.virtual_memory()
    available_mem = mem.available
    total_mem = mem.total
    
    available_pct = 100 * available_mem / total_mem
    
    if(float(available_pct) <= float(threshold)):
        return False
    else:
        return True
        
def cpu_checker():
    cpu_use = psu.cpu_percent(interval = 0.0001, percpu = True)
    for item in cpu_use:
        if(item == 0):
            return True
        
    return False

def make_folder(folder_path):
    if not (os.path.exists(folder_path)):
        os.makedirs(folder_path)

def get_ec_to_cutoff(mapping_file, beta):
    """Return the mapping of EC to cutoff from the file with the mapping."""
    ec_to_cutoff = {}
    #open_file = open(mapping_file)
    with open(mapping_file, "r") as open_file:
        for i, line in enumerate(open_file):
            line = line.strip()
            if line == "":
                continue
            # which column to use?
            if i == 0:
                headers = line.split("\t")
                for col_of_interest, header in enumerate(headers):
                    if header.find("beta=" + str(beta)) != -1:
                        break
                continue
            split = line.split()
            ec, cutoff = split[0], float(split[col_of_interest])
            ec_to_cutoff[ec] = cutoff
    #open_file.close()
    return ec_to_cutoff

class PairwiseAlignment:
    """Structure to store the numbers associated with a Needleman Wunsch pairwise alignment.
    Keeps track of two sequence ids (represented by Sequence object) and associated alignment score (float)"""
    def __init__(self, query, hit, score):  
        self.query = query
        self.hit = hit
        self.score = score

class Sequence: 
    """Represents a FASTA sequence string with its header"""
    def __init__(self, header, data):
        self.header=header
        self.data=data
    
    """Return an indetifier from the fasta sequence
    First non-whitespace string, excluding the first character (>)"""
    def name (self):
        return self.header.split()[0][1:]

    """Return the complete FASTA sequence with both header and sequence data
    """
    def fasta (self):
        return "\n".join([self.header,self.data])

class Identification:
    """Represents a functional identification of a sequence toward an EC number
    Hypotheses is a possibly redundant list list of Hypothesis objects.
    The probability of a hypothesis being correct is calculated using the Bayes theorem.
    Please address Hung et al. (2010) for more details.
        Hung S, Wasmuth J, Sanford C & Parkinson J.
        DETECT - A Density Estimation Tool for Enzyme ClassificaTion and its application to Plasmodium falciparum.
        Bioinformatics. 2010 26:1690-1698
    This probability represents a singular alignment match event.
    Predictions is a non-redundant set of ec numbers associated with cumulative probabilities.
    The probability of a prediction is a cumulative probability of all hypotheses with the same EC number.
    """
    def __init__(self, query_id):
        self.query_id = query_id
        self.hypotheses = list()
        self.predictions = defaultdict(self.__one)
        self.prediction_count = defaultdict(int)

    """A callable function to initiate new values in a defaultdict to float 1
    """
    def __one(self):
        return 1.0

class Hypothesis:
    """Represents a single alignemnt result with an associated probability, as calculted using the Bayes theorem.
    EC is retrieved from the swiss-to-EC mapping database.
    """
    def __init__(self, swissprot_id,score):
        self.swissprot_id = swissprot_id
        self.score = score
        self.ec = "unknown"
        self.probability= 0.0

top_predictions_count = 5
probability_cutoff = 0.2
verbose=False
zero_density = 1e-10
"""Small number that is used as zero"""

def run_pair_alignment (seq, blast_db, num_threads, e_value_min, bitscore_cutoff, uniprot_df, blastp, needle, dump_dir):
    
    
    """Core alignment routine.
    1) Takes a single sequence, acquires multiple BLASTp alignemnts to the swissprot enzyme database.
    2) Canonical sequences of the results from (1) are retrieved from dictionary of ids to swissprot records derived
    from a swissprot fasta
    3) Original query is globally aligned versus sequences from (2)
    """
    
    #First pass cutoff with BLAST alignments
    if verbose: print( "[DETECT]: Running BLASTp for {} ...".format(seq.name()))

    invalid_chars = ["?","<",">","\\",":","*","|"]
    valid_seq_name = seq.name()
    for char in invalid_chars:
        valid_seq_name = valid_seq_name.replace(char, "_")
    try:
        p = subprocess.Popen((blastp, "-query", "-", 
                        "-out", "-",
                        "-db", blast_db,
                        "-outfmt", "6 sseqid bitscore",
                        "-max_target_seqs", "100000",
                        "-num_threads",str(num_threads),
                        "-evalue", str(e_value_min)),
                    stdin=subprocess.PIPE,  
                    stdout=subprocess.PIPE,
                    encoding='utf8')
        stdout,stderr = p.communicate(seq.data)
    except Exception as e:
        print(dt.today(), "BLASTp FAILED", e)
        sys.exit()
    
    blast_hits_path = dump_dir + "blast_hits_" + valid_seq_name
    #print("creating:", blast_hits_path)
    #blast_hits = open(blast_hits_path,"w")
    blast_hit_list = list() 
    for line in stdout.split("\n"):
        if not line in whitespace:
            swissprot_id,bitscore = line.split("\t")
            if float(bitscore) > bitscore_cutoff:
                blast_hit_list.append(">" + str(swissprot_id))
                
    #don't continue if there's nothing from BLAST
    if(len(blast_hit_list) == 0):
        #print(dt.today(), "No BLAST hits found. ending process")
        return list()
    else:
        blast_df = uniprot_df[uniprot_df.index.isin(blast_hit_list)]
        #print("BLAST LIST:", blast_hit_list)
        #print("blast_df:", blast_df)
        if(blast_df.empty):
            print(blast_hits_path, "BLAST_DF is empty... ERROR")
            return list()
        blast_df.to_csv(blast_hits_path, mode="w+", sep='\n', header=False, quoting = 3)
        

    #Run Needleman-Wunsch alignment on the results of the BLAST search
    
    gapextend_value = "0.5"
    if(needle.endswith("stretcher")):
        gapextend_value = "2"
        
    try:
        p = subprocess.Popen((needle, "-filter",
                        "-bsequence", blast_hits_path, #"blast_hits_" + valid_seq_name,
                        "-gapopen", "10",
                        "-gapextend", gapextend_value,
                        "-sprotein", "Y",
                        "-aformat_outfile", "score"),
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    encoding='utf8')
        stdout,stderr = p.communicate(seq.fasta())
    except Exception as e:
        print(dt.today(), "NEEDLE FAILED:", e)
        #print("deleting:", blast_hits_path)
        os.remove(blast_hits_path)
        sys.exit()
    #print("deleting:", blast_hits_path)
    os.remove(blast_hits_path)
    return parse_needle_results(stdout)

def import_fasta(file_name_in):
    #this imports the fasta into a pandas dataframe
    fasta_df = pd.read_csv(file_name_in, error_bad_lines=False, header=None, sep="\n")  # import the fasta
    fasta_df.columns = ["row"]
    #There's apparently a possibility for NaNs to be introduced in the raw fasta.  We have to strip it before we process (from DIAMOND proteins.faa)
    fasta_df.dropna(inplace=True)
    new_df = pd.DataFrame(fasta_df.loc[fasta_df.row.str.contains('>')])  # grab all the IDs
    new_df.columns = ["names"]
    new_data_df = fasta_df.loc[~fasta_df.row.str.contains('>')]  # grab the data
    new_data_df.columns = ["data"]
    fasta_df = new_df.join(new_data_df, how='outer')  # join them into a 2-col DF
    fasta_df["names"] = fasta_df.fillna(method='ffill')  # fill in the blank spaces in the name section
    fasta_df["names"] = fasta_df["names"].apply(lambda x: x.split(" ")[0])
    fasta_df.dropna(inplace=True)  # remove all rows with no sequences
    fasta_df.index = fasta_df.groupby('names').cumcount()  # index it for transform
    temp_columns = fasta_df.index  # save the index names for later
    fasta_df = fasta_df.pivot(values='data', columns='names')  # pivot
    fasta_df = fasta_df.T  # transpose
    fasta_df["sequence"] = fasta_df[fasta_df.columns[:]].apply(lambda x: "".join(x.dropna()), axis=1)  # consolidate all cols into a single sequence
    fasta_df.drop(temp_columns, axis=1, inplace=True)
    #fasta_df["names"] = fasta_df.index
    return fasta_df
    
    
"""Split a fasta file into separate sequences, 
    return a list of Sequence objects.
    See class definition below"""

def split_fasta(input_file):
    #resultant array of peptide sequences
    sequences=list()

    #temporary array of lines for each sequence as it is being read
    seq_buffer=list()
    
    header = ""
    for line in open(input_file):
            
        #if we find a header line and there are already lines in sequence buffer
        if ">" in line and seq_buffer:
            #flush lines from the buffer to the sequences array
            sequences.append(Sequence(header,"".join(seq_buffer)))
            seq_buffer=list()
            header = line.strip()

        #skip the ID line
        elif ">" in line and not seq_buffer:
            header = line.strip()
        
        #add next line to sequence buffer
        else:
            seq_buffer.append(line.strip())

    #dont forget the last sequence
    sequences.append(Sequence(header,"".join(seq_buffer)))

    return sequences

"""Parse tab-delimited BLAST results,
    return only hit IDs.
    Output generated with blastp -outfmt 6 (NCBI BLAST 2.2.26)
    BLAST docs http://www.ncbi.nlm.nih.gov/books/NBK1763 
        output format arguments: section 4.2.26"""
def get_blast_hit_identifiers (input_file):
    
    #resultant array of hit identifiers
    hit_IDs=list()
    
    #results are stored in-string as <query_id>\t<hit_id>\t...

    for line in open(input_file):
        hit_ID = line.strip().split("\t")[1]
        hit_IDs.append(hit_ID)

    return hit_IDs

"""Parse EMBOSS-needle output, 
    return list of structs.
    Output generated with needle -auto -aformat_outfile score (EMBOSS suite 6.4.0.0)
    Details of SCORE format: http://emboss.sf.net/docs/themes/AlignFormats.html
    Needle docs at http://emboss.sourceforge.net/apps/cvs/emboss/apps/needle.html"""
def parse_needle_results (needle_results):
    results=list()

    #results are stored in-string as <query> <hit> <alignment_length> (<score>)
    for line in needle_results.split("\n"):
        #ignore comment lines
        if not "#" in line and not line in whitespace:
            fields = line.strip().split()
            query = fields[0]
            hit = fields[1]
            score = float(fields[3][1:-1])
            h = Hypothesis(hit,score)

            results.append(h)

    return results

def calculate_probability (hypothesis, db_connection):
    score = hypothesis.score
        
    cursor = db_connection.cursor()

    #Fetch the data from tables
    
    #Fetch the EC mapped to the swissprot ID
    cursor.execute("SELECT ec FROM swissprot_ec_map WHERE swissprot_id = '{}'".format(hypothesis.swissprot_id))

    #sqlite3.fetchone() returns a tuple. Since only one value <ec> was requested, this is a one-member tuple. Still,
    # it is important to subset [0]
    mapping = cursor.fetchone()
    if mapping:
        ec = mapping[0]
        hypothesis.ec = ec

        #Get Prior probabilities for that EC
        cursor.execute("SELECT probability FROM prior_probabilities WHERE ec = '{}'".format(ec))
        prior = cursor.fetchone()[0]

        #Get positive density for the given score and EC
        cursor.execute("SELECT density FROM positive_densities WHERE ec = '{}' "
                       "AND score < {} ORDER BY score DESC LIMIT 1 ".format(ec,score))
        previous_point = cursor.fetchone()

        cursor.execute("SELECT density FROM positive_densities WHERE ec = '{}' "
                       "AND score > {} ORDER BY score ASC LIMIT 1".format(ec,score))
        next_point = cursor.fetchone()
        
        if previous_point and next_point:
            positive = (previous_point[0] + next_point[0])/2
        else:
            positive = 0

        #Get negative density for the given score and EC
        cursor.execute("SELECT density FROM negative_densities WHERE ec = '{}' "
                       "AND score < {} ORDER BY score DESC LIMIT 1".format(ec,score))

        previous_point = cursor.fetchone()
    
        cursor.execute("SELECT density FROM negative_densities WHERE ec = '{}' "
                       "AND score > {} ORDER BY score ASC LIMIT 1".format(ec,score))

        next_point = cursor.fetchone()
        
        if previous_point and next_point:
            negative = (previous_point[0] + next_point[0])/2
        else:
            negative = 0

        if positive == 0 and negative == 0:
            probability = zero_density
        else:
            positive_hit = prior * positive
            probability = positive_hit / (positive_hit + ((1.0-prior) * negative ))

        hypothesis.probability = probability
    else:
        probability = 0

    return probability
    
def main_loop(seq, blast_db,num_threads, e_value, bit_score, uniprot_df, blastp, needle, dump_dir, top_pred_flag, fbeta_flag, script_path, out_dict):
    
    identification = Identification(seq.name())
    identification.hypotheses = run_pair_alignment(seq, blast_db,num_threads, e_value, bit_score, uniprot_df, blastp, needle, dump_dir)
    if not identification.hypotheses:
        if verbose: 
            print( "[DETECT]: No BLASTp hits for {}".format(seq.name()))
        
        return None
        
    
    connection = sqlite3.connect(script_path + "/data/detect.db")
    
    if verbose: 
        print( "[DETECT]: Running density estimations for {} ...".format(seq.name()))
    for hypothesis in identification.hypotheses:
        probability = calculate_probability(hypothesis, connection)
        if not (hypothesis.ec == "unknown" or hypothesis.ec == "2.7.11.1" or hypothesis.ec == "2.7.7.6" or hypothesis.ec == "2.7.13.3"):
            identification.predictions[hypothesis.ec] *= (1.0-probability)  
            identification.prediction_count[hypothesis.ec] += 1
    
    low_density = []
    for ec,probability in identification.predictions.items():
        cumulative = 1.0 - probability
        if (cumulative > zero_density):
            identification.predictions[ec] = cumulative
        else:
            low_density.append(ec)
    for ec in low_density:
        del identification.predictions[ec]
    
    if (top_pred_flag or fbeta_flag):
        #sort
        identification.predictions = OrderedDict(sorted(identification.predictions.items(), key=itemgetter(1), reverse=True))
    
    if (top_pred_flag):
        top_predictions = list()

    if (fbeta_flag):
        fbeta_predictions = list()
    id_entry_list = []
    
    for ec in identification.predictions:
        identification_entry = "{seq_id}\t{pred_ec}\t{prob:.3e}\t{pos_hits}\t{neg_hits}\n".format(
                    seq_id=identification.query_id,
                    pred_ec=ec,
                    prob=identification.predictions[ec],
                    pos_hits=identification.prediction_count[ec],
                    neg_hits= len(identification.hypotheses)-identification.prediction_count[ec])
        #write them all to this one place
        id_entry_list.append(identification_entry)
        
        #These are all just alternate write points.
        if top_pred_flag and identification.predictions[ec] > probability_cutoff and len(top_predictions) < top_predictions_count:
            top_predictions.append(identification_entry)

        if fbeta_flag and identification.predictions[ec] > ec_to_cutoff[ec]:
            fbeta_predictions.append(identification_entry)

    final_pack = dict()
    final_pack["id_entry"] = id_entry_list
    if(fbeta_flag):
        final_pack["fbeta"] = fbeta_predictions
    if(top_pred_flag):
        final_pack["top_pred"] = top_predictions
        
    if not(id_entry_list is None):
        out_dict[seq.name()] = final_pack
    
    connection.close()
        # ignore these writes
        #output.write(identification_entry)

        #if (args.top_predictions_file):
        #    for entry in top_predictions:
        #        top_predictions_file.write(entry)

        #if (args.fbeta_file):
        #    for entry in fbeta_predictions:
        #        fbeta_file.write(entry)
def merge_dicts(x, y):
    z = x.copy()
    z.update(y)
    return z

if __name__=="__main__":
    parser = ArgumentParser(description="DETECT - Density Estimation Tool for Enzyme ClassificaTion. "
                                        "Version 2.0. May 2016")
    
    parser.add_argument("target_file",type=str,
                        help="Path to the file containing the target FASTA sequence(s)")
    parser.add_argument("--output_file",type=str,
                        help="Path of the file to contain the output of the predictions")
    parser.add_argument("--verbose",
                        help="print verbose output",action="store_true")
    parser.add_argument("--num_threads",type=int,
                        help="Number of threads used by BLASTp")
    parser.add_argument("--bit_score",type=float,
                        help="The cutoff for BLASTp alignment bitscore")
    parser.add_argument("--e_value",type=float,
                        help="The cutoff for BLASTp alignment E-value")
    parser.add_argument("--top_predictions_file",type=str,
                        help="Path to the file that enumerates predictions with probability over 0.2")
    parser.add_argument("--fbeta_file",type=str,
                        help="Path to the file that enumerates predictions that pass EC-specific cutoffs")
    parser.add_argument("--beta",type=float,choices=[1.0, 0.5, 2.0], default=1.0,
                        help="Value of beta in Fbeta: 1 (default), 0.5 or 2. Fbeta is maximized along EC-specific "
                             "precision-recall curves to derive EC-specific score cutoffs")
    parser.add_argument("--db",type=str,
                        help="Location of the Detect databases")
    parser.add_argument("--blastp",type=str,
                        help="Path for the blastp binary")
    parser.add_argument("--needle",type=str,
                        help="Path for the Needleman-Wunsch search binary")
    parser.add_argument("--dump_dir", type = str, help = "Path for interim files")                    
    
    parser.add_argument("--n_count", type = int, help = "number of concurrent processes")    
    
    parser.add_argument("--mem_limit", type = int, help = "percentage of memory to maintain during the run")
    
    parser.add_argument("--job_delay", type = int, help = "time (seconds) to wait in between parallel jobs")
    
    parser.add_argument("--cpu_interval_time", type = float, help = "time(seconds) to poll the cpu interval time")
    
    args = parser.parse_args()
    script_path = args.db if args.db else os.path.dirname(os.path.realpath(__file__))

    verbose = args.verbose
    num_threads = args.num_threads if args.num_threads else 1
    bit_score = args.bit_score if args.bit_score else 50
    e_value = args.e_value if args.e_value else 1
    blastp = args.blastp if args.blastp else "blastp"
    needle = args.needle if args.needle else "needle"
    DETECT_mem_limit = int(args.mem_limit) if args.mem_limit else 50
    DETECT_job_delay = int(args.job_delay) if args.job_delay else 5
    DETECT_cpu_interval = float(args.cpu_interval_time) if args.cpu_interval_time else 0.001
    
    print(dt.today(), "mem limit used:", DETECT_mem_limit, "%")
    
    
    dump_dir = args.dump_dir
    if not(dump_dir.endswith("/")):
        dump_dir += "/"
    make_folder(dump_dir) #make if doesn't exist
    if(args.n_count):
        n_count = int(args.n_count)
        
    else:
        n_count = 80
    print(dt.today(), "concurrency used:", n_count)    
    
        
    sequences = split_fasta(args.target_file)
    if verbose: print( "Found {} sequences in file.".format(len(sequences)))
    blast_db = script_path+"/data/uniprot_sprot.fsa"
    
    final_predictions = list()

    #connection = sqlite3.connect(script_path + "/data/detect.db")
    if (args.output_file):
        output = open(args.output_file, "w")
    else:
        output = stdout
    
    header = "ID\tEC\tprobability\tpositive_hits\tnegative_hits\n"
    output.write(header)
    
    fbeta_flag = False
    top_pred_flag = False
    
    if (args.top_predictions_file):
        top_predictions_file = open(args.top_predictions_file, "w")
        top_predictions_file.write(header)
        top_pred_flag = True
        
    if (args.fbeta_file):
        fbeta_file = open(args.fbeta_file, "w")
        fbeta_file.write(header)
        fbeta_flag = True
        mapping_file = script_path + "/ec_to_cutoff.mappings"
        ec_to_cutoff = get_ec_to_cutoff(mapping_file, args.beta)

    #ids_to_recs = SeqIO.index(script_path + "/data/uniprot_sprot.fsa", "fasta")
    print(dt.today(), "importing uniprot db")
    uniprot_df = import_fasta(script_path + "/data/uniprot_sprot.fsa")
    print(dt.today(), "finished importing uniprot db")
    #print(uniprot_df)
    print("number of cores:", mp.cpu_count())
    process_counter = 0
    wait_counter = 0
    process_list = []
    
    #collects all the results from the mp threads
    #note:  a regular dict, and Manager-dict is needed because of the limits we place on the concurrent processes.
    #SeqIO has issues with too many threads. In the future, we should get rid of this.
    #but for now, the mechanism is: manager-dicts only survive for as long as there's active processes.
    #There exists points in time where we have no processes running, and so in order for the data to persist, we
    #must 
    final_dict = mp.Manager().dict() #this work because there's no duplicate in DETECT.  
    outer_dict = dict()
    number_of_seqs = len(sequences)
    true_job_count = 0
    max_open_file_limit = rs.getrlimit(rs.RLIMIT_NOFILE)[0]
    max_open_file_limit = int(max_open_file_limit *0.85) #don't go crazy
    #limiter bypass
    #if n_count == 0:
    #    n_count = number_of_seqs
    #if n_count > mp.cpu_count():
    #    print("n_count is NERF'd:", mp.cpu_count())
    #    n_count = mp.cpu_count()
    print("Number of concurrent processes:", n_count)    
    #========================================================================
    # launch the jobs
    for i,seq in enumerate(sequences):
        
    #    if verbose: 
    #        print( "[DETECT]: Analyzing {} ({}/{}) ...".format(seq.name(), i + 1, len(sequences)))
        job_submitted = False
        while (not job_submitted):
            #if(len(process_list) < n_count):
            if(cpu_checker()):
                if(mem_checker(DETECT_mem_limit)):
                    identification = Identification(seq.name())
                    process = mp.Process(
                        target = main_loop,
                        args = (seq, blast_db,num_threads, e_value, bit_score, uniprot_df, blastp, needle, dump_dir, top_pred_flag, fbeta_flag, script_path, final_dict)
                    )
                    process.start()
                    process_list.append(process)
                    process_counter += 1
                    true_job_count += 1
                    job_submitted = True
                    if(true_job_count % 100 == 0):
                        print(dt.today(), true_job_count, "jobs launched!")
                    #print(dt.today(), process_counter, " process launch!")
                    
                    if(process_counter > max_open_file_limit):
                        print(dt.today(), "limit of open jobs reached:" + str(max_open_file_limit) +  "|closing before launching more.")
                        for item in process_list:
                            item.join()
                        process_list[:] = []
                        process_counter = 0
                        
                else:
                    #print(dt.today(), "mem limit reached")
                    time.sleep(0.001)#DETECT_job_delay)
                
            else:
                #print(dt.today(), "process launch paused due to cpu limit")
                time.sleep(0.001)
                
            
        #outer_dict = merge_dicts(outer_dict, final_dict) #merge the 2
        #final_dict = mp.Manager().dict()    #empty it. start a new collection
            
    #======================================================================
    # Final sync
    print(dt.today(), "Final batch of jobs")
    for item in process_list:
        item.join()
    process_list[:] = []
    process_counter = 0
    outer_dict = dict(final_dict)
    #outer_dict = merge_dicts(outer_dict, final_dict) #merge the 2

                        
    
    for item in outer_dict:
        if(fbeta_flag):
            if not(outer_dict[item]['fbeta'] is None):
                for inner_item in outer_dict[item]['fbeta']: 
                    fbeta_file.write(inner_item)
        if(top_pred_flag):
            if not(outer_dict[item]['top_pred'] is None):
                for inner_item in outer_dict[item]['top_pred']:
                    top_predictions_file.write(inner_item)
        for inner_item in outer_dict[item]['id_entry']:
            output.write(inner_item)
        #show_counter += 1
        #if(show_counter > 10):
        #    break
        #print(item, outer_dict[item])
    #    identification.hypotheses = run_pair_alignment(seq, blast_db,num_threads, e_value, bit_score, ids_to_recs, blastp, needle, dump_dir)
        
    
    if (args.top_predictions_file):
        top_predictions_file.close()

    if (args.fbeta_file):
        fbeta_file.close()

    output.close()
    #connection.close()
