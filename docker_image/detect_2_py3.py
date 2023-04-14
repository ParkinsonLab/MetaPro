#!/usr/bin/env python

import subprocess
import os
import sqlite3  
from string import whitespace
from sys import stdout
from argparse import ArgumentParser
from operator import itemgetter
from collections import defaultdict, OrderedDict
from Bio import SeqIO

def get_ec_to_cutoff(mapping_file, beta):
    """Return the mapping of EC to cutoff from the file with the mapping."""
    ec_to_cutoff = {}
    open_file = open(mapping_file)
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
    open_file.close()
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

def run_pair_alignment (seq, blast_db, num_threads, e_value_min, bitscore_cutoff, ids_to_recs):
    """Core alignment routine.
    1) Takes a single sequence, acquires multiple BLASTp alignemnts to the swissprot enzyme database.
    2) Canonical sequences of the results from (1) are retrieved from dictionary of ids to swissprot records derived
    from a swissprot fasta
    3) Original query is globally aligned versus sequences from (2)
    """
    
    #First pass cutoff with BLAST alignments
    if verbose: print( "[DETECT]: Running BLASTp for {} ...".format(seq.name()))
    p = subprocess.Popen(("blastp", "-query", "-", 
                    "-out", "-",
                    "-db", blast_db,
                    "-outfmt", "6 sseqid bitscore",
                    "-max_target_seqs", "100000",
                    "-num_threads",str(num_threads),
                    "-evalue", str(e_value_min)),
                stdin=subprocess.PIPE,  
                stdout=subprocess.PIPE)
    stdout,stderr = p.communicate(seq.data)
    
    with open("blast_hits","w") as blast_hits:
        blast_hit_list = list() 
        for line in stdout.split("\n"):
            if not line in whitespace:
                swissprot_id,bitscore = line.split("\t")
                if float(bitscore) > bitscore_cutoff:
                    blast_hit_list.append(swissprot_id)
        blast_hit_ids = "\n".join(frozenset(blast_hit_list))
            
        if verbose: print( "[DETECT]: Found {} hits for {} ...".format(len(blast_hit_ids),seq.name()))
        
        #stop if nothing found
        if len(blast_hit_ids) == 0:
            return list()

        SeqIO.write((ids_to_recs[hid] for hid in blast_hit_ids.split("\n")), blast_hits, "fasta")

    if verbose: print( "[DETECT]: Running Needleman-Wunch alignments for {} ...".format(seq.name()))

    #Run Needleman-Wunsch alignment on the results of the BLAST search
    p = subprocess.Popen(("needle", "-filter",
                    "-bsequence", "blast_hits",
                    "-gapopen", "10",
                    "-gapextend", "0.5",
                    "-sprotein", "Y",
                    "-aformat_outfile", "score"),
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE)
        
    stdout,stderr = p.communicate(seq.fasta())

    return parse_needle_results(stdout)

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
    
    args = parser.parse_args()
    script_path = os.path.dirname(os.path.realpath(__file__))

    verbose = args.verbose
    num_threads = args.num_threads if args.num_threads else 1
    bit_score = args.bit_score if args.bit_score else 50
    e_value = args.e_value if args.e_value else 1
    
    sequences = split_fasta(args.target_file)
    if verbose: print( "Found {} sequences in file.".format(len(sequences)))
    blast_db = script_path+"/data/uniprot_sprot.fsa"
    
    final_predictions = list()

    connection = sqlite3.connect(script_path + "/data/detect.db")
    if (args.output_file):
        output = open(args.output_file, "w")
    else:
        output = stdout
    
    header = "ID\tEC\tprobability\tpositive_hits\tnegative_hits\n"
    output.write(header)

    if (args.top_predictions_file):
        top_predictions_file = open(args.top_predictions_file, "w")
        top_predictions_file.write(header)

    if (args.fbeta_file):
        fbeta_file = open(args.fbeta_file, "w")
        fbeta_file.write(header)
        mapping_file = script_path + "/ec_to_cutoff.mappings"
        ec_to_cutoff = get_ec_to_cutoff(mapping_file, args.beta)

    ids_to_recs = SeqIO.index(script_path + "/data/uniprot_sprot.fsa", "fasta")

    for i,seq in enumerate(sequences):
        
        if verbose: 
            print( "[DETECT]: Analyzing {} ({}/{}) ...".format(seq.name(), i + 1, len(sequences)))

        identification = Identification(seq.name())
        identification.hypotheses = run_pair_alignment(seq, blast_db,num_threads, e_value, bit_score, ids_to_recs)
        
        if not identification.hypotheses:
            if verbose: 
                print( "[DETECT]: No BLASTp hits for {}".format(seq.name()))

            continue

        if verbose: 
            print( "[DETECT]: Running density estimations for {} ...".format(seq.name()))
        for hypothesis in identification.hypotheses:

            probability = calculate_probability(hypothesis, connection)
            if not (hypothesis.ec == "unknown" or hypothesis.ec == "2.7.11.1" or hypothesis.ec == "2.7.7.6" or hypothesis.ec == "2.7.13.3"):
                identification.predictions[hypothesis.ec] *= (1.0-probability)  
                identification.prediction_count[hypothesis.ec] += 1
        
        for ec,probability in identification.predictions.items():
            cumulative = 1.0 - probability
            if (cumulative > zero_density):
                identification.predictions[ec] = cumulative
            else:
                del identification.predictions[ec]
        
        if (args.top_predictions_file or args.fbeta_file):
            #sort
            identification.predictions = OrderedDict(sorted(identification.predictions.iteritems(), key=itemgetter(1), reverse=True))
        
        if (args.top_predictions_file):
            top_predictions = list()

        if (args.fbeta_file):
            fbeta_predictions = list()

        for ec in identification.predictions:
            identification_entry = "{seq_id}\t{pred_ec}\t{prob:.3e}\t{pos_hits}\t{neg_hits}\n".format(
                        seq_id=identification.query_id,
                        pred_ec=ec,
                        prob=identification.predictions[ec],
                        pos_hits=identification.prediction_count[ec],
                        neg_hits= len(identification.hypotheses)-identification.prediction_count[ec])

            
            if args.top_predictions_file and identification.predictions[ec] > probability_cutoff and len(top_predictions) < top_predictions_count:
                top_predictions.append(identification_entry)

            if args.fbeta_file and identification.predictions[ec] > ec_to_cutoff[ec]:
                fbeta_predictions.append(identification_entry)

            output.write(identification_entry)
            if (args.top_predictions_file):
                for entry in top_predictions:
                    top_predictions_file.write(entry)

            if (args.fbeta_file):
                for entry in fbeta_predictions:
                    fbeta_file.write(entry)

    if (args.top_predictions_file):
        top_predictions_file.close()

    if (args.fbeta_file):
        fbeta_file.close()

    output.close()
    connection.close()