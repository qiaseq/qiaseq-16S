# standard library
import itertools
import logging
import functools
import sys
import os
import collections
import multiprocessing
# 3rd party
import edlib
# our modules
import utils
#import RemoteException

# Initialize some global constants
_KMER_SIZE_ = 8 # k-mer size
_SEED_SIZE_ = 30 # first _SEED_SIZE_ bases from read sequence will be used for generating kmers and looking up primer k-mer index
_PADDING_ = 5 # num bases to pad to read sequence past the primer length when computing editdistance
_MISMATCH_THRESHOLD_ = 0.10 # editdist/len(primer) ; only primers >= to this value will be considered valid for trimming
_CHUNK_SIZE_ = 40000000 # max lines to read into memory (will read this many lines into memory for each R1 and R2 fastq file, should be divisible by 4 for accurate parsing of fastq)

# Setup logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler(sys.stdout)
logger.addHandler(ch)

def parse_primer_file(primer_file):
    ''' Parse primer file
    '''

def return_kmers(seq,k=_KMER_SIZE_):
    ''' Return all k-mers for a sequence
    :param str seq : The nucleotide sequence to kmerize
    :param int k   : The k-mer size
    '''
    kmers = set(''.join(itertools.islice(seq,i,i+k)) for i in range(len(seq)+1-k))
    return kmers

def create_primer_search_datastruct(primer_file,cache=False,cache_file=None):
    ''' Create datastructures used for searching primers
    Overlapping 8-mer index for each primer
    Information for each primer including closely clustered primers from cd-hit

    :param str  primer_file : The PrimerFile for this readset
    :param bool cache       : Whether to cache the results to a file
    :param str  cache_file  : If cache is True , file path to cache to

    :rtype If cache is False, a tuple of : (defaultdict(set),defaultdict(list),dict)
    '''
    forward_primer_kmer = collections.defaultdict(set)
    forward_primers = collections.defaultdict(list)
    reverse_primer_kmer = collections.defaultdict(set)
    reverse_primers = collections.defaultdict(list)

    with open(primer_file,'r') as IN:
        for line in IN:
            contents = line.strip('\n\r').split('\t')
            if line.startswith("Region"):
                continue
            amplicon,forward_primer_name,forward_primer_id,forward_primer_seq,reverse_primer_name,reverse_primer_id,reverse_primer_seq = contents            
            # forward primer            
            primer_len = len(forward_primer_seq)
            if forward_primer_seq not in forward_primers:
                if forward_primer_seq != "":
                    info = [amplicon,forward_primer_name,forward_primer_id,primer_len]
                    forward_primers[forward_primer_seq] = info
            # create k-mer index
            for oligo in return_kmers(forward_primer_seq):
                forward_primer_kmer[oligo].add(forward_primer_seq)
            # reverse primer
            primer_len = len(reverse_primer_seq)
            if reverse_primer_seq not in reverse_primers:
                if reverse_primer_seq != "":
                    info = [amplicon,reverse_primer_name,reverse_primer_id,primer_len]
                    reverse_primers[reverse_primer_seq] = info
            # create k-mer index
            for oligo in return_kmers(reverse_primer_seq):
                reverse_primer_kmer[oligo].add(reverse_primer_seq)

    datastruct = {
        "fwd_primer_kmer":forward_primer_kmer,"fwd_primer_info":forward_primers,
        "rev_primer_kmer":reverse_primer_kmer,"rev_primer_info":reverse_primers
        }
    if cache:
        with open(cache_file,"wb") as OUT:
            pickle.dump(datastruct,OUT)
    else:
        return datastruct

def trim_primer(primer_datastruct,read):
    ''' Trim appropriate Primer from the read sequence
    :param tuple primer_datastruct  : Object returned by the function 
                                      create_primer_search_datastruct
    :param str   read               : The read sequence

    :returns : The best primer hit and the trimming site
    :rtype tuple
    '''
    best_editdist = None
    best_score = None
    best_primer = None
    best_plen = None
    best_alignment = None
    primer_kmer , primers = primer_datastruct

    candidates = set()
    for oligo in set(''.join(itertools.islice(read[0:_SEED_SIZE_],i,i+_KMER_SIZE_)) for i in range(len(read[0:_SEED_SIZE_])+1-_KMER_SIZE_)): # all 8-mers of the first 30 bp of the read(length will be truncated to read length if read is less than 30 bp)
        if oligo in primer_kmer: # check if the 8-mer is in the index
            for c in primer_kmer[oligo]: # iterate over all primers corresponding to this 8-mer
                candidates.add(c)
    if len(candidates) == 0: # no hits in the index, exhaustive search over all primers
        candidates = primers
    num_candidates = len(candidates)
    to_log = {}    
    for p in candidates:        
        p_len = primers[p][-1]
        if read[0:p_len] == p and num_candidates == 0: # exact match and no other primers to check
            amplicon = primers[p][0]
            primer_name = primers[p][1]
            return (p_len+1,p,'0',amplicon,primer_name)
        else:
            alignment = edlib.align(p,read[0:p_len + _PADDING_],mode="SHW")
            editdist = alignment['editDistance']
            score =  float(editdist)/p_len
            if best_score == None or score <= best_score:
                if score == best_score:
                    to_log[best_primer] = "*** Primers have same score *** Score : {score} : Primer1 : {p1} ; Primer2 : {p2} Read : {read}".format(score=score,p1=best_primer,p2=p,read=read)

                if best_plen == None or best_plen < p_len: # when same score with same primer length; the first hit will be chosen
                    best_alignment = alignment
                    best_primer = p
                    best_score = score
                    best_editdist = editdist                

    assert best_score != None, "Bug in logic ! Primer could not be scored correctly !"

    if best_score <= _MISMATCH_THRESHOLD_:
        if best_primer in to_log:
            utils.log_stuff(to_log[best_primer],logger)
            
        amplicon = primers[best_primer][0]
        primer_name = primers[best_primer][1]
        return (best_alignment['locations'][-1][1]+1,best_primer,str(best_editdist),amplicon,primer_name)
    else:
        return (-1,"-1","-1","-1","-1")

def iterate_fastq(R1_fastq,R2_fastq):
    ''' Iterate over R1 and R2 fastq file in chunks
    :param str R1_fastq : path to R1 fastq file
    :param str R2_fastq : path to R2 fastq file

    :yields [[[R1_record],[R2_record],[[R1_record],[R2_record]].....]
    '''
    with utils.open_by_magic(R1_fastq) as IN1 , utils.open_by_magic(R2_fastq) as IN2:
        for R1,R2 in itertools.izip(utils.grouper(IN1,_CHUNK_SIZE_),utils.grouper(IN2,_CHUNK_SIZE_)):
            to_yield = []
            count = 1
            R1_info = []
            R2_info = []
            for R1_lines,R2_lines in zip(R1,R2):
                R1_info.append(R1_lines.strip('\n'))
                R2_info.append(R2_lines.strip('\n'))
                if count%4==0:
                    to_yield.append((R1_info,R2_info))
                    count = 0
                    R1_info = []
                    R2_info = []                    
                count+=1
            yield to_yield            

def wrapper_trimmer(primer_datastruct,read_info):
    '''
    '''
    R1_info, R2_info = read_info
    R1_id,R1_seq,R1_t,R1_qual = R1_info
    R2_id,R2_seq,R2_t,R2_qual = R2_info
    R1_trim_pos,R1_primer,R1_primer_err,R1_amplicon,R1_primer_name = trim_primer((primer_datastruct["fwd_primer_kmer"],primer_datastruct["fwd_primer_info"]),R1_seq)
    R2_trim_pos,R2_primer,R2_primer_err,R2_amplicon,R2_primer_name = trim_primer((primer_datastruct["rev_primer_kmer"],primer_datastruct["rev_primer_info"]),R2_seq)

    return ((R1_id,R1_seq,R1_t,R1_qual,R1_trim_pos,R1_primer,R1_primer_err,R1_amplicon,R1_primer_name),
            (R2_id,R2_seq,R2_t,R2_qual,R2_trim_pos,R2_primer,R2_primer_err,R2_amplicon,R2_primer_name))
        
    
def main(sample_name,output_dir,R1_fastq,R2_fastq,primer_file,primer_3_bases,num_cores,load_cache=False,cache_file=None):
    ''' Main function
    :param str sample_name          : Sample Name , will be used as in primer count file
    :param str output_dir           : Output file path for demux results, will be created if it does not exist
    :param str R1_fastq             : Path to Input R1 fastq file
    :param str R2_fastq             : Path to Input R2 fastq file
    :param str primer_file          : Path to primer file for this readset
    :param int primer_3_bases       : Number of bases to keep on the 3' end of the primers
    :param int num_cores            : The number of cores to use
    :param bool load_cache          : Whether to load the primer search datastructure (useful when splitting a fastq and running in parallel)
    :param str cache_file           : If load_cache is True, file path to load from
    '''
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    num_cores = int(num_cores)
    # counters
    num_R1=0
    num_trimmed_R1=0
    num_trimmed_R2=0
    num_not_trimmed_R1=0
    num_not_trimmed_R2=0
    fwd_primer_counts = collections.defaultdict(int)
    rev_primer_counts = collections.defaultdict(int)
    overall_metrics = collections.defaultdict(int)    

    # primer search datastruct
    if load_cache: # pickle and dill do not seem to be working. will debug this later. do not use cache for now.
        #raise UserWarning("Feature not implemented yet, cannot serialize primer data structure objects.")
        with open(cache_file,"rb") as IN:
            primer_datastruct = pickle.load(IN)
            utils.log_stuff("Loaded Primer Datastructure",logger)
    else:
        primer_datastruct = create_primer_search_datastruct(primer_file,cache=False)
        utils.log_stuff("Created Primer Datastructure",logger)

    # store all primer names , fwd and rev
    fwd_primer_names = set()
    rev_primer_names = set()
    # store all amplicon names
    amplicon_names = set()
    for primer in primer_datastruct["fwd_primer_info"]:
        fwd_primer_names.add(primer_datastruct["fwd_primer_info"][primer][1])
        amplicon_names.add(primer_datastruct["fwd_primer_info"][primer][0])
    for primer in primer_datastruct["rev_primer_info"]:
        rev_primer_names.add(primer_datastruct["rev_primer_info"][primer][1])

    # dict for file handles
    OUT_R1 = {}
    OUT_R2 = {}
    
    # spawn threads
    p = multiprocessing.Pool(num_cores)
    # function to be parallelized
    func = functools.partial(wrapper_trimmer,primer_datastruct)    
    

    for chunks in iterate_fastq(R1_fastq,R2_fastq):
        for R1_info,R2_info in p.map(func,chunks):            
            R1_id,R1_seq,R1_t,R1_qual,R1_trim_pos,R1_primer,R1_primer_err,R1_amplicon,R1_primer_name = R1_info
            R2_id,R2_seq,R2_t,R2_qual,R2_trim_pos,R2_primer,R2_primer_err,R2_amplicon,R2_primer_name = R2_info            

            if R1_trim_pos == -1 or R2_trim_pos == -1: # No primer found
                overall_metrics['read fragments dropped, primer not found']+=1
                out_file = "unknown"
            elif R1_amplicon != R2_amplicon: # Found primer but not on same amplicon
                overall_metrics['read fragments dropped, R1 and R2 not on same amplicon']+=1
                outsuffix = "unknown"
            else: # Found primer and on same amplicon
                # trim R1
                if primer_3_bases ==  -1 or primer_3_bases > len(R1_seq): # keep R1 and R1_qual to be as is
                    pass
                elif primer_3_bases == 0: # trim all bases belonging to the primer
                    R1_seq = R1_seq[R1_trim_pos:]
                    R1_qual = R1_qual[R1_trim_pos:]
                else: # keep said bases belonging to the primer
                    R1_seq = R1_seq[R1_trim_pos - primer_3_bases:]
                    R1_qual = R1_qual[R1_trim_pos - primer_3_bases:]

                # trim R2
                if primer_3_bases ==  -1 or primer_3_bases > len(R2_seq): # keep R2 and R2_qual to be as is
                    pass
                elif primer_3_bases == 0: # trim all bases belonging to the primer
                    R2_seq = R2_seq[R2_trim_pos:]
                    R2_qual = R2_qual[R2_trim_pos:]
                else: # keep said bases belonging to the primer
                    R2_seq = R2_seq[R2_trim_pos - primer_3_bases:]
                    R2_qual = R2_qual[R2_trim_pos - primer_3_bases:]

                fwd_primer_counts[R1_primer_name]+=1
                rev_primer_counts[R2_primer_name]+=1                    
                outsuffix = R1_amplicon

            overall_metrics['read fragments, total']+=1
            if outsuffix not in OUT_R1: # initialize file handle 
                OUT_R1[outsuffix] = open(os.path.join(output_dir,sample_name+"_"+outsuffix+"_R1.fastq"),"w")
                OUT_R2[outsuffix] = open(os.path.join(output_dir,sample_name+"_"+outsuffix+"_R2.fastq"),"w")
            OUT_R1[outsuffix].write("{r_id}\n{r_seq}\n{plus}\n{r_qual}\n".format(r_id=R1_id,r_seq=R1_seq,plus=R1_t,r_qual=R1_qual))
            OUT_R2[outsuffix].write("{r_id}\n{r_seq}\n{plus}\n{r_qual}\n".format(r_id=R2_id,r_seq=R2_seq,plus=R2_t,r_qual=R2_qual))                
    p.close()
    p.join()
    for outsuffix in OUT_R1:
        OUT_R1[outsuffix].close()
    for outsuffix in OUT_R2:
        OUT_R2[outsuffix].close()
    
    # fwd primer counts
    with open(os.path.join(output_dir,"{sample_name}.forward_primer_counts.txt".format(sample_name=sample_name)),"w") as OUT:
        OUT.write("#{sample_name}\n".format(sample_name=sample_name))
        for primer_name in fwd_primer_names:
            if primer_name not in fwd_primer_counts:
                count = 0
            else:
                count = fwd_primer_counts[primer_name]                
            OUT.write("{primer}\t{count}\n".format(primer=primer_name,count=count))
            
    # rev primer counts
    with open(os.path.join(output_dir,"{sample_name}.reverse_primer_counts.txt".format(sample_name=sample_name)),"w") as OUT:
        OUT.write("#{sample_name}\n".format(sample_name=sample_name))
        for primer_name in rev_primer_names:
            if primer_name not in rev_primer_counts:
                count = 0
            else:
                count = rev_primer_counts[primer_name]                
            OUT.write("{primer}\t{count}\n".format(primer=primer_name,count=count))
            
    # metrics
    with open(os.path.join(output_dir,"{sample_name}.metrics.txt".format(sample_name=sample_name)),"w") as OUT:
        for met in ["read fragments, total","read fragments dropped, primer not found","read fragments dropped, R1 and R2 not on same amplicon"]:
            OUT.write("{metric}\t{value}\n".format(metric=met,value=overall_metrics[met]))

if __name__ == '__main__':
    assert len(sys.argv) == 8, "Incorrect command line params specified !"
    main(*sys.argv[1:])
