import ConfigParser
import datetime
from collections import defaultdict
import gzip
import glob
import itertools
import io
import os
import subprocess

def log_stuff(stuff,logger):
    ''' Log given string with time
    :param str stuff: string to log
    '''
    logger.info("{time}:{stuff}".format(time=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),stuff=stuff))

def parse_config_file(samples_cfg):
    ''' Parse a config file. This is specific
    to the file being used for the 16S pipeline
    '''
    parser = ConfigParser.ConfigParser()
    parser.read(samples_cfg)
    for section in parser.sections():
        sample_name = section
        R1_fastq = parser.get(section,'R1_fastq')
        R2_fastq = parser.get(section,'R2_fastq')
        yield (sample_name,R1_fastq,R2_fastq)

def compress_fastqs(output_dir,samples_cfg):
    ''' Compress fastq files using pigz
    :param str output_dir  : base output directory
    :param str samples_cfg : config file for samples
    '''
    pigz_cmd = "pigz {fname}"
    for sample_name,R1_fastq,R2_fastq in parse_config_file(samples_cfg):
        sample_dir = os.path.join(output_dir,sample_name)
        fastqs_to_gzip = glob.glob(os.path.join(sample_dir,"*.fastq"))        
    for f in fastqs_to_gzip:
        subprocess.check_call(pigz_cmd.format(fname=f),shell=True)


def accumulate_metrics(metric_files,outfile):
    ''' Accumulate read level metrics
    :param dict metric_files : Metric files keyed by sample name
    :param str  outfile      : The output file path
    '''
    metric_names = []
    samples = []
    first = True
    metric_accumulate = defaultdict(lambda:defaultdict(int))
    
    for sample_name in metric_files:
        with open(metric_files[sample_name],"r") as IN:
            samples.append(sample_name)
            for line in IN:
                metric,value = line.strip("\n").split("\t")
                if first:
                    metric_names.append(metric)
                metric_accumulate[metric][sample_name]+=int(value)
        first = False
    header = "Metric\t{samples}\n".format(samples="\t".join(samples))
    with open(outfile,"w") as OUT:
        OUT.write(header)
        for metric in metric_names:
            out = []
            out.append(metric)
            for sample_name in samples:
                out.append(str(metric_accumulate[metric][sample_name]))
            out.append("\n")
            OUT.write("\t".join(out))

def accumulate_counts(count_files,outfile):
    ''' Accumualate primer level counts
    :param dict count_files : primer count files keyed by sample name
    :param str  outfile     : The output file path
    '''
    primer_counts = defaultdict(lambda:defaultdict(int))
    samples = []
    primer_names = []
    amplicons = []
    primer_sequences = []
    first = True
    for sample_name in count_files:
        with open(count_files[sample_name],"r") as IN:
            samples.append(sample_name)
            for line in IN:
                if line.startswith("#"): # header
                    sample_name_in_file = line.strip("\n").strip("#")
                    assert sample_name == sample_name_in_file,"Mixup in sample names !" # sanity check
                    continue
                amplicon,primer_name,primer_sequence,count = line.strip("\n").split("\t")
                primer_counts[primer_name][sample_name] = count
                if first:
                    primer_names.append(primer_name)
                    amplicons.append(amplicon)
                    primer_sequences.append(primer_sequence)
        first = False

    header = "Region\tPrimerName\tPrimerSequence\t{samples}\n".format(samples="\t".join(samples))
    with open(outfile,"w") as OUT:
        OUT.write(header)
        for i,primer_name in enumerate(primer_names):
            amplicon = amplicons[i]
            primer_seq = primer_sequences[i]
            out = []
            out.extend([amplicon,primer_name,primer_seq])
            for sample_name in samples:                
                out.append(primer_counts[primer_name][sample_name])
            out.append("\n")
            OUT.write("\t".join(out))
        
def aggregate_metrics(output_dir,samples_cfg):
    ''' Aggregate metric files
    :param str output_dir  : base output directory
    :param str samples_cfg : config file for samples
    '''
    read_metric_files = {}
    primer_metric_files = {}    
    for sample_name,R1_fastq,R2_fastq in parse_config_file(samples_cfg):
        sample_dir = os.path.join(output_dir,sample_name)
        read_metric_files[sample_name] = os.path.join(sample_dir,"{sample_name}.metrics.txt".format(sample_name=sample_name))
        primer_metric_files[sample_name] = os.path.join(sample_dir,"{sample_name}.primer_counts.txt".format(sample_name=sample_name))
    accumulate_metrics(read_metric_files,os.path.join(output_dir,"QIAseq16S.ITS.summary.txt"))
    accumulate_counts(primer_metric_files,os.path.join(output_dir,"QIAseq16S.ITS.primer.counts.txt"))
    
def open_by_magic(filename):
    '''
    Adapted from : http://stackoverflow.com/questions/18367511/how-do-i-automatically-handle-decompression-when-reading-a-file-in-python
    with modifications
    Uses the initial bytes of a file to detect the file compression.

    :param str filename: path to the input file
    :return: the appropriate file handle for reading
    :rtype: file object
    '''

    ## Add more magic strs here for various compressions
    magic_dict = {"\x1f\x8b\x08":gzip.open}
    max_len = max(len(x) for x in magic_dict)
    with open(filename) as f:
        file_start = f.read(max_len)
        for magic,fn in magic_dict.items():
            if file_start.startswith(magic):
                return io.BufferedReader(fn(filename))
            return open(filename,'r') ## Otherwise just a regular file

        
def grouper(iterable,n=10000000):
    ''' 
    Groups and returns n elements of an iterable at a time

    :param iterator iterable : The iterator to chop up
    :param int n: the chunks to group the iterator into
    :yields: iterator, i.e. the n1,n2,...n elements of the iterable
    '''
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it,n))
        if not chunk:
            return
        yield chunk    
