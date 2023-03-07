import os
import sys
import logging
import luigi

# Internal
import core.utils
import core.demux

# Set up logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler(sys.stdout)
logger.addHandler(ch)


class config(luigi.Config):
    ''' Initialize values from configuration file
    '''
    primer_file = luigi.Parameter(description="The primer file used for the experiment")
    num_cores = luigi.IntParameter(description="Number of cores to use for primer finding")


class MyExtTask(luigi.ExternalTask):
    ''' Checks whether the file specified exists on disk
    '''
    file_loc = luigi.Parameter()
    def output(self):
        return luigi.LocalTarget(self.file_loc)


class TrimPrimersAndDemux(luigi.Task):
    ''' Task for identifying primers, trimming them and
    demultiplexing reads by amplicon
    '''
    # Parameters
    R1_fastq = luigi.Parameter()
    R2_fastq = luigi.Parameter()
    output_dir = luigi.Parameter()
    sample_name = luigi.Parameter()

    def __init__(self,*args,**kwargs):
        ''' class constructor
        '''
        super(TrimPrimersAndDemux,self).__init__(*args,**kwargs)
        self.verification_dir = os.path.join(self.output_dir,"verification")
        self.verification_file = os.path.join(self.verification_dir,"{task}.{sample}.verification.txt".format(
            task=self.__class__.__name__,sample=self.sample_name))


    def requires(self):
        ''' task dependency
        R1 and R2 fastq file must be present
        '''
        yield MyExtTask(self.R1_fastq)
        yield MyExtTask(self.R2_fastq)


    def run(self):
        ''' work to be done
        run the primer trimming and demux function
        '''
        core.demux.main(self.sample_name,os.path.join(self.output_dir,self.sample_name),self.R1_fastq,self.R2_fastq,config().primer_file,0,config().num_cores,load_cache=False,cache_file=None)
        with open(self.verification_file,"w") as OUT:
            OUT.write("task_verified\n")


    def output(self):
        ''' output from this task
        Check for existence of the verification file
        '''
        return luigi.LocalTarget(self.verification_file)


class AggregateResults(luigi.Task):
    ''' Task to aggregate metrics and create archive
    '''
    # Parameters
    output_dir = luigi.Parameter()
    samples_cfg = luigi.Parameter()

    def __init__(self,*args,**kwargs):
        ''' class constructor
        '''
        super(AggregateResults,self).__init__(*args,**kwargs)
        self.verification_dir = os.path.join(self.output_dir,"verification")
        if not os.path.exists(self.verification_dir):
            os.makedirs(self.verification_dir)
        self.verification_file = os.path.join(self.verification_dir,"{task}.verification.txt".format(
            task=self.__class__.__name__))


    def requires(self):
        ''' task dependency
        TrimPrimersAndDemux task should be done for all samples
        '''
        dependencies = []
        for sample_name,R1_fastq,R2_fastq in core.utils.parse_config_file(self.samples_cfg):
            dependencies.append(
                TrimPrimersAndDemux(
                    R1_fastq=R1_fastq,R2_fastq=R2_fastq,
                    output_dir=self.output_dir,sample_name=sample_name
                ))
        yield dependencies


    def run(self):
        ''' work to be done
        aggregate metrics , use pigz to compress fastqs and zip to gather files into an archive
        '''
        core.utils.compress_fastqs(self.output_dir,self.samples_cfg,config().num_cores)
        core.utils.aggregate_metrics(self.output_dir,self.samples_cfg)
        with open(self.verification_file,"w") as OUT:
            OUT.write("task_verified\n")


    def output(self):
        ''' output from this task
        Check for existence of the verification file
        '''
        return luigi.LocalTarget(self.verification_file)
