# qiaseq-16S

### Quick Start Example
```
# get the docker image
sudo docker pull qiaseq/qiaseq-16s
# start an interactive shell
sudo docker run -i qiaseq/qiaseq-16s sh
# get code
cd /home/qiagen/
git clone https://github.com/qiaseq/qiaseq-16S.git
# setup luigi
luigid --background &
export LUIGI_CONFIG_PATH=/home/qiagen/qiaseq-16S/pipeline.cfg
# run code on the example fastq
cd qiaseq-16S/
export PYTHONPATH=$PYTHONPATH:""
luigi --module workflow AggregateResults --samples-cfg samples.cfg --output-dir /home/qiagen/data/output/ --workers 1 --worker-wait-interval 30

# Output files produced are : 
QIAseq16S.ITS.primer.counts.txt : Read counts for each primer
QIAseq16S.ITS.summary.txt : A summary file
S26_S2/*.fastq.gz : The demultiplexed fastq files
```

Please modify the samples.cfg file to match your own fastq files.

You can tweak the num_cores parameter in the pipeline.cfg file and the --workers flag to tweak parallelization.

To run this code in a non interactive way , create an appropriate samples.cfg file , assuming samples.cfg and your input fastq files are in the folder /home/foo/my_fav_folder/ :
```
sudo docker run -i -v /home/foo/my_fav_folder/:/home/foo/my_fav_folder/ qiaseq/qiaseq-16S \
/bin/sh -c "luigid --background & cd /home/qiagen/; git clone https://github.com/qiaseq/qiaseq-16S.git; \
cd qiaseq-16S; export LUIGI_CONFIG_PATH=/home/qiagen/qiaseq-16S/pipeline.cfg; export PYTHONPATH=$PYTHONPATH:""; \
luigi --module workflow AggregateResults --samples-cfg /home/foo/my_fav_folder/samples.cfg \
--output-dir /home/foo/my_fav_folder/output/ --workers 1 --worker-wait-interval 30;"
```

You can run this code on your own machine as well since the dependencies are pretty light. You would need :

pigz : https://zlib.net/pigz/

pip install edlib luigi
