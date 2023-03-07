# QiaSEQ 16S

## Quick Start (Interactive) Example

Pull the latest docker image and launch an interactive shell

```bash
sudo docker pull qiaseq/qiaseq-16s
sudo docker run -i qiaseq/qiaseq-16s sh
```

Clone the 16S repository and setup luigi

```bash
cd /home/qiagen/
git clone https://github.com/qiaseq/qiaseq-16S.git

luigid --background &
export LUIGI_CONFIG_PATH=/home/qiagen/qiaseq-16S/pipeline.cfg
```

Running the code on the example fastq

```bash
cd qiaseq-16S/
export PYTHONPATH=$PYTHONPATH:""
luigi --module workflow AggregateResults --samples-cfg samples.cfg --output-dir /home/qiagen/data/output/ --workers 1 --worker-wait-interval 30
```

### Output Files

|||
|-|-|
QIAseq16S.ITS.primer.counts.txt |   Read counts for each primer
QIAseq16S.ITS.summary.txt       |   A summary file
S26_S2/*.fastq.gz               |   The demultiplexed fastq files

*Remember to modify `samples.cfg` to match your fastq files*

You can tweak the num_cores parameter in the pipeline.cfg file and the --workers flag to tweak parallelization.

## Non-Interactive Runs

Create an appropriate samples.cfg file  
It is assumed that samples.cfg and your fastq files are in the folder /home/foo/my_fav_folder/ :

```bash
sudo docker run -i -v /home/foo/my_fav_folder/:/home/foo/my_fav_folder/ qiaseq/qiaseq-16S /bin/sh -c "luigid --background & cd /home/qiagen/; git clone https://github.com/qiaseq/qiaseq-16S.git; cd qiaseq-16S; export LUIGI_CONFIG_PATH=/home/qiagen/qiaseq-16S/pipeline.cfg; export PYTHONPATH=$PYTHONPATH:""; luigi --module workflow AggregateResults --samples-cfg /home/foo/my_fav_folder/samples.cfg --output-dir /home/foo/my_fav_folder/output/ --workers 1 --worker-wait-interval 30;"
```

## Options

The `--workers` flag and the `num_cores` parameter in *pipeline.cfg* can be tuned to your parallelization needs

## Dependencies

[pigz](https://zlib.net/pigz/) - Parallel gzip  
[edlib](https://github.com/Martinsos/edlib) - Levenshtein distance sequence alignment  
[luigi](https://github.com/spotify/luigi) - Pipeline and job management interface
