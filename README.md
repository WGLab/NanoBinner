# NanoBinner: A demultiplexing tool for Oxford Nanopore long-read sequencing data. 

## Features
- NanoBinner uses the sequence upstream/downstream of the barcode (i.e. anchor) to help locate barcode position and eliminates random matching due to sequencing error.
- NanoBinner does not require a pretrained model for demultiplexing and supports officially provided barcodes as well as custom-designed barcodes.


## Table of Contents

- [Requirements](#Requirements)
- [Installation](#Installation)
- [Usage](#Usage)
  - [Quick start](#Quick_start)
  - [Demultiplexing regular amplicons](#Demultiplexing_regular_amplicons)
    - [Case 1: The barcode is next to the forward primer](#case1)
    - [Case 2: The barcode is next to the reverse primer](#case2)
    - [Case 3: The barcodes are on both ends. One sample have the same barcodes on both ends. Only one barcode is required to bin the reads.](#case3)
    - [Case 4: The barcodes are on both ends. One sample may or may not have the same barcodes on both ends. Two barcodes are required to bin the reads.](#case4)
- [Citation](#Citation)

## <a name="Requirements"></a>Requirements
- Operating system: Linux or macOS
- [Python](https://www.python.org/) 3.6 or later
- [minimap2](https://github.com/lh3/minimap2) 2.8 or later

## <a name="Installation"></a>Installation

NanoBinner calls `minimap2` to do sequence alignment. If you don't have `minimap2` in your system, you can install it following the instructions [here](https://github.com/lh3/minimap2#install).  
If you are using Linux, you can acquire precompiled binaries using the following commands:

```
wget https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2
tar -jxvf minimap2-2.24_x64-linux.tar.bz2
./minimap2-2.24_x64-linux/minimap2
```

Next, you can clone the repository of NanoBinner using the following command.
```
git clone https://github.com/WGLab/NanoBinner.git
```
The scripts in the `./NanoBinner` can run directly without additional compilation or installation.

## <a name="Usage"></a>Usage


### <a name="Quick_start"></a> Quick start

```
# The barcode is beside forward primer
path/to/NanoBinner/nanoBinner.py --in_fq example_data.fastq.gz --amp_seq_fasta example_amplicon_seq.fasta --out_dir . --exp_name testing --num_threads 4 --fwd_barcode_fasta example_barcodes.fasta --minimap2 path/to/minimap2

# The barcode is beside reverse primer
path/to/NanoBinner/nanoBinner.py --in_fq example_data.fastq.gz --amp_seq_fasta example_amplicon_seq.fasta --out_dir . --exp_name testing --num_threads 4 --rev_barcode_fasta example_barcodes.fasta --minimap2 path/to/minimap2

# The barcodes are on both ends. One sample have the same barcodes on both ends. Only one barcode is required to bin the reads.
path/to/NanoBinner/nanoBinner.py --in_fq example_data.fastq.gz --amp_seq_fasta example_amplicon_seq.fasta --out_dir . --exp_name testing --num_threads 4 --fwd_barcode_fasta example_barcodes.fasta --rev_barcode_fasta example_barcodes.fasta --minimap2 path/to/minimap2

# The barcodes are on both ends. One sample may or may not have the same barcodes on both ends. Two barcodes are required to bin the reads.
path/to/NanoBinner/nanoBinner.py --in_fq example_data.fastq.gz --amp_seq_fasta example_amplicon_seq.fasta --out_dir . --exp_name testing --num_threads 4 --fwd_barcode_fasta example_barcodes.fasta --rev_barcode_fasta example_barcodes.fasta --require_two_barcodes --minimap2 path/to/minimap2

# Input DNA is from a 10X Genomics single cell library
/home/fangl/NanoBinner/nanoBinner_10X.py --in_fq example.fastq.gz --barcode_list barcodes.txt --barcode_upstream_seq AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT --out_prefix testing --num_threads 8 
```

### <a name="Demultiplexing_regular_amplicons"></a> Demultiplexing regular amplicons

```
$ ./nanoBinner.py --help 
usage: nanoBinner.py [-h] [--in_fq FILE] [--in_fq_list FILE] --amp_seq_fasta
                    FILE --out_dir PATH --exp_name STRING
                    [--fwd_barcode_fasta FILE] [--rev_barcode_fasta FILE]
                    [--require_two_barcodes] [--num_threads INT]
                    [--minimap2 FILE] [--version]

A barcode demultiplexer for Oxford Nanopore long-read sequencing data

optional arguments:
  -h, --help            show this help message and exit
  --in_fq FILE          input sequencing reads in one FASTQ(.gz) file
  --in_fq_list FILE     a list file specifying all input FASTQ(.gz) files, one
                        file per line
  --amp_seq_fasta FILE  reference amplicon sequence in FASTA format
  --out_dir PATH        output directory
  --exp_name STRING     experimental name, used as prefix of output files
  --fwd_barcode_fasta FILE
                        barcode sequences of the forward primer (in FASTA
                        format)
  --rev_barcode_fasta FILE
                        barcode sequences of the reverse primer (in FASTA
                        format)
  --require_two_barcodes
                        require matched barcodes on both ends (default:
                        False). Notice: this option is valid only if both '--
                        fwd_barcode_fasta' and '--rev_barcode_fasta' are
                        supplied.
  --num_threads INT     number of threads (default: 1)
  --minimap2 FILE       path to minimap2 (default: using environment default)
  --version             show program's version number and exit

```
If you have one single input fastq file, you can supply the input with `--in_fq`. 
If you have multiple fastq files, you can supply a list file with `--in_fq_list`. The list file contains all input fastq files, one file per line. 

`--amp_seq_fasta` is the reference amplicon sequence in FASTA format. Sometimes the barcode sequence is not at the very begining of the long read. Sometimes the first a few bases of a read is truncated. Due to the sequencing error, the barcode matching is flexible and allows some mismatches. **`nanoBinner.py` assumes the reference amplicon sequence is known and uses it to distinguish amplicon sequence and barcode sequence, thus eliminates random fuzzy matching inside the amplicon.** 

`--fwd_barcode_fasta` and `--rev_barcode_fasta` are barcode sequences in FASTA format. If you use the same barcodes on both ends, you can supply `--fwd_barcode_fasta` and `--rev_barcode_fasta` with the same file. An example of `--fwd_barcode_fasta` is shown below. We supplied FASTA files of official barcodes in the `NanoBinner/ONT_barcodes` folder.  

```
>BC01
CACAAAGACACCGACAACTTTCTT
>BC02
ACAGACGACTACAAACGGAATCGA
>BC03
CCTGGTAACTGGGACACAAGACTC
>BC04
TAGGGAAACACGATAGAATCCGAA
>BC05
AAGGTTACACAAACCCTGGACAAG
```

`nanoBinner.py` supports different barcoding strategies.  
#### <a name="case1"></a> Case 1. The barcode is next to the forward primer

In this case, the amplicon structure is shown below. 
<p align="center"><img src="images/fwd_only.png" width="100%"></p>

You can use the `--fwd_barcode_fasta` argument to supply the barcode FASTA file and use the `--amp_seq_fasta` argument to supply reference amplicon FASTA file. Please note that **the `--amp_seq_fasta` file should INCLUDE the primer sequence but EXCLUDE the barcode sequence**.  An example command is shown below: 

```
/home/fangl/NanoBinner/nanoBinner.py --in_fq example_data.fastq.gz --amp_seq_fasta example_amplicon_seq.fasta --out_dir . --exp_name testing --num_threads 4 --fwd_barcode_fasta example_barcodes.fasta --minimap2 /home/fangl/software/minimap2-2.8_x64-linux/minimap2
```

#### <a name="case2"></a> Case 2. The barcode is next to the reverse primer
In this case, the amplicon structure is shown below. 
<p align="center"><img src="images/rev_only.png" width="100%"></p>

Similar to case 1, you can use the `--rev_barcode_fasta` argument to supply the barcode FASTA file and use the `--amp_seq_fasta` argument to supply reference amplicon FASTA file. Please note that **the `--amp_seq_fasta` file should INCLUDE the primer sequence but EXCLUDE the barcode sequence**.  An example command is shown below: 

```
/home/fangl/NanoBinner/nanoBinner.py --in_fq example_data.fastq.gz --amp_seq_fasta example_amplicon_seq.fasta --out_dir . --exp_name testing --num_threads 4 --rev_barcode_fasta example_barcodes.fasta --minimap2 /home/fangl/software/minimap2-2.8_x64-linux/minimap2
```

#### <a name="case3"></a> Case 3. The barcodes are on both ends. One sample have the same barcodes on both ends. Only one barcode is required to bin the reads.
This might be the most common case. In this case, the amplicon structure is shown below. 
<p align="center"><img src="images/both_sides.png" width="100%"></p>

You can supply `--fwd_barcode_fasta` and `--rev_barcode_fasta` with the same file, and use the `--amp_seq_fasta` argument to supply reference amplicon FASTA file. Please note that **the `--amp_seq_fasta` file should INCLUDE the primer sequence but EXCLUDE the barcode sequence**.  An example command is shown below: 

```
/home/fangl/NanoBinner/nanoBinner.py --in_fq example_data.fastq.gz --amp_seq_fasta example_amplicon_seq.fasta --out_dir . --exp_name testing --num_threads 4 --fwd_barcode_fasta example_barcodes.fasta --rev_barcode_fasta example_barcodes.fasta --minimap2 /home/fangl/software/minimap2-2.8_x64-linux/minimap2
```

#### <a name="case4"></a> Case 4. The barcodes are on both ends. One sample may or may not have the same barcodes on both ends. Two barcodes are required to bin the reads.
This might be the most common case. In this case, the amplicon structure is shown below. 
<p align="center"><img src="images/both_sides.png" width="100%"></p>

You can supply `--fwd_barcode_fasta` and `--rev_barcode_fasta` with the barcode FASTA file. The `--fwd_barcode_fasta` and `--rev_barcode_fasta` file may or may not be the same. **You want to use the `--require_two_barcodes` option to specify that two barcodes are required to bin the reads.**  You can use the `--amp_seq_fasta` argument to supply reference amplicon FASTA file. Please note that **the `--amp_seq_fasta` file should INCLUDE the primer sequence but EXCLUDE the barcode sequence**.  An example command is shown below: 

```
/home/fangl/NanoBinner/nanoBinner.py --in_fq example_data.fastq.gz --amp_seq_fasta example_amplicon_seq.fasta --out_dir . --exp_name testing --num_threads 4 --fwd_barcode_fasta example_barcodes.fasta --rev_barcode_fasta example_barcodes.fasta --require_two_barcodes --minimap2 /home/fangl/software/minimap2-2.8_x64-linux/minimap2
```

## <a name="Citation"></a>Citation

If you use NanoBinner, please cite:

Li Fang, Alex Mas Monteys, Alexandra Dürr, Megan Keiser, Congsheng Cheng, Akhil Harapanahalli, Pedro Gonzalez-Alegre, Beverly L. Davidson, Kai Wang. Haplotyping SNPs for allele-specific gene editing of the expanded huntingtin allele using long-read sequencing. Human Genetics and Genomics Advances, 2022, DOI: https://doi.org/10.1016/j.xhgg.2022.100146.

BibTeX format:

```
@article{FANG2022100146,
title    = {Haplotyping SNPs for allele-specific gene editing of the expanded huntingtin allele using long-read sequencing},
journal  = {Human Genetics and Genomics Advances},
pages    = {100146},
year     = {2022},
issn     = {2666-2477},
doi      = {https://doi.org/10.1016/j.xhgg.2022.100146},
url      = {https://www.sciencedirect.com/science/article/pii/S266624772200063X},
author   = {Li Fang and Alex Mas Monteys and Alexandra Dürr and Megan Keiser and Congsheng Cheng and Akhil Harapanahalli and Pedro Gonzalez-Alegre and Beverly L. Davidson and Kai Wang},
keywords = {Huntington's disease, long-read sequencing, CRISPR, SNP, repeat detection}
}
```
