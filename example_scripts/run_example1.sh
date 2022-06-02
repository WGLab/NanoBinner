#!/bin/bash

MINIMAP2=path/to/minimap2

NANO_BINNER=path/to/NanoBinner/nanoBinner.py

BARCODE_FASTA=path/to/example1/example1.barcodes.fasta

AMP_SEQ_FASTA=path/to/example1/example1.amplicon_seq.fasta

FASTQ_FILE=path/to/example1/example1.fastq.gz

OUT_DIR=./example1_output/

EXP_NAME=example1 # This is the prefix of output files

time python $NANO_BINNER \
  --in_fq $FASTQ_FILE \
  --amp_seq_fasta $AMP_SEQ_FASTA \
  --fwd_barcode_fasta $BARCODE_FASTA \
  --rev_barcode_fasta $BARCODE_FASTA \
  --require_two_barcodes \
  --num_threads 8 \
  --minimap2 $MINIMAP2 \
  --exp_name $EXP_NAME \
  --out_dir $OUT_DIR

