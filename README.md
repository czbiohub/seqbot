# seqbot

A library for demultiplexing and sequencing automation.

### `seqbot.bcl2fastr`

Functions for processing CBCL (NovaSeq) data files and extracting reads and quality scores, as an alternative to Illumina's `bcl2fastq` software. This is a work in progress!

### `seqbot.bcl2fastr.demuxer.py` or `demuxer` 

A script to watch local storage for new sequencing runs. If a sample-sheet can be found on S3 it will download it, split it into batches (if necessary) and demux it. If a NovaSeq run doesn't have a sample-sheet it will provide a preliminary count of the index combinations, for QC purposes. Uses reflow to process non-NovaSeq sequencing data.

### `seqbot.bcl2fastr.index_count.py` or `nova_index`

A script to return the top N most common index combinations from a NovaSeq run, by processing the CBCL files.

### `seqbot.bcl2fastr.write_fastq.py` or `nova_demux`

A script to write fastq.gz files from a NovaSeq run using a samplesheet, just like Illumina's `bcl2fastq` software but which much better support for large sample numbers.

