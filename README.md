# seqbot

A library for demultiplexing and sequencing automation.

`demuxer` is a script to watch local storage for new sequencing runs. If a sample-sheet can be found on S3 it will download it, split it into batches (if necessary) and demux it. If a NovaSeq run doesn't have a sample-sheet it will provide a preliminary count of the index combinations, for QC purposes.
