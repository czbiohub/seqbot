# seqbot
Scripts for sequencing automation

### `watch_flexo.py`

Currently runs on a local VM and uploads BCL files to AWS S3

### `demuxer.py`

Similarly watched the local storage for new sequencing runs. If a sample-sheet can be found on S3 it will download it, split it into batches (if necessary) and demux it. A work in progress.

