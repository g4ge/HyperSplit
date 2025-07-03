# Shape-Agnostic Largest Table Overlap Detection via Hypergraphs

## Packages

1. python==3.8.0
2. pandas==2.0.3
3. timeout-decorator==0.5.0

## Datasets

The two datasets WikiTables and GitTables are available [here](https://my.hidrive.com/share/6tuees3os3#$/). For more info, please see [Armadillo](https://github.com/HPI-Information-Systems/Armadillo?tab=readme-ov-file#datasets).

## File Description

1. `hypersplit.py`: Implementation of the HyperSplit algorithm for finding the shape-agnostic largest table overlap.
2. `config.py`: Configuration settings for timeout values.

## Usage

`python hypersplit.py --x <str> --y <str> [--xh <int>] [--yh <int>] [--d <float>] [--n <int>]`

- `--x <str>`: Required. CSV filename for table X.
- `--y <str>`: Required. CSV filename for table Y.
- `--xh <int>`: Optional. Number of header rows to skip in table X (default: 0).
- `--yh <int>`: Optional. Number of header rows to skip in table Y (default: 0).
- `--d <float>`: Optional. Delta value for tolerance-based pruning (default: None).
- `--n <int>`: Optional. Number of search tree levels to apply tolerance-based pruning (default: None).
