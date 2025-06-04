# fisher-metabarcode

fisher-metabarcode is the in-house script used by Fisher Lab to perform ITS sequence analysis on metabarcode samples.

## Requirements

### Software

- R (≥ 4.0)
- Internet connection

### R Packages

The following R packages are required:

- `biomaRt`
- `seqinr`
- `optparse`

You can install them using the following commands in R:

```r
install.packages("optparse")
install.packages("seqinr")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("biomaRt")
```

## Installation

fisher-metabarcode can be "installed" by simply downloading the set of scripts inside ```scripts``` .


## Usage

fisher-metabarcode can be run using our test dataset which can be downloaded from [Zenodo](https://zenodo.org/records/15594328).
Unpack and save the demo data in a convenient location that you can use for the remaining usage tutorial.
Install relevant pages using the ```.yaml``` files saved in ```envs```.

```bash
conda env create -f envs/fisher-metabarcode.yaml
```

fisher-metabarcode comes with a Bash script wrapper can be used on the command line like so:

```bash
bash fisher-metabarcode.sh -f your-sequence-input.tsv -o your-sequence-metadata.tsv
```

If you are using this for a different species/strain you can supply the appropriate dataset using the following options (in full transparency I have not extensively tested this):

```bash
# Aspergillus fumigatus Af293
Bash how to run XXX
```

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.


## Acknowledgements

## License

[MIT](https://choosealicense.com/licenses/mit/)

