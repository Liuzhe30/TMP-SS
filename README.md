# TMP-SS
A Deep Learning-Based Predictor for Secondary Structure & Topology Structure Prediction of Alpha-helical Transmembrane Proteins.

<p align="center"><img width="100%" src="images/pipeline.png" /></p>

## Download data
We provide the test dataset used in this study,  you can download TEST.fasta to evaluate our method.

## Quick Start

### Requirements
- Python ≥ 3.6
- Tensorflow and Keras
- HH-suite for generating HHblits files (with the file suffix of .hhm)

### Testing & Evaluation in Command Line
We provide run.py that is able to run pre-trained models. Run it with:
```
python run.py -f sample/sample.fasta -p sample/hhblits/ -o results/
```

* To set the path of fasta file, use `--fasta` or `-f`.
* To set the path of generated HHblits files, use `--hhblits_path` or `-p`.
* To save outputs to a directory, use `--output` or `-o`.

## Progress
- [x] README for running TMP-SS.