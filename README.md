# HIVprofiler

An interactive tool to recommend appropiate antiretroviral treatment based on genomic HIV subtype sequences. It is implemmented in python and relies on LANL databases for HIV subtype assignation. Binding DB is fetched to obtain binding affinity on th viral target. 

## Current version

v1.0

## Installation

HIVprofiler runs properly under Python 3.10

### Dependencies 

HIV_profiler requires the following programs and packages:

- pandas
- pathlib
- subprocess
- json
- rdkit
- plotly
- flash
- dask

We recommend installing the program under a conda environment:

```
bash
conda create -n deeparg_env python=3.10
source activate deeparg_env
```
Download the data required by the program available in this folder.


## Usage

To start the program run the app.py in your terminal. Click the local host link that displays the GUI.
Then specify the directory where your fasta file is located.

## Input Format

Input files must be formatted as a fasta file.

```
bash
>seq_name
ATTTTGGGGATTGAGCCGCGCGGG
```

## Theory


## Limitations
What users should know…


