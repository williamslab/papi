# PAPI

PAPI is a tool for inferring the admixture proportions and admixture times for the parents of an admixed sample given unphased local ancestry tracts.

## Installation

We recommend using anaconda to create a virtual environment using the included papi.yaml file
```
conda env create -f papi.yaml
```

## Usage
Clone and navigate to the papi directory

```
cd src
source activate papi
usage: inference.py [-h] --inputfile INPUTFILE --ind IND [--tracefile [TRACEFILE]] --outfile OUTFILE [--mode [MODE]] [--typ [TYP]] [--err]

optional arguments:
  -h, --help            show this help message and exit
  --inputfile INPUTFILE, -i INPUTFILE
                        input tracts file
  --ind IND, -ind IND   individual on which to run papi
  --tracefile [TRACEFILE], -t [TRACEFILE]
                        optional argument, used to store trace output when run in mcmc mode. If not provided, MCMC solver will find MAP estimate
  --outfile OUTFILE, -o OUTFILE
                        required default output file
  --mode [MODE], -m [MODE]
                        inference mode-'pymc' or 'scipy-optimize'
  --typ [TYP], -typ [TYP]
                        model to use-'bino','hmm', or 'full'
  --err, -err
```

## Input

## Options

## Output




