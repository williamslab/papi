# PAPI

PAPI is a tool for inferring the admixture proportions and admixture times for the parents of an admixed sample given unphased local ancestry tracts.

## Installation

We recommend using anaconda to create a virtual environment using the included file ```papi_spec-file.txt```.
```
conda create --name papi --file papi_spec-file.txt
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

An example input tracts file is provided in ```examples/tracts.txt``` that has the following structure

```
[[('10', 0.004568105), ('00', 0.46384804), ('10', 42.318381695), ('00', 27.1541), ... ], ...]
[[('10', 33.363797840000004), ('11', 18.6777), ('10', 8.2969), ('11', 14.64119999999999), ... ], ...]
...
...
...
```
Each line represents the tracts of an individual as a nested list of lists; each nested list corresponds to a chromosome. The first element of each tuple e.g ```('10', 0.004568105)``` represents a single tract that is, in this case, heterozygous for the two ancestry states ```1``` and ```0```, while the second element represents the length of the tract in centiMorgans.


## Example

Using the included tracts file, the simplest way of running PAPI is as follows:
```
python src/inference.py -i examples/tracts.txt -ind 1 -o test -m scipy-optimize -typ full
```
which will output GD estimates for the first line of the tracts file in `test.scipy.map` The parameters under which these tracts were simulated can be found in the corresponding ```headers.txt``` file.


## Output




