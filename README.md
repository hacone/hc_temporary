### General
Here are code fragments for processing short/long reads in centromeric reads.

format of pickled data

```
EncodedRead = namedtuple("EncodedRead", ("name", "mons", "length")) # string, list(AssignedMonomer), Int 
	AssignedMonomer = namedtuple("AssignedMonomer", ("begin", "end", "monomer")) # Int, Int, Monomer
	Monomer = namedtuple("Monomer", ("name", "snvs")) # string, list(SNV) or id: Int, list(SNV)
	SNV = namedtuple("SNV", ("pos", "base")) # Int, char

variant_sites = { "MonomerName" : { "Monomer" : freq_as_mon, (pos, base) : freq } }
```

### Workflow

Assume, monomers were already defined, and alignment was done for centromeric long reads, thus you have SAM files.
`EncodedRead.py` is then used to encode reads properly (i.e., assign monomers onto regions in reads, in a manner almost-non-overlapping). `EncodedRead.py` also defines related data, so it's referenced by other scripts later.

Read segregation (into distinctive arrays; n.e. chromosome segregation) based on assigned monomers, is then done using (TODO: which script?).

Then you use a collection of commands in `HOR_work.sh`; you may want to see k-monomer frequencies, and define some identification of similar monomers and then tentative HOR unit... You will HOR-encode the reads by one of the command there.

HOR encoded reads are then processed by the latest `Layout.py` to make/investigate layouts.

### TODO
Too many to be written down here.

### Data Source

- Monomer database (edit distance matrix is in local)
https://raw.githubusercontent.com/volkansevim/alpha-CENTAURI/nhmmscan/example/MigaKH.HigherOrderRptMon.fa

- squeakr by splatlab for k-mer counting, which I extended in read_squeakr/
https://github.com/splatlab/squeakr

### Requirements

`HOR_work.sh shor-hor` uses `fasta_formatter` from FASTX-Toolkit.

Our `venv` is not maintained by git, but you should be able to replicate the environment as follows.

```
python3 -m venv venv
source venv/bin/activate
pip3 install -r requirements.txt
```

As you expect, `requirements.txt` is the output from `pip3 freeze`.

```requirements.txt
biopython==1.70
cycler==0.10.0
kiwisolver==1.0.1
matplotlib==2.2.2
numpy==1.14.2
pyparsing==2.2.0
pysam==0.14
python-dateutil==2.7.2
pytz==2018.3
scikit-learn==0.19.1
scipy==1.0.1
six==1.11.0
```
