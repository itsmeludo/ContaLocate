# ContaLocate
Bioinformatics / Genomics 
A tool to extract DNA segments with homogeneous olignonucleotide composition from a genome assembly.
* Learns a compositional profile for the host and the contaminant, previoulsy identified with i.e. [PhylOligo](https://github.com/itsmeludo/PhylOligo).
* Generates a GFF3 map of the contaminant positions in the genome.
 
 
Once you have explored your assembly's oligonucleotide composition, identified and selected -potentially partial- contaminant material, use ContaLocate to target species-specific contaminant DNA according to a double parametrical threshold.

```bash
contalocate.R -i genome.fasta -r genome_host.fa -c genome_conta_1.fa 
```


The training set for the host genome can be omitted if the amount of contaminant is negligible. In this case, the profile of the host will be trained on the whole genome, including the contaminant.
```bash
contalocate.R -i genome.fasta -c genome_conta_1.fa 
```


The set up of the thresholds can be manually enforced. The user will interactively prompted to set the thresholds given the distribution of windows divergence.
```bash
contalocate.R -i genome.fasta -c genome_conta_1.fa -m
```

Parameters:
* -i    --genome            Multifasta of the genome assembly (required)
* -r    --host_learn        Host training set (optional)
* -c    --conta_learn       Contaminant training set (optional, but recommended^^) if missing and sliding window parameters are given, the sliding windows composition will be compared to the whole genome composition to contrast potential HGTs (prokaryotes and simple eukaryotes only)
* -t    --win_step          Step of the sliding windows analysis to locate the contaminant (optional) default: 500bp or 100bp
* -w    --win_size          Length of the sliding window to locate the contaminant (optional) default: 5000bp 
* -d    --dist              Divergence metric used to compare profiles: (KL), JSD or Eucl
* -m    --manual_threshold  You will be asked to manually set the thresholds
* -h    --help              What it says





Install
-------

* Dependencies:
    * Python 3.x
        * [BioPython](biopython.org)
        * [Numpy](numpy.org)
        * [Cython](http://cython.org/)
    * R 3.x
        * [ape](http://ape-package.ird.fr/)
        * [getopt](https://cran.r-project.org/web/packages/getopt/getopt.pdf)
        * [gtools](https://cran.r-project.org/web/packages/gtools/index.html)

* Install python3 and the latest R version [according to your system](https://xkcd.com/1654/) 

In the Bash/Shell, as root/admin if wanted installed globally.
```Bash
#Ubuntu/Debian-based
apt-get install python3-dev python3-setuptools r-base
easy_install3 -U setuptools
pip3 install biopython 
pip3 install cython
pip3 install numpy
```

in R, as root or user
```R
install.packages(c("ape","getopt","gtools"))
```

* clone the repo

```Bash
https://github.com/itsmeludo/ContaLocate.git
```

* Link the programs into a directory listed in your $PATH

```Bash
cd ContaLocate
sudo ln -s `pwd`/{contalocate.R,Kount.py} /usr/local/bin/
```