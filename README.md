# ContaLocate
Bioinformatics / Genomics 
A tool to extract DNA segments with homogeneous olignonucleotide composition from a genome assembly.
 * Learns a compositional profile for the host and the contaminant, previoulsy identified with i.e. [PhylOligo](https://github.com/itsmeludo/PhylOligo).
 * Generate a GFF3 map of the contaminant positions in the genome.
 
 
Once you have explored your assembly's oligonucleotide composition, identified and selected -potentially partial- contaminant material, use ContaLocate to target species-specific contaminant DNA according to a double parametrical threshold.

----------------------------------------------
```bash
contalocate.R -d Eucl -i genome.fasta -o genome.Eucl.mat -u 64
```

Regroup contigs by compositional similarity on a tree and explore branching
---------------------------------------------------------------------------

```bash
phyloselect.R -d -v  -i genome.Eucl.mat -a genome.fasta -o genome_conta
```

matrix(c(
        'genome'         , 'i', 1, "character", "file from fastq-stats -x (required)",
        'host_learn'     , 'r', 2, "character", "input gc content file (optional)",
        'conta_learn'    , 'c', 1, "character", "output filename (optional)",
        'win_step'       , 't', 2, "int", "output filename (optional)",
        'win_size'       , 'w', 2, "int", "output filename (optional)",
        'dist'           , 'd', 2, "character", "Divergence metric used to compare profiles: (KL), JSD or Eucl",
        'manual_threshold' , 'm', 0, "logical", "You will be asked to manually set the thresholds",
        'help'           , 'h', 0, "logical",   "this help"





Install
-------

Dependencies:
* Python 3.x
 * BioPython
 * numpy
 * Cython
* R 3.x
 * ape
 * getopt
 * gtools

In the Bash/Shell, as root/admin if wanted installed globally.
Install python3 and the latest R version [according to your system](https://xkcd.com/1654/) 
```Bash
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

