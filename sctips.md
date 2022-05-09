# Required software

## Python - https://www.python.org/downloads/
## R - https://cran.r-project.org/

### Packages for R
```
install.packages('ape')
install.packages('adegenet')
install.packages('phangorn')
install.packages('dendextend')
install.packages('seqinr')


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")
BiocManager::install("tidyverse")
BiocManager::install("treeio")
BiocManager::install("Biostrings")
```
## FastTree
```
curl -O http://www.microbesonline.org/fasttree/FastTree.c
gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm
```
## Mafft
```
https://mafft.cbrc.jp/alignment/software/source.html
gunzip -cd mafft-x.x-src.tgz | tar xfv -
cd mafft-x.x/core/

DASH_CLIENT = dash_client
make clean
make
su
make install

```
## Muscle
```
conda install muscle

```
## T-Coffee
```
1. Download the installer package corresponding to your system from:
   http://tcoffee.org/Packages/Stable/Latest/

2. Grant execution permission to the downloaded file with the following command:
   ##: chmod +x T-COFFEE_installer_"version_x".bin

3. Launch the installation wizard with:
   ##: ./T-COFFEE_installer_"version_x".bin

4. Follow the wizard instructions and complete the installation

5. Open a new terminal session to be sure that your environment is updated

6. Type the following command to verify the installation was successful:
   $$: t_coffee -version

```

## Mr Bayes
```
git clone https://github.com/NBISweden/MrBayes.git
./configure
make
```




1. Manually picked up sequences of the orthologs [genes of interest] from this database https://www.genenames.org/tools/hcop/#!/
 - organized them in .fa format
 - put all fastas in Data directory

2. Wrote a tiny Python script for alignment with different algorithms (MAFFT, T-Coffee and Muscle):
```
import os

fastas = os.listdir('path')
algothrim = 'mafft'
path = "path"
for i in fastas:
    command = "mafft --globalpair {pat}data/{gene} > {pat}/aligs/mafft/{gene}_{algo}.fa".format(pat = path, gene = i, algo=algothrim)
    print(command)
    os.system(command)


algothrim = 'muscle'
for i in fastas:
    command = "muscle -maxiters 10000 -in {pat}data/{gene} -out {pat}/aligs/muscle/{gene}_{algo}.fa".format(pat = path, gene = i, algo=algothrim)
    print(command)
    os.system(command)


    # t_coffee
path = "path"
for i in fastas:
    command = "t_coffee {pat}data/{gene}".format(pat = path, gene = i)
    print(command)
    os.system(command)
```

I compared all 3 set of alignments among each other, the t-coffee's results where significantly different and also I didn't like t-coffee format so I chose mafft output for further analysis.

3. Parsimony trees
For Parsimony algorithm I used script provided in the class, though I put it into the loop.
```
library(ape)
library(adegenet)
library(phangorn)

files <- list.files(path="path")
for (i in files){
  path <- paste("/home/ns/Documents/phil/aligs/mafft/", i, sep="")
  dna <- fasta2DNAbin(file=path)
  dna2 <- as.phyDat(dna)
  tre.ini <- nj(dist.dna(dna, model="raw"))
  parsimony(tre.ini, dna2)
  tre.pars <- optim.parsimony(tre.ini, dna2)
  out <- paste("path", i, sep="")
  out2 <- paste(out, '.pdf', sep="")
  pdf(out2)
  plot(tre.pars, cex=1.6)
  dev.off()
}

```

4. Maximum likelihood trees, instead of True ML I used approximate ML from FastTree.
```
FastTree -nt -gtr name.fa > tree
```
All trees are located in directory 'trees'
For Tree visualization I used R script:
```
library(tidyverse)
library(ggtree)
library(treeio)

tree <- read.tree(nwk)
tree <- read.tree("path_to_tree")
ggplot(tree) + geom_tree() + theme_tree() + geom_tiplab()
ggtree(tree)
```



5. Finally, for Bayesian algorithms I used MrBayes:

But before that I had to covert fasta alignmetns to Nexus format with R script:
```
library(seqinr)
library(Biostrings)

data = read.fasta("/path/18s.fa_mafft.fa")
write.nexus.data(data, file="path/18s.fa_mafft.nexus", format="dna")
```

I run MrBayes in -i mode:

```
./mb
Execute 18S.nexus # new alignment
lset nst=6
mcmc ngen=20000 samplefreq=100
sump
sumt
```
All output files were moved to tree/mafft_mb folder
