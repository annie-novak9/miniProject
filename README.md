# miniProject
---
## Languages and Import Statements for Packages
---
**Python**
```
import os
import csv
import argparse
from Bio import Entrez
from Bio import SeqIO
```
**R**
```
library(sleuth)
library(dplyr)
```

## Software Required ##
---
- Kallisto
    1. Manual:  [Link](https://pachterlab.github.io/kallisto/manual)
    2. Installation:  [Link](https://pachterlab.github.io/kallisto/download)
- Sleuth
    1. Manual:  [Link](https://pachterlab.github.io/sleuth/manual)
    2. Installation:  [Link](https://pachterlab.github.io/sleuth/download)
- Bowtie2
    1. Manual:  [Link](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
    2. Installation:  [Link](https://sites.google.com/site/wiki4metagenomics/tools/bowtie2/install)
- SPAdes
    1. Manual:  [Link](https://cab.spbu.ru/files/release3.13.0/manual.html)
    2. Installation:  Click '2. Installation' in manual link above
- Blast+
    1. Manual:  [Link](https://www.ncbi.nlm.nih.gov/books/NBK279691/)
    2. Installation:  [Link](https://www.ncbi.nlm.nih.gov/books/NBK279671/)
 
 ## To Run ##
 ---
* To run pipeline from your command line, clone this repository to your workspace.
    ```
    git clone https://github.com/annie-novak9/miniProject.git
    ```

* Then you need to manually change your directory into the miniProject directory you just cloned using this command:
    ```
    cd ~/miniProject
    ```

* When running script in command line, the pipeline script takes one required argument 'dataset'. The options are either 'full' (in order to test the full set of input reads) or 'test' (which only runs the reduced input read files and runs in about 5 minutes).

    **Example to run full set:**
    ```
    python3 pipeline.py full
    ```
    **Example to run test set:**
    ```
    python3 pipeline.py test
    ```

* When running **full** set, all outputs will be written to the miniProject_Annie_Novak folder.

* When running **test** set, all outputs will be written to the test_outputs folder.
