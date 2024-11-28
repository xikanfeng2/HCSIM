# Tutorials

## Contents ##

1. [Overview](#overview)
    - [Implementation](#implement)
    - [CNA Types](#types)
2. [Installation](#setup)
    - [Creation of python virtual env](#env)
    - [Standard installation](#standard)
    - [Custom installation](#custom)
    - [Dependencies](#depend)
3. [Usage](#usage)
    - [Required data](#requireddata)
    - [Commands](#commands)
    - [Python Class](#class)

<a name="overview"></a>
## Overview

<a name="implement"></a>
### Implementation


<a name="types"></a>
### CNA Types

**HCSIM** have provided **nine types of hcCNAs**. The two basic hcCNA types are **deletion (DEL)**, indicating CN loss, and **duplication (DUP)**, indicating CN gain. Next, HCSIM summarizes complex hcCNA associated with **loss of heterozygosity (LOH)** with three types. LOH means the loss of one parental allele in a region, which removes genetic diversity. 1) **LOH with copy mumber loss (CNL-LOH)** means the remaining allele reduces the copy number below two, i.e. 1m0. 2) **CN-neutral LOH (CNN-LOH)** refers to a genetic state where one allele has a copy number of two, while the other allele has zero copies, resulting in a normal total copy number but the loss of heterozygosity~\cite{chen2015allele, shen2016facets}, i.e., 2m0. 3) **CN-gain LOH (CNG-LOH)** refers to the remaning allele has a copy number larger than two, e.g., 3m0, 4m0. Besides, we have the **whole chromosome loss (WCL)**, the complete loss of an entire chromosome, resulting in substantial genetic material loss, which can drive tumor progression by affecting many genes at once. We summarize the complex hcCNA associated with CN gain with two types. 1) **Gain of heterozygosity (GOH)** is an increase in the number of different alleles due to extra copies, potentially facilitating adaptation by introducing new mutations. 2) **Whole genome doubling (WGD)** is a process where the entire genome is duplicated, leading to an overall increase in copy number that promotes further chromosomal instability. Lastly, we have complex hcCNA related to haplotype status, **mirror CNAs** involve reciprocal changes in copy numbers between homologous chromosomes, e.g., the copy numbers of two genomic regions (A and B) being CN(A)=1 and CN(B)=2 in one haplotype, while in another haplotype, the copy numbers of those two bins are reversed, with CN(A)=2 and CN(B)=1.

<a name="setup"></a>
## Installation

[HCSIM](https://hcsim.readthedocs.io/en/latest/index.html) is distributed as a [PyPI](https://pypi.org/project/HCSIM/) package, and can be installed either from pip or directly from source.

<a name="env"></a>
### Creation of python virtual env

We recommend creating a virtual environment to run the HCSIM(This step is optional!). You can use the following command to create a virtual python env:

```shell
# create a python virtual env in scsilicon2 folder
python3 -m venv my_env

# activate the virtual env
source path-to-my_env/bin/activate

# deactivate the virtual env
deactivate
```

Also, you can create virtual env with `Conda`. 

```shell
# create a python virtual env in scsilicon2 folder
conda create --name myenv python=3.11

# activate the virtual env
conda activate myenv

# deactivate the virtual env
conda deactivate
```

<a name="standard"></a>
### Standard installation

HCSIM can be installed using an existing installation of `pip`.

```shell
pip install hcsim
```

<a name="custom"></a>
### Custom installation

To clone the repository and install manually, run the following from a terminal:

```shell
git clone https://github.com/xikanfeng2/HCSIM.git
cd HCSIM
python setup.py install
```

<a name="depend"></a>
### Dependencies

HCSIM is written in `python3` and requires few standard python packages and several additional standard tools.

#### > Python packages

HCSIM depends on the following standard python packages, which must be available in the python environment where the user runs HCSIM. These Python packages will be automatically installed when you install hcsim.

| Package | Tested version | Comments |
|---------|----------------|----------|
| [numpy](https://numpy.org/) | 2.1.0 | Efficient scientific computations |
| [pandas](https://pandas.pydata.org/) | 2.2.2 | Dataframe management |
| [matplotlib](https://matplotlib.org/) | 3.9.1 | Basic plotting utilities |
| [networkx](https://networkx.org) | 3.3 | Drawing cell-lineage tree

#### > Additional software

HCSIM also requires few standard additional tools, which should be included in `PATH` (e.g. by executing `export PATH=${SAMTOOLS_HOME}/bin/:${PATH}`). You can also specify the path to these tools by passing it as an argument, such as `--samtools`.

| Software | Tested version | Comments |
|----------|----------------|----------|
| [wgsim](https://github.com/lh3/wgsim) | 0.3.2 | Wgsim is a small tool for simulating sequence reads from a reference genome. |
| [SAMtools and BCFtools](http://www.htslib.org/doc/)  | 1.3.1 | Suite of programs for interacting with high-throughput sequencing data. |
| [bwa](http://bio-bwa.sourceforge.net/)  | 0.7.17-r1188 | Alignment of DNA sequencing reads. |
| [Picard Tools](https://broadinstitute.github.io/picard/)  | 2.27.5 | Picard is a set of command line tools for manipulating high-throughput sequencing (HTS) data and formats such as SAM/BAM/CRAM and VCF. |

<a name="usage"></a>
## Usage

1. [Required data](#requireddata)
2. [Commands](#commands)
3. [Python Class](#class)

<a name="requireddata"></a>
### Required data

HCSIM requires 3 input data:

1. **A reference genome file with fasta format (required).**  The most-used human reference genome are available at [GRC](https://www.ncbi.nlm.nih.gov/grc/human) or [UCSC](https://hgdownload.soe.ucsc.edu/downloads.html). Moreover, the reference genome `${REF}.fa` must be index by SAMtools (i.e. `samtools faidx ${REF}.fa`) and a dictionary must be created (`i.e. samtools dict ${REF}.fa > ${REF}.dict`) in the same folder as `${REF}`. 

2. **A list of known germline SNPs (optional).** This list is optional and SNPs can be introduced in arbitrary positions of the genome. However, we strongly suggest to obtain the lsit from known databases. One of the most used database of germline SNPs is [dbSNP](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=2388304805_pDcvxm8EAdslKzcosQNlKBXZWgUG&clade=mammal&org=Human&db=hs1&hgta_group=varRep&hgta_track=hub_3671779_censat&hgta_table=0&hgta_regionType=genome&position=chr9%3A145%2C458%2C455-145%2C495%2C201&hgta_outputType=primaryTable&hgta_outFileName=) and the common germline SNPs can be downloaded from the most recent realease [here](https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/). In order to consider only a subset of the positions to reduce the number of SNPs, the user should sample only some of the positions in this list according to estimated frequency of SNPs (typically 1 SNP every 1000 genomic positions along the entire genome. HCSIM requires the list to be given in a tab-separate file with the following fields:

| Field | Comment |
|-------|---------|
| `CHR` | Name of a chromosome, required to be present also in the reference genome |
| `POS` | A genomic position corresponding to a germline SNP in `CHR` |
| `REF_ALLELE` | Allele of the reference genome in `POS` |
| `ALT_ALLELE` | Alternate allele of the germline SNP in `POS` |

An example of this list is given [here](https://github.com/xikanfeng2/HCSIM/blob/main/example/dbsnp.tsv). If this file is not provided, HCSIM will randomly introduce SNPs into the reference genome based on the snp-ratio parameter.

3. **A list of contig to exclude (optional).**. This list is optional but highly reccomended. This is a list containing all the contigs in the given reference genome that should be excluded. An example of this list is given [here](https://github.com/xikanfeng2/HCSIM/blob/main/example/ignore.txt). HCSIM requires the list to be given in a file with every excluded contig in a new line.

<a name="commands"></a>
### HCSIM Commands

HCSIM offers different sub commands to run either the entire pipeline with all the steps or only some specific steps. In particular, the latter commands are useful when user wants to re-run some specific steps by varying some of the default parameters.
Every sub-command can be run directly when HCSIM has been correctly installed, such as `hcsim sim`.

| SubCommand | Description | Required input | Output |
|--------------|-------------|----------------|--------|
| [sim](man/hcsim-sim.md) | Running the complete HCSIM pipeline | The first required input data | [Final outputs](doc/chisel.md) |
| [gprofile](man/hcsim-gprofile.md) | Generating CNA profile | One or more running directories of previous runs of CHISEL nonormal preprocessing | [Final outputs](doc/chisel.md) |
| [gfasta](man/hcsim-gfasta.md) | Generating clone FASTA file | One or more running directories of previous runs of CHISEL nonormal preprocessing | [Final outputs](doc/chisel.md) |
| [gfastq](man/hcsim-gfastq.md) | Generating clone FASTQ file | One or more running directories of previous runs of CHISEL nonormal preprocessing | [Final outputs](doc/chisel.md) |
| [align](man/hcsim-align.md) | Aligning clone FASTQ file | One or more running directories of previous runs of CHISEL nonormal preprocessing | [Final outputs](doc/chisel.md) |
| [downsam](man/hcsim-downsam.md) | Downsampling clone BAM | One or more running directories of previous runs of CHISEL nonormal preprocessing | [Final outputs](doc/chisel.md) |
| [pbam](man/hcsim-pbam.md) | Processing cell BAMs | One or more running directories of previous runs of CHISEL nonormal preprocessing | [Final outputs](doc/chisel-calling.md) |
| [bcbam](man/hcsim-bcbam.md) | Generating barcode BAM file | One or more running directories of previous runs of CHISEL nonormal preprocessing | [Final outputs](doc/chisel-cloning.md) |


Click on the name of each command to obtain a description of all the available parameters.

## 3. Quick start
The following code runs SCSilicon.

```
import scsilicon2 as scs

# create SCSilicon2 object: ref_genome and snp_file are required, and outdir, clone_no, and cell_no are optional.
simulator = scs.SCSilicon2(ref_genome='your reference fasta file here', snp_file='your snp list file here', outdir='your output directory here', clone_no=4, cell_no=10)

# simulate dataset
simulator.sim_dataset()
```

## 4. Input file required

1. **A reference genome file with fasta format.**  
Please refer to the example fasta file `example/input/chr22.fa`.
2. **A list of SNPs.**   
The SNPs in this list can be introduced in arbitrary positions of the genome. Please refer to the example snp list file `example/input/dbsnp.tsv`.

## 5. Output files of SCSilicon2
The output directory contains three subfolders: fastq folder, fasta folder and profile folder. The structure of one example output directory is listed as follows (the clone no is 3 and the cell no is 10 in this example):

```
output
 |-fastq
 | |-normal_r2.fq
 | |-clone2
 | | |-cell0_r1.fq
 | | |-cell0_r2.fq
 | |-normal_r1.fq
 | |-clone2_r2.fq
 | |-clone1_r1.fq
 | |-clone0_r2.fq
 | |-clone0
 | | |-cell2_r1.fq
 | | |-cell3_r2.fq
 | | |-cell2_r2.fq
 | | |-cell1_r2.fq
 | | |-cell1_r1.fq
 | | |-cell3_r1.fq
 | | |-cell0_r1.fq
 | | |-cell0_r2.fq
 | |-normal
 | | |-cell2_r1.fq
 | | |-cell2_r2.fq
 | | |-cell1_r2.fq
 | | |-cell1_r1.fq
 | | |-cell0_r1.fq
 | | |-cell0_r2.fq
 | |-clone2_r1.fq
 | |-clone1
 | | |-cell1_r2.fq
 | | |-cell1_r1.fq
 | | |-cell0_r1.fq
 | | |-cell0_r2.fq
 | |-clone0_r1.fq
 | |-clone1_r2.fq
 |-fasta
 | |-clone2.fasta
 | |-normal_paternal.fasta
 | |-clone2_paternal.fasta
 | |-clone0.fasta
 | |-clone1.fasta
 | |-clone0_paternal.fasta
 | |-normal.fasta
 | |-clone1_paternal.fasta
 | |-clone2_maternal.fasta
 | |-clone1_maternal.fasta
 | |-clone0_maternal.fasta
 | |-normal_maternal.fasta
 |-profile
 | |-changes.csv
 | |-tree.pdf
 | |-maternal_cnv_matrix.csv
 | |-paternal_cnv_matrix.csv
 | |-phases.csv
 | |-cnv_profile.csv
 | |-tree.newick
```

* `fasta folder`: stores all the fasta file for each clone.

* `fastq folder`: stores all the paired-reads with fastq format for each clone and each cell.

*  `profile folder`: stores all the profile file which is related to the simulation process. The detailed explanation of the format for each file in this folder is as follows.

    1. `changes.csv`: stores the evlution path for each clone. One example is listed below:

        |Parent|Child |Haplotype|Type|Segment                |Change|
        |------|------|---------|----|-----------------------|------|
        |normal|clone0|paternal |dup |chr22:500001-1000000   |1->3  |
        |normal|clone0|maternal |del |chr22:3500001-4000000  |1->0  |
        |normal|clone0|maternal |dup |chr22:4000001-4500000  |1->2  |
        |normal|clone0|maternal |dup |chr22:5000001-5500000  |1->2  |
        |normal|clone0|maternal |dup |chr22:8000001-8500000  |1->4  |
 

    2. `cnv_profile.csv`: stores the cnv ground truth for ech clone with maternal|paternal format. One example is listed below:

        |Chromosome|Start |End     |clone0|clone1                 |clone2|
        |----------|------|--------|------|-----------------------|------|
        |chr22     |1     |500000  |1&#124;1   |3&#124;1                    |3&#124;1   |
        |chr22     |500001|1000000 |1&#124;3   |1&#124;3                    |3&#124;5   |
        |chr22     |1000001|1500000 |1&#124;1   |3&#124;2                    |3&#124;2   |
        |chr22     |1500001|3000000 |1&#124;1   |1&#124;1                    |1&#124;1   |
        |chr22     |3000001|3500000 |1&#124;1   |3&#124;2                    |3&#124;2   |
 

    3. `maternal_cnv_matrix.csv` and `paternal_cnv_matrix.csv`: store the cnv matrix of each clone seperated by maternal haplotype and paternal haplotype. One example is listed below:

        |Index|clone0_maternal_cnvs|clone1_maternal_cnvs|clone2_maternal_cnvs|
        |------|--------------------|--------------------|--------------------|
        |chr22:1-500000|1                   |3                   |3                   |
        |chr22:500001-1000000|1                   |1                   |3                   |
        |chr22:1000001-1500000|1                   |3                   |3                   |
        |chr22:1500001-3000000|1                   |1                   |1                   |
        |chr22:3000001-3500000|1                   |3                   |3                   |
    
    4. `phases.csv`: stores the SNPs in maternal|paternal haplotype. One example is listed below:

        |chr22 |16578327|1&#124;0     |
        |------|--------|--------|
        |chr22 |17307398|1&#124;0     |
        |chr22 |18025718|1&#124;0     |
        |chr22 |21416314|0&#124;1     |
        |chr22 |22418251|1&#124;0     |

    5. `tree.newick` and `tree.pdf`: the cnv elution tree with newick format and pdf format.

    The example profile folder can be found in `data/profile` folder.

## 6. `SCSilicon2` object
All the general parameters for the SCSilicon2 simulation are stored in a `SCSilicon2` object. Let’s create a new one.

```Python
simulator = scs.SCSilicon2()
```

### 6.1 All parameters in `SCSilicon2` object

* `ref_genome`: str, required<br>
    The reference genome file path
        
* `snp_file`: str, required<br>
    The snp list file

* `outdir`: str, optional, default: './'<br>
    The output directory

* `clone_no`: int, optional, default: 1<br>
    The random clone number contained in evolution tree

* `cell_no`: int, optional, default: 2<br>
    The total cell number for this simultion dataset. Please make sure the `cell_no` is large than `clone_no`. At least one cell is geneated for nomal case.

* `max_cnv_tree_depth`: int, optional, default: 4<br>
    The maximum depth of random evolution tree

* `bin_len`: int, optional, default: 500000<br>
    The fixed bin length

* `HEHOratio`: float, optional, default: 0.5<br>
    Ratio of heterozygous SNPs

* `cnv_prob_cutoff`: float, optional, default: 0.8<br>
    The cutoff probability of a bin undergoing CNV, if random probability is larger than cutoff, CNV happens

* `clone_coverage`: float, optional, default: 30<br>
    The coverage for clone fastq file

* `cell_coverage`: float, optional, default: 0.5<br>
    The coverage for each cell in a clone

* `reads_len`: int, optional, default: 150<br>
    The reads length in fastq file

* `insertion_size`: int, optional, default: 350<br>
    The outer distance between the two ends

* `error_rate`: float, optional, default: 0.02<br>
    The base error rate

