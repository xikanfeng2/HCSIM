# Tutorials

## Contents ##

1. [Overview](#overview)
    - [Implementation](#implement)
    - [CNA Types](#cna-types)
2. [Installation](#setup)
    - [Creation of python virtual env](#env)
    - [Standard installation](#standard)
    - [Custom installation](#custom)
    - [Dependencies](#depend)
3. [Usage](#usage)
    - [Required data](#requireddata)
    - [Commands](#hcsim-commands)
    - [Class](#hcsim-class)
    - [Outputs](#outputs)

<a name="overview"></a>
## Overview

<a name="implement"></a>
### Implementation


<a name="cna-types"></a>
### CNA Types

**HCSIM** have provided **nine types of hcCNAs**. The two basic hcCNA types are **deletion (DEL)**, indicating CN loss, and **duplication (DUP)**, indicating CN gain. Next, HCSIM summarizes complex hcCNA associated with **loss of heterozygosity (LOH)** with three types. LOH means the loss of one parental allele in a region, which removes genetic diversity. 1) **LOH with copy mumber loss (CNL-LOH)** means the remaining allele reduces the copy number below two, i.e. 1m0. 2) **CN-neutral LOH (CNN-LOH)** refers to a genetic state where one allele has a copy number of two, while the other allele has zero copies, resulting in a normal total copy number but the loss of heterozygosity, i.e., 2m0. 3) **CN-gain LOH (CNG-LOH)** refers to the remaning allele has a copy number larger than two, e.g., 3m0, 4m0. Besides, we have the **whole chromosome loss (WCL)**, the complete loss of an entire chromosome, resulting in substantial genetic material loss, which can drive tumor progression by affecting many genes at once. We summarize the complex hcCNA associated with CN gain with two types. 1) **Gain of heterozygosity (GOH)** is an increase in the number of different alleles due to extra copies, potentially facilitating adaptation by introducing new mutations. 2) **Whole genome doubling (WGD)** is a process where the entire genome is duplicated, leading to an overall increase in copy number that promotes further chromosomal instability. Lastly, we have complex hcCNA related to haplotype status, **mirror CNAs** involve reciprocal changes in copy numbers between homologous chromosomes, e.g., the copy numbers of two genomic regions (A and B) being CN(A)=1 and CN(B)=2 in one haplotype, while in another haplotype, the copy numbers of those two bins are reversed, with CN(A)=2 and CN(B)=1.

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
2. [Commands](#hcsim-commands)
3. [Class](#hcsim-class)
4. [Outputs](#outputs)

<a name="requireddata"></a>
### Required data

HCSIM requires 3 input data:

1. **A reference genome file with fasta format (required).**  The most-used human reference genome are available at [GRC](https://www.ncbi.nlm.nih.gov/grc/human) or [UCSC](https://hgdownload.soe.ucsc.edu/downloads.html). Moreover, the reference genome `${REF}.fa` must be index by SAMtools (i.e. `samtools faidx ${REF}.fa`) and a dictionary must be created (`i.e. samtools dict ${REF}.fa > ${REF}.dict`) in the same folder as `${REF}`. 

2. **A list of known germline SNPs (optional).** This list is optional and SNPs can be introduced in arbitrary positions of the genome. However, we strongly suggest to obtain the list from known databases. One of the most used database of germline SNPs is [dbSNP](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=2388304805_pDcvxm8EAdslKzcosQNlKBXZWgUG&clade=mammal&org=Human&db=hs1&hgta_group=varRep&hgta_track=hub_3671779_censat&hgta_table=0&hgta_regionType=genome&position=chr9%3A145%2C458%2C455-145%2C495%2C201&hgta_outputType=primaryTable&hgta_outFileName=) and the common germline SNPs can be downloaded from the most recent realease [here](https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/). In order to consider only a subset of the positions to reduce the number of SNPs, the user should sample only some of the positions in this list according to estimated frequency of SNPs (typically 1 SNP every 1000 genomic positions along the entire genome. HCSIM requires the list to be given in a tab-separate file with the following fields:

| Field | Comment |
|-------|---------|
| `CHR` | Name of a chromosome, required to be present also in the reference genome |
| `POS` | A genomic position corresponding to a germline SNP in `CHR` |
| `REF_ALLELE` | Allele of the reference genome in `POS` |
| `ALT_ALLELE` | Alternate allele of the germline SNP in `POS` |

An example of this list is given [here](https://github.com/xikanfeng2/HCSIM/blob/main/example/dbsnp.tsv). If this file is not provided, HCSIM will randomly introduce SNPs into the reference genome based on the snp-ratio parameter.

3. **A list of contig to exclude (optional).** This list is optional but highly reccomended. This is a list containing all the contigs in the given reference genome that should be excluded. An example of this list is given [here](https://github.com/xikanfeng2/HCSIM/blob/main/example/ignore.txt). HCSIM requires the list to be given in a file with every excluded contig in a new line.

<a name="hcsim-commands"></a>
### HCSIM Commands

HCSIM offers different sub commands to run either the entire pipeline with all the steps or only some specific steps. In particular, the latter commands are useful when user wants to re-run some specific steps by varying some of the default parameters. Every sub-command can be run directly when HCSIM has been correctly installed, such as `hcsim sim`.

```{note}
The complete HCSIM pipeline will sequentially execute the following modules in order: gprofile -> gfasta -> gfastq -> align -> downsam -> pbam -> bcbam. If you want to re-run one of these steps, please make sure the previous commands have been successfully executed.
```

```note
Click on the name of each command to obtain a description of all the available parameters.
```

| SubCommand | Description | Required input |
|--------------|-------------|----------------|
| [sim](man/hcsim-sim.md) | Running the complete HCSIM pipeline | a reference genome file |
| [gprofile](man/hcsim-gprofile.md) | Generating CNA profile | a reference genome file |
| [gfasta](man/hcsim-gfasta.md) | Generating clone FASTA file | One or more running directories of previous runs of `gprofile` |
| [gfastq](man/hcsim-gfastq.md) | Generating clone FASTQ file | One or more running directories of previous runs of `gfasta` |
| [align](man/hcsim-align.md) | Aligning clone FASTQ file | One or more running directories of previous runs of `gfastq` |
| [downsam](man/hcsim-downsam.md) | Downsampling clone BAM | One or more running directories of previous runs of `align` |
| [pbam](man/hcsim-pbam.md) | Processing cell BAMs | One or more running directories of previous runs of `downsam` | [Final outputs](doc/chisel-calling.md) |
| [bcbam](man/hcsim-bcbam.md) | Generating barcode BAM file | One or more running directories of previous runs of `bcbam` | [Final outputs](doc/chisel-cloning.md) |

<a name="hcsim-class"></a>
### HCSIM Class

In addition to HCSIM commands, HCSIM also provides the HCSIM class, allowing users to import HCSIM for use in the Python console. The following is an example of using the HCSIM class in the Python console.


```Python
from hcsim import HCSIM

# init HCSIM object with params
hcsim = HCSIM(ref_genome='hg38.fa', clone_no=10, cell_no=100)

# equivalent to the hcsim sim subcommand
hcsim.sim()

# equivalent to the hcsim gprofile subcommand
hcsim.gprofile()

# equivalent to the hcsim gfasta subcommand
hcsim.gfasta()

# equivalent to the hcsim gfastq subcommand
hcsim.gfastq()

# equivalent to the hcsim align subcommand
hcsim.align()

# equivalent to the hcsim downsam subcommand
hcsim.downsam()

# equivalent to the hcsim pbam subcommand
hcsim.pbam()

# equivalent to the hcsim bcbam subcommand
hcsim.bcbam()
```

The parameters supported by the HCSIM class are the same as those supported by the HCSIM command. Below are all the supported parameters. 

```Python
class HCSIM:
    def __init__(self, 
        ref_genome: str, 
        snp_list: str = None, 
        ignore: str = None, 
        outdir: str = './', 
        clone_no: int = 2, 
        cell_no: int = 2, 
        max_tree_depth: int = 4, 
        bin_size: str = '5Mb', 
        snp_ratio: float = 0.001, 
        thread: int = None, 
        heho_ratio: float = 0.67, 
        cna_prob_cutoff: float = 0.8, 
        clone_coverage: float = 30, 
        cell_coverage: float = 0.5, 
        reads_len: int = 150, 
        insertion_size: int = 350, 
        error_rate: float = 0.0, 
        wgd_cna_no: int = 0, 
        wcl_cna_no: int = 0, 
        loh_cna_no: int = 30, 
        goh_cna_no: int = 10, 
        mirror_cna_no: int = 10, 
        barcode_len: int = 12,
        wgsim: str = 'wgsim', 
        samtools: str = 'samtools', 
        bwa: str = 'bwa', 
        picard: str = 'picard.jar',
        bcftools: str = 'bcftools'):
```

<a name="outputs"></a>
### HCSIM Outputs

The command `hcsim sim` runs the entire HCSIM pipeline starting from the required reference genome. During the executiong, the command creates eight folders which contain the temporary and final outputs produced by the each step of HCSIM.


#### `profile` folder

This folder contains all files related to CNA profiles, which are primarily generated by the `hcsim gprofile` command.

1. `cna_profile.csv`: a CSV dataframe containing the haplotype-specific CNAs for each clone. More specifically, the fields are:
   1. `Chromosome`: the chromosome number
   2. `Start`: the start position of bin
   3. `End`: the end position of bin
   4. `normal`: the allele-specific CNAs of normal clone, such as 1|1
   5. `clone1`: the allele-specific CNAs of clone1 clone, such as 3|0
2. `maternal_cna_matrix.csv` and `paternal_cna_matrix.csv`:  the CSV matrix containing the haplotype-specific CNAs of maternal and paternal genome for each clone. These two matrices essentially break down the above CNV profile file, making it more convenient for subsequent plotting and analysis workflows.
3. `changes.csv`: a CSV dataframe containing all the CNA changes. More specifically, the fields are:
	1. `Parent`: the parental clone
	2.	`Child`: the child clone
	3.	`Haplotype`: maternal or paternal genome
	4.	`Type`: type of CNA, all possible types includes `WGD`, `WCL`, `CNL_LOH`, `CNN_LOH`, `CNG_LOH`, `GOH`, `Mirror CNA`, `DUP` and `DEL`
	5.	`Segment`: the genomic bin where CNA occurs
	6.	`Change`: detailed changes in CNA

    The example line following in `changes.csv` means that, transitioning from the normal clone to clone1, a duplication CNA occurs in the maternal genome within the genomic bin chr1:105000001-110000000, with the CNA changing from 1 to 3.

    ```CSV
    normal,clone1,maternal,DUP,chr1:105000001-110000000,1->3
    ```
4. `snp_phases.csv`: a CSV dataframe containing all the SNPs inserted to reference genome. Each line represents a SNP. For the phasing result, 0 means the allele is the same as the allele in the reference genome, while 1 means it is the opposite. Below is an example.

    ```CSV
    chr1,844308,C,A,1|0
    chr1,849963,A,G,1|0
    chr1,850556,A,G,1|0
    chr1,851587,C,A,1|0
    chr1,852718,T,A,1|0
    chr1,855396,G,A,0|1
    chr1,857116,A,C,1|0
    chr1,859092,A,T,1|0
    chr1,859133,T,G,0|1
    chr1,859277,C,T,0|1
    ```
5. `tree.newick`, `tree.pdf` and `tree.json`: the cell-lineage tree in different formats.
6. `barcodes.txt`: the cell list file, with each line representing a cell name.
7. `reference.csv`: a CSV dataframe containing the genomic position information for all bins in the reference genome. This file is used to re-run specific sub-commands.

#### `fasta` folder

This folder contains all maternal and paternal fasta files for all clones.

#### `fastq` folder

This folder contains all maternal and paternal fastq files for all clones.

#### `clone_bams` folder

This folder contains all BAM files for all clones, such as `normal.bam`, `clone1.bam`.

#### `cell_bams` folder

This folder contains all BAM files for all cells in each clone, such as `normal_cell1.bam`, `clone1_cell1.bam`.

#### `barcode_bam` folder

This folder contains barcode bam generating from all cell bams.

#### `log` folder

This folder contains all log files for all HCSIM steps. 

1. `hcsim_log.txt`: the log information generated by HCSIM piepline.
2. `samtools_log.txt`: the log information generated by samtools.
3. `bwa_log.txt`: the log information generated by bwa.
4. `wgsim_log.txt`: the log information generated by wgsim.
5. `picard_log.txt`: the log information generated by picard tools.
6. `barcode_bam_log.txt`: the log information generated by `hcsim bcbam` command.

### `tmp` folder

This folder contains all temporary files for all HCSIM steps.