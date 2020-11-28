# BIOF501A Term Project: Bacterial Variant Calling and Annotation Pipeline for Quorum Sensing Genes
## By: Janet Li

![Image of Pseudomonas aeruginosa](Pseudomonas_aeruginosa.png "Pseudomonas aeruginosa")

### Background and Rationale

The purpose of this pipeline is to call and annotate variants in the genome sequence of an AiiA-lactonase treated isolate of *Pseudomonas aeruginosa* strain PAO1 against the reference genome.

*Pseudomonas aeruginosa* is a pathogenic Gram-negative bacterium that causes serious infection in humans. It is often found in the airways of cystic fibrosis patients and is the second most common cause of infections in hospitals. [(Driscoll, Brody & Kollef, 2007)](https://doi.org/10.2165/00003495-200767030-00003). *P. aeruginosa* is also capable of forming biofilms on medical instruments and implants. Antimicrobial resistance is a major public health issue, because infections by resistant bacteria are difficult or even impossible to treat. There are many antimicrobial resistant strains of *P. aeruginosa*, and these mechanisms can arise through a variety biological processes such as modifications to the cell wall or gene regulation through quorum sensing [(Mohanty, Baliyarsignh & Nayak, 2020)](https://doi.org/10.5772/intechopen.88706). Studying the evolution of antimicrobial resistance mechanisms through mutations to genes involved in these biological processes can allow researchers to identify more effective methods of treating *P aeruginosa* infections, as well as new antibacterial targets.

A major aspect of bacterial infection is cell-to-cell communication. Quorum sensing (QS) is the process of multi-cellular gene regulation through the production of signalling molecules. It is used by many bacteria, including *P. aeruginosa*. The production of these signalling molecules increases with cell density, and once cell density reaches a certain threshold, the QS system alters the expression of a subset of genes. QS genes differ between species and strains, but in *P. aeruginosa*, we are interested in the QS genes involved in virulence. Acylated homoserine lactone (AHL) is the most common QS signalling molecule in Gram-negative bacteria [(Venturi, 2006)](https://doi.org/10.1111/j.1574-6976.2005.00012.x). 

AiiA-lactonase is an enzyme that degrades AHL. The isolate used in this pipeline was treated with Aiia-lactonase - more information about the sample can be found [here](https://www.ncbi.nlm.nih.gov//bioproject/667949). The lack of AHL should theoretically reduce the abilities of a bacterial population to communicate among one another, therefore also reducing virulence. This treatment puts a selective pressure on *P. aeruginosa*, essentially selecting for individuals with mutations that allow them to evade the lack of AHL-mediated QS gene regulation. My hypothesis was that there would be at least one mutation in a quorum sensing-related gene. I created this pipeline to call, annotate and parse the genomic variants in the AiiA-lactonase treated isolate.

More specifically, this pipeline:
1. Downloads the reference genome and sequencing reads
2. Indexes the reference genome with `samtools`
3. Aligns the reads to the reference genome with `bwa mem`
4. Calls and filters low-quality variants with `bcftools`
5. Annotates the variants with `snpEff`, and
6. Parses and plots the variant data with a small script, `parse_variants.py`

The pipeline is implemented with `Snakemake`.

Here is a visualization of the workflow:

![Directed acyclic graph of Snakemake workflow](snakemake_workflow_dag.png "Snakemake DAG")

The major package dependencies for this pipeline are:
```bash
python v3.6
snakemake
sra-tools
bwa
samtools v1.9
bcftools
snpEff
graphviz
```

And the Python-specific dependencies are `pandas` and `seaborn` (which include `numpy` and `matplotlib`). Information about installation can be found in the [usage](#usage) section of this README.

---

### Usage
- step by step usage details
- installation (if necessary)
- exact step by step usage with descriptive comments on what action is being performed in each step

1. In the directory where you want to run this pipeline, clone this repository, then navigate into the directory:
```bash
git clone https://github.com/janetxinli/biof501a-term-project.git
cd biof501a-term-project
```

2. Create a `conda` environment with the provided environment file (this will likely only work on a Mac/Unix environment):
```bash
conda env create --file environment.yml
conda activate term_project
```
This will create and activate `conda` environment called `term_project`.

If you're using a different operating system, you most likely won't able to create the environment from the file directly. In that case, check out the [supplementary steps](#supplementary-steps). `conda` environments don't translate well across operating systems, and you'll have to set up the environment manually (but don't worry, it's pretty simple!).

3. Run the pipeline with `Snakemake`. You'll need to set the number of cores that `Snakemake` uses. For example, if you want to run the pipeline with two cores, specify `--cores 2`.

```bash
snakemake --cores 2
```

I like to use `nohup` so I can run other processes at the same time or log out without stopping the process:

```bash
nohup snakemake --cores 2 &> snakemake.out &
```

This will print all stdout and stderr messages to a file called `snakemake.out`.

4. Compare your outputs with the files in the `expected_outputs` directory. More details in the [outputs](#outputs) section.

#### Supplementary steps
To create and activate a `conda` environment called `term_project` with all of the requirements manually:

```bash
conda env create --name term_project python=3.6
conda activate term_project
```

Some of the other dependencies are not compatible with higher versions of Python, so make sure to specify the version. Next, to install the dependencies:

```bash
conda install -c bioconda snakemake sra-tools bwa samtools=1.9 bcftools snpEff openssl=1.0
conda install -c anaconda graphviz seaborn
```

Your environment is now set up with all of the dependencies required for this pipeline.

---

### Input

### Ouput

