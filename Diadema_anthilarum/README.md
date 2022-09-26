# Assembly of *Diadema antillarum* mtDNA

## List of files:

Genome assembly: **diadema_anthilarium_mtDNA.fasta**

Genome annotation with MITOS: **diadema_anthilarium_mtDNA.gff**

Genome visualization: **genome.jpeg**

Genome alignment: **diadema_anthilarium_mtDNA.sorted.bam**

Coverage: **diadema_anthilarium_mtDNA.cov**

![Fig1_revised](https://user-images.githubusercontent.com/142793/192270699-20a0e533-365b-4ec9-9441-3fa5efefa9b0.jpg)

The raw fastq data that supports the assembly were extracted with [Cookiecutter](https://github.com/ad3002/Cookiecutter): **mtDNA_1.fastq.gz** and **mtDNA_2.fastq.gz**

## Reproducibility instruction:

You can download the extracted mtDNA-positive reads and assemble them using any prokaryotic assembler such as SPAdes.

Requirements:

- x64 Linux system
- git

Step 1.

Download cookiecutter2 toolset:

```bash
git clone git@github.com:aglabx/Tools.git
cd Tools
```

Step 2. Trim data with v2trim

```bash
./V5_trim.exe <raw_reads_1.fastq> <raw_reads_2.fastq> <trimmed> <threads> 0 fastq illumina_ext.data
```

(illumina_ext.data file you can find in Tools folder)

Step 3. Filter out mtDNA positive reads:

```bash
./Cookiecutter2.exe <trimmed_1.fastq> <trimmed_1.fastq> <mtDNA_reads> <threads> fastq 23 extract 5 kmers_from_the_reference.kmers
```

**kmers_from_the_reference.kmers** file contains one 23-mer per line data computed from the any close reference fasta file (we used diadema.kmers from this repository).

Step 4. Assembly extracter

```bash
spades.py -1 mtDNA_1.fastq.gz -2 mtDNA_2.fastq.gz -o assembly
```

Resulting assembly with partially assembled NUMTs in the file: **scaffolds.spades.fasta**

After that, we manually rotated it to the usual beginning and do the reverse complement (the final assembly available in the file: **diadema_anthilarium_mtDNA.fasta**).



