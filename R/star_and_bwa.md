---
title: "Module 5 Exercise: Aligning FASTQ reads and interpreting the output"
author: "Jake Sauter"
date: "2/12/2021"
output: 
  html_document: 
    toc: true
    toc_float: true
    keep_md: true
---



## **Running `BWA` and `STAR`**

To begin, I will present scripts that will:

1. Retrieve and organize the **Gierlinski** dataset on the user's system.
2. Execute `BWA` on all data collected from the `WT_1` sample from this dataset..<br>
3. Execute `STAR` on the same sample

These scripts will ensure that the resultant files will:

* Be in the `.bam` compressed alignment file format.
* The respective `.bam` file will be sorted.
* The respective `.bam` file will be indexed and have a corresponding `.bai` index file.

### **Retrieving the Gierlinkski Dataset**

The Gierlinkski Dataset can be downloaded via [this slurm-submittable script](https://raw.githubusercontent.com/jakesauter/ANGSD_STAR_BWA/main/scripts/slurm_gierlinski_download_and_qc_script.sh).


### **Running BWA on WT_1**

We will choose **Wild-Type biological replicate 1** (`WT_1`) as
our sample to run as input to `BWA`. 

In order to run `BWA`, we first need the **.fasta** file that
contains the reference genome we are aligning to. For this
exercise this has been obtained for us and we can copy 
from `/home/luce/angsd/referenceGenomes/sacCer3.fa`

The following script when executed will retrieve
this reference into the project data-lake architecture. It will then build our `BWA` reference using `bwa index`, followed
by performing the alignment process with `bwa mem`. This script has also been setup as so it is slurm-submittable. 


```bash
#! /bin/bash -l

#SBATCH  --partition=angsd_class
#SBATCH  --nodes=1
#SBATCH  --ntasks=1
#SBATCH  --job-name=gierlinski_bwa
#SBATCH  --time=06:00:00    # HH/MM/SS
#SBATCH  --mem=16G 
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err


spack load bwa@0.7.15% gcc@6.3.0
spack load samtools@1.9% gcc@6.3.0

user=$(whoami)
gierlinski_dataset="/athena/angsd/scratch/${user}/gierlinski"

fastq_prefix="WT_1_combined"
fastq_file=$(find ${gierlinski_dataset}/fastqs/ -name ${fastq_prefix}.fastq.gz)

if [ ! -f "$fastq_file" ] ; then
  echo -e "\n\nA fastq file with the given prefix ${fastq_prefix} could not be"
  echo -e "found in the datset directory: ${gierlinski_dataset}! Aborting!\n\n"
  exit
else
  echo -e "\n\nFastq file ${fastq_file} found! Now beginging alignment ...\n\n"
fi

sam_alignment_file="${bwa_alignment_dir}/${fastq_prefix}.bwa.sam"
bwa_index_file="${gierlinski_dataset}/sacCer3_BWAindex/sacCer3"
mkdir -p `dirname $bwa_index_file`

reference_file="${gierlinski_dataset}/sacCer3_ref/sacCer3.fa"
bwa_alignment_dir="${gierlinski_dataset}/bwa_alignments"
samtools_qc_dir="${gierlinski_dataset}/samtools_qc"
mkdir -p `dirname $reference_file`
mkdir -p $bwa_alignment_dir
mkdir -p $samtools_qc_dir
cp /home/luce/angsd/referenceGenomes/sacCer3.fa $reference_file 

bwa index -p $bwa_index_file $reference_file
bwa mem $bwa_index_file $fastq_file > $sam_alignment_file

bam_output_file="${bwa_alignment_dir}/${fastq_prefix}.bwa.bam"
bwa_sorted_bam="${bwa_alignment_dir}/${fastq_prefix}.bwa.sorted.bam"
samtools view -b $sam_alignment_file -o $bam_output_file
samtools sort $bam_output_file -o $bwa_sorted_bam
samtools index $bwa_sorted_bam 

rm $sam_alignment_file

samtools stats $bwa_sorted_bam  > $samtools_qc_dir/${fastq_prefix}.bwa.bam.stats
samtools flagstat $bwa_sorted_bam > $samtools_qc_dir/${fastq_prefix}.bwa.bam.flagstats
```


### **Running `STAR` on WT_1**


In order to run `STAR` on this same sample, along side our previous
**.fasta** reference file we will also need a **.gtf** feature
file. This file can be copied from another location already on 
the system: `/home/luce/angsd/referenceGenomes/sacCer3.sgd.gtf`.

The following script when executed will retrieve this reference into the project data-lake architecture. It will then build our `STAR` reference using `STAR --runMode  genomeGenerate`, followed by performing the alignment process with `STAR --runMode  alignReads`. This script has also been setup as so it is slurm-submittable.


```bash
#! /bin/bash -l

#SBATCH  --partition=angsd_class
#SBATCH  --nodes=1
#SBATCH  --ntasks=1
#SBATCH  --job-name=gierlinski_star
#SBATCH  --time=06:00:00    # HH/MM/SS
#SBATCH  --mem=16G 
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

spack load star@2.7.0e
spack load samtools@1.9% gcc@6.3.0

user=$(whoami)
gierlinski_dataset="/athena/angsd/scratch/${user}/gierlinski"

fastq_prefix="WT_1_combined"
fastq_file=$(find ${gierlinski_dataset}/fastqs/ -name ${fastq_prefix}.fastq.gz)

if [ ! -f "$fastq_file" ] ; then
  echo -e "\n\nA fastq file with the given prefix ${fastq_prefix} could not be"
  echo -e "found in the datset directory: ${gierlinski_dataset}! Aborting!\n\n"
  exit
else
  echo -e "\n\nFastq file ${fastq_file} found! Now beginging alignment ...\n\n"
fi

star_index_dir="${gierlinski_dataset}/sacCer3_STARindex"
reference_file="${gierlinski_dataset}/sacCer3_ref/sacCer3.fa"
gtf_reference_file="${gierlinski_dataset}/sacCer3_ref/sacCer3.sgd.gtf"
alignment_file="${gierlinski_dataset}/star_alignments/${fastq_prefix}.star.bam"
samtools_qc_dir="${gierlinski_dataset}/samtools_qc"

mkdir -p $star_index_dir
mkdir -p $samtools_qc_dir
mkdir -p `dirname $reference_file`
mkdir -p `dirname $alignment_file`

cp /home/luce/angsd/referenceGenomes/sacCer3.sgd.gtf $gtf_reference_file 

cd $gierlinski_dataset

STAR --runMode  genomeGenerate \
     --runThreadN 6 \
     --genomeDir  $star_index_dir \
     --genomeFastaFiles $reference_file \
     --sjdbGTFfile  $gtf_reference_file \
     --sjdbOverhang  99
  
STAR --runMode  alignReads \
  --runThreadN 6 \
  --genomeDir $star_index_dir \
  --readFilesIn $fastq_file \
  --readFilesCommand zcat \
  --outFileNamePrefix $alignment_file \
  --outSAMtype BAM SortedByCoordinate
  
  
alignment_path=$(dirname $alignment_file)
star_sorted_bam=$(ls ${alignment_path}/*.bam)

samtools index $star_sorted_bam

samtools stats $star_sorted_bam  > $samtools_qc_dir/${fastq_prefix}.star.bam.stats
samtools flagstat $star_sorted_bam > $samtools_qc_dir/${fastq_prefix}.star.bam.flagstats
```


### **UCSD Integrated Genomics Viewer (IGV)**

As seen below, the UCSD [Integrated Genomics Viewer (IGV)](https://igv.org/) can be used to view our **sorted**, **indexed** `.bam` files. The image below shows only a small portion of the output from the previous `BWA` alignment, though illustrates how each read alignment can be interacted with to see all the relevant techincal details about it.

![](/home/x1/Documents/Weill_Cornell/ANGSD/homework/homework_4/imgs/IGV_verification_bwa.png)

**Multiple Tracks**

Due to the ability of **IGV** to allow for multiple "Tracks" to be 
uploaded at once, we can use this functionality to check that the
alignments we have performed with `BWA` and `STAR` respectively
are in general agreement as a first rough QC measure and sanity check.

![](/home/x1/Documents/Weill_Cornell/ANGSD/homework/homework_4/imgs/bwa_star_concordance.png)


## **Subsetting Chromosomal Reads**

In order to subset reads only aligned to particular 
regions of the reference, we can use the `samtools` command
line utility tool. In the usage example that we see below, 
the region of interest can be specified as the last 
positional argument to the function.


```bash
samtools view [options] in.sam|in.bam|in.cram [region...]
```


```bash
spack load samtools@1.9% gcc@6.3.0

bwa_sorted_bam="/athena/angsd/scratch/jns4001/gierlinski/bwa_alignments/WT_1_combined.bwa.sorted.bam"
star_sorted_bam="/athena/angsd/scratch/jns4001/gierlinski/star_alignments/WT_1_combined.star.bamAligned.sortedByCoord.out.bam"

bwa_filtered_output_file="/athena/angsd/scratch/jns4001/gierlinski/bwa_alignments/WT_1_combined.bwa.chrI.sam"
star_filtered_output_file="/athena/angsd/scratch/jns4001/gierlinski/star_alignments/WT_1_combined.star.bamAligned.sortedByCoord.chrI.sam"

samtools view -h $bwa_sorted_bam chrI > $bwa_filtered_output_file
samtools view -h $star_sorted_bam chrI > $star_filtered_output_file
```

To confirm that this filtering has actually worked, we can
cut the 3rd field (being the chromosomal region that
the query sequence was mapped too) and check that this
always falls on `chrI`.
    

```bash
samtools view $bwa_filtered_output_file | cut -f3 | uniq
```


```bash
chrI
```


```bash
samtools view $star_filtered_output_file | cut -f3 | uniq
```


```bash
chrI
```


## **Comparing `BWA` and `STAR` Output**

When determining the differences between the output
alignment `.bam` files produced from popular alignment
tools such as `BWA` and `STAR`, some good places 
to start could include the **header** as well as 
the **optional SAM fields** that the programs
choose include. 

**SAM Headers**
    
In order to view the differences in headers to begin with, 
we can use the `samtools view` command, passing the `-H` flag
to only include the SAM header.
  

```bash
spack load samtools@1.9% gcc@6.3.0

bwa_sorted_bam="/athena/angsd/scratch/jns4001/gierlinski/bwa_alignments/WT_1_combined.bwa.sorted.bam"
star_sorted_bam="/athena/angsd/scratch/jns4001/gierlinski/star_alignments/WT_1_combined.star.bamAligned.sortedByCoord.out.bam"

samtools view -H $bwa_sorted_bam
samtools view -H $star_sorted_bam
```

In order to more easily compare these outputs, I have placed
them side with the `BWA` header on the **left** and the 
`STAR` header on the **right**.
  
![](/home/x1/Documents/Weill_Cornell/ANGSD/homework/homework_4/imgs/bwa_vs_star_header.png)

We can see that the headers do not contain many differences, 
though if parsing these files manually it is important to keep
in mind that the headers are not **exactly** the same. 

**SAM Optional Fields**

We will now compare how `BWA` and `STAR` both use SAM
optional fields to provide more information about each
alignment the program has picked.
  

```bash
samtools view $bwa_sorted_bam | head -n 1
```


```bash
ERR458496.275402	0	chrI	1800	0	51M	*	0	0TTTGTCTCTAGTTTGCGATAGTGTAGATACCGTCCTTGGATAGAGCACTGG	?+1ADFFFGDHHHJJIIGGIGHGCHGGIIEIIIIJJJBHB@DHHIIIJJGI	NM:i:0	MD:Z:51	AS:i:51	XS:i:51
```



```bash
samtools view $star_sorted_bam | head -n 1
```


```bash
ERR458493.552967	16	chrI	140	255	12M61232N37M2S	*	00	CCACTCGTTCACCAGGGCCGGCGGGCTGATCACTTTATCGTGCATCTTGGC	BB?HHJJIGHHJIGIIJJIJGIJIJJIIIGHBJJJJJJHHHHFFDDDA1+B	NH:i:1	HI:i:1	AS:i:41	nM:i:2
```

In order to dissect the differences in these optional field outputs, 
we must refer to the the manual of each program. 

**Reference Manuals**: 

* [BWA Manual](http://bio-bwa.sourceforge.net/bwa.shtml)
*  [STAR Manual](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)

In these documents, we can find the excerpts below. Again, information regarding `BWA` is on the **left** and information regarding `STAR`
is on the **right**. 


![](/home/x1/Documents/Weill_Cornell/ANGSD/homework/homework_4/imgs/bwa_vs_star_attributes.png)

From these excerpts, we can now dissect the differences in 
optional SAM field output. From previously printing an 
entry from `BWA`, we can see that the program listed the following output during out analysis: 

* **NM** -- Edit distance 
* **MD** -- Mismatching positions/bases
* **AS** -- Alignment Score
* **XS** -- Suboptimal Alignment Score ; X indicates that the
                tag is specific to` BWA`.

And `STAR` Listed: 

* **NH** -- Number of reported alignments that contain the query in the current record 
* **HI** -- Query hit index 
* **AS** -- Alignment score generated by aligner 
* **nM** -- Number of mismatches per (paired) alignment, not 
    to be confused with NM, which is the number of mismatches
    in each mate.


## **`BamQC`**

In order to determine more closely the quality of our mappings, 
we can use the `BamQC` tool it via `/softlib/apps/EL7/BamQC/bin/bamqc` after logging on to an `SCU` compute node. As is the first step
in learning any new tool or command, we will ask for `--help`. 


```bash
/softlib/apps/EL7/BamQC/bin/bamqc --help 
```

As apart of this output, we can read the **DESCRIPTION** section
to understand more about the purpose and logic of this tool. 


```bash
DESCRIPTION

    BamQC reads a set of mapped BAM files and produces from each one a quality control report consisting of a number of different modules, each one of which will help to identify a different potential type of problem in your data.
    
    If no files to process are specified on the command line then the program will start as an interactive graphical application.  If files are provided on the command line then the program will run with no user interaction required.  In this mode it is suitable for inclusion into a standardised analysis pipeline.

```


We can see that this tool could be helpful for understanding
the quality of our BAM files that we have just produced
from `BWA` and `STAR`. In order to use `BamQC`, we must pass
appropriate parameters to get reasonble results. 
Garbage in, garbage out. 

Looking more closely into the result of the `--help` command, 
we find that we can list all reference genomes already 
available to use through the Babraham Server. Since we are
looking for the genome of Saccharomyces cerevisiae, we 
will filter the output for  containing `"cerevisiae"`. 


```bash
/softlib/apps/EL7/BamQC/bin/bamqc --available 
```

We will find the following line in the output. 


```bash
Saccharomyces cerevisiae [ EF4 | R64-1-1 | SGD1.01 | SGD1 | SGD1_new ]
```

This line indicates to us that we can query the `--species` parameter
with `Saccharomyces cerevisiae`, with any of the listed `--assembly`
options. Upon inspection of the [website hosting the genomes](https://www.bioinformatics.babraham.ac.uk/seqmonk/genomes/Saccharomyces%20cerevisiae/), **EF4** appears to be the newest. 

![](/home/x1/Documents/Weill_Cornell/ANGSD/homework/homework_4/imgs/yeast_index_seqmonk.png)

In order to download this genome for use, we can use the following command. 


```bash
/softlib/apps/EL7/BamQC/bin/bamqc -s "Saccharomyces cerevisiae" -a "EF4"
```

We can then see this saved genome with the `--saved` command.


```bash
/softlib/apps/EL7/BamQC/bin/bamqc --saved
```

We can see where our reference has downloaded to. 


```bash
Downloaded genomes:
/home/jns4001/BamQC_files/genomes/Saccharomyces cerevisiae
```


We will now call `BamQC` with our sorted `BAM` files that we have previously generated through both `BWA MEM` and `STAR`.

**Note: Although we have seen how to download reference 
genomes and assemblies through this tool, the use of any 
I have tried results in empty feature coverage graphs leading me to believe incompatibility issue with specific versions of genomes
and features exists. Therefore I will opt to not providing
and of these parameters explicitly as this does not cause
the formation of empty graphs.**


```bash
bwa_sorted_bam="/athena/angsd/scratch/jns4001/gierlinski/bwa_alignments/WT_1_combined.bwa.sorted.bam"
star_sorted_bam="/athena/angsd/scratch/jns4001/gierlinski/star_alignments/WT_1_combined.star.bamAligned.sortedByCoord.out.bam"

/softlib/apps/EL7/BamQC/bin/bamqc $bwa_sorted_bam
/softlib/apps/EL7/BamQC/bin/bamqc $star_sorted_bam
```

The result of running this command is a `BamQC` generated `.html` file in the same directory as the input files. In our case, these files are 

* `BWA`: WT_1_combined.bwa.sorted_bamqc.html
* `STAR`: WT_1_combined.star.bamAligned.sortedByCoord.out_bamqc.html

Three differences observed between STAR and BWA: 

**1.)** **Differences in Basic Statistics**

As seen from the table below, `BAM` and `STAR` appeared
to have mapped a different number of **Total sequences**, 
as well as the amount of sequences that were assigned
a **primary alignment** is different. We also see from 
the data below that `STAR` left no read unmapped.


Metric | BWA | STAR 
------ | --- | ----
Total sequences | 7014727 | 7585553
Percent primary alignments | 99.998 | 88.718
Percent sequences unmapped | 2.967 |  0.000
Percent sequences discarded for variant call detection |  2.967 | 0.000

**2.)** **`STAR`'s Increased Soft Clipping** 

Another difference seen from the plots below is the fact that
`STAR` has seemed to allow for more **soft-clipping** on the 
5' end of the reads of **length 2** compared to the soft
clips generated in the alignments of `BWA`. 

**Left: `BWA`, Right: `STAR`**

![](/home/x1/Documents/Weill_Cornell/ANGSD/homework/homework_4/imgs/soft_clip_len_dist_comparison.png)

**3.)** **Mapping Quality Implementations** 

Finally, the third difference that I have spotted is
the value and distribution of **Mapping Quality** values
generated by both tools. Mapping quality values for 
confidently mapped reads default to **60** for `BWA`, 
however the same confident mapped reads will be 
assigned a **255** when using `STAR`.

**Left: `BWA`, Right: `STAR`**

![](/home/x1/Documents/Weill_Cornell/ANGSD/homework/homework_4/imgs/mapq_dist_comparison.png)


## **Alignment Score vs Mapping Quality**

Useful articles: 

* [**SEQanswers forum question**](http://seqanswers.com/forums/showthread.php?t=66634)
* [**Biostars forum question**](https://www.biostars.org/p/179457/)

The difference between **alignment score** and **mapping quality** in SAM/BAM files is an important one. 
As indicated by the **Biostars forum answers**, *alignment score* (*AS*) tells you how similar the read is to the reference, though this should not be taken at face value. *Mapping quality* (*MAPQ*) is a more dependable metric that tells you how confident you can be that the read comes from the reported position.

A situation in which this difference is important: You can have high *AS* and low *MAPQ* if the read aligns perfectly at multiple positions, and you can have low *AS* and high *MAPQ* if the read aligns with mismatches but still the reported position is still much more probable than any other.

Furthermore, interpretation of the **mapping quality** field is different between `BWA` and `STAR`!

We can learn more about how `STAR` implements *mapping quality* 
from the [STAR manual](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf).

> The mapping quality MAPQ (column 5) is 255 for uniquely mapping reads, and int(-10*log10(1-1/Nmap)) for multi-mapping reads.  This scheme is same as the one used by TopHat and is com-patible  with  Cufflinks.    The  default  MAPQ=255  for  the  unique  mappers  maybe  changed  with--outSAMmapqUniqueparameter (integer 0 to 255) to ensure compatibility with downstream toolssuch as GATK.

The [BWA Manual](http://bio-bwa.sourceforge.net/bwa.shtml) does not
appear to allude to the specifics of *Mapping Quality* calculation, 
however we can spot from the **Mapping Quality Distribution** images
from above that `BWA` must assign unique mapping reads a value of
**60**, just like how `STAR` assigns uniquely mapping reads a value
of **255**. We can also assume that reads with low mapping qualities
(around 0) may have equally good *alignment scores* over many 
locations in the referece.


## **Multi-Mapping Reads vs Split-Reads**

In this section the difference between a **mutli-mapping read** and
**split read** will first be defined, then we will delve deeper 
and show how the same `STAR` multi-mapping read was handled in `BWA`.

A **multi-mapping** read is a read that has mapped to multiple locations in the reference, while a **split read** is a read that has been mapped to two different portions of the reference, with a large *split* in the 
middle of the read alignment. 

**Useful Resources:**

* [Forum answer on split-reads](https://bioinformatics.stackexchange.com/questions/5233/what-are-split-reads-and-intron-clusters#5234)


**Finding a split read in STAR**

We can identify a split read in star by finding a read that
has a large `"N"` placement to indicate an intron has been 
spliced out. 


```bash
ERR458497.799793	16	chrIII	177875	255	32M307N19M	*	0	0	AACTTGGGAATTGTCACGAGCTTGAACAACGTTAGACATGGCGGGTTCTTG	FCGGIJJHEIJJJJJIJJJJIJJJIGIJIJIJJJJIJIHHGHGFFEFFC@C	NH:i:1	HI:i:1	AS:i:51	nM:i:0
```

In this read on Chromosome 3, it has matched 32 bases, followed by 307 introns spliced, with the remaining 19 matches occurring after.


```bash
samtools view ERR458497.bwa.sorted.bam | egrep AACTTGGGAATTGTCACGAGCTTGAACAACGTTAGACATGGCGGGTTCTTG
ERR458497.51047	16	chrIII	177875	60	32M19S	*	0	0	AACTTGGGAATTGTCACGAGCTTGAACAACGTTAGACATGGCGGGTTCTTG	
```

It appears that `BWA` has identified the same 32 matching nucleotides, however due to `BWA` not being splice aware, it calls the remainder of the bases in the sequence alignment as a **soft clip (S)**.

## **Removing Unampped Reads From BWA**

In order to remove unmapped reads from `BWA` output, the `u` flag can be used to filter out sequences that "the query sequence itself is unmapped". This flag is mapped to integer `4`, and thus can be used in `samtools view`


```bash
samtools view -F 4 -b -o ERR458497.bwa.sorted.mapped.bam ERR458497.bwa.sorted.bam
```

Now if we filter to "only include reads with all  of the FLAGs in INT present " with 
the `-f` option, our new file with supposedly only mapped reads should 
return no reads at all, and the output of the following command
should be blank!


```bash
samtools view -f 4 ERR458497.bwa.sorted.mapped.bam
```
















