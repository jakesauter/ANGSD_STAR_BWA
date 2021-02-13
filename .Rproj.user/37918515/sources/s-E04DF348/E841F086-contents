#! /bin/bash -l

#SBATCH  --partition=angsd_class
#SBATCH  --nodes=1
#SBATCH  --ntasks=1
#SBATCH  --job-name=gierlinski_data_gen
#SBATCH  --time=06:00:00    # HH/MM/SS
#SBATCH  --mem=16G 
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

spack load fastqc
spack load -r py-multiqc

user=$(whoami)
dataset_dir="/athena/angsd/scratch/${user}/gierlinski"


if [ -d "$dataset_dir" ] ; then
  echo -e '\n\n\n${dataset_dir} already exists! Please remove and re-run script!\n\n\n'
  exit
fi

rm -r $dataset_dir ; mkdir -p $dataset_dir ; cd $dataset_dir

if [ $dataset_dir != $PWD ] ; then
  echo "Aborting, not in the correct directory"
  exit
fi


echo -e "\n\n\nDownloading Data to ${dataset_dir}...\n\n\n"

wget "https://ndownloader.figshare.com/files/2194841" --output-document=sample_mapping.tsv

wget "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB5348&result=read_run&fields=run_accession,fastq_ftp&format=tsv&download=true" --output-document=run_accession_and_fastq_ftp_table.txt

mkdir fastqs

for exp in WT SNF2; do
for exp_num in 1 2 3; do
  accession_nums=$(cat sample_mapping.tsv | grep $exp | grep -P "\s$exp_num$" | cut -f1)
  for acc_num in $accession_nums ; do
    ftp=$(egrep $acc_num run_accession_and_fastq_ftp_table.txt | cut -f2)
    wget ftp://$ftp -P "fastqs/${exp}_${exp_num}"
  done
done
done


fastq_files=$(find fastqs/ -name *.fastq.gz | xargs readlink -f)

for file in $fastq_files ; do
  fastqc $file --noextract 
done

multiqc -dd 1 -n gierlinski.multiqc.html fastqs


dirs=$(ls -d  fastqs/* | xargs realpath)

for dir in $dirs ; do
  basename=$(basename $dir)
  filename=$(echo -e "${dir}/${basename}_combined.fastq")
  zcat *.fastq.gz > $filename ; gzip $filename
done
