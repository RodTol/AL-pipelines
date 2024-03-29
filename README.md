<p align="center">
  <img src="docs/assets/logo-area.png" alt="Area logo" width="200"/>
</p>

# AL-pipelines 🧬 
This project contains some Jenkins-pipeline created to perform Alignment on the Orfeo cluster at AreaSciencePark.

## Notes on alignment

Next step in the processing of nanopore data after basecalling is the alignment. I could use `minmap2` but this is already implemented inside `dorado`. So, to launch it, look at the GitHub documentation of the lattest, but more information are available at `minmap2` documentation.  
I need a reference file, and looking at the [dataset](https://42basepairs.com/browse/s3/ont-open-data/cliveome_kit14_2022.05/gdna/flowcells/ONLA29134/20220510_1127_5H_PAM63974_a5e7a202/aligned?file=read_processor_log-2022-05-16_09-11-04.log&preview=contents) page, I decided to use the GRCh38, which is downloadable from [here] (https://www.ncbi.nlm.nih.gov/genome/guide/human/).  
Now let's look at the commands:
```bash 
dorado aligner Basecalled_10G_dataset/hac/pass/ GRCh37_latest_genomic.fna > 10G_aligned.bam
```
Where I have the `.fastq` files inside the `pass` directory and the GRCh37_latest_genomic.fna is the reference genome.  
To visualize/analyze the output I loaded the samtools (module load samtools) and then runned
```bash
samtools index 10G_aligned.bam #to create index .bam.bai file
samtools tview 10G_aligned.bam #for a rough visualization
```

## Observations
- `samtools` is the software for managing and using BAM/CRAM files
- I was suggested to use CRAM files instead of normal BAM since they are faster