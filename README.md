<p align="center">
  <img src="docs/assets/logo-area.png" alt="Area logo" width="200"/>
</p>

# AL-pipelines 🧬 
This project contains some Jenkins-pipeline created to perform Alignment on the Orfeo cluster at AreaSciencePark.

## Notes on alignment

Next step in the processing of nanopore data after basecalling is the alignment. I could use `minmap2` but this is already implemented inside `dorado`. So, to launch it, look at the GitHub documentation of the lattest, but more information are available at `minmap2` [documentation](https://github.com/lh3/minimap2?tab=readme-ov-file#general), but also at [samtools](https://www.htslib.org/workflow/fastq.html).  
I need a reference file, and looking at the [dataset](https://42basepairs.com/browse/s3/ont-open-data/cliveome_kit14_2022.05/gdna/flowcells/ONLA29134/20220510_1127_5H_PAM63974_a5e7a202/aligned?file=read_processor_log-2022-05-16_09-11-04.log&preview=contents) page, I decided to use the GRCh38, which is downloadable from [here] (https://www.ncbi.nlm.nih.gov/genome/guide/human/).  

>[!CAUTION]
> I need to have all my **fastq** files merged into one to perform the alignment using all my data. Otherwise I can align individually each file and then merge them all toghether using *samtools*.  
> In order to merge the file you need to something like:
> ```
> cat file*.fastq > bigfile.fastq
> ```

### dorado
Now let's look at the commands with dorado:
```bash 
dorado aligner Basecalled_10G_dataset/hac/pass/ GRCh37_latest_genomic.fna > 10G_aligned.bam
```
Where I have the `.fastq` files inside the `pass` directory and the GRCh37_latest_genomic.fna is the reference genome.  

### minimap2 
If I want a more fine control I need to use the `minimap2` tool. The following command to the alignment for one file:
```
minimap2 -t 8 -a -x map-ont GRCh37_latest_genomic.fna  Basecalled_CliveOME/hac/pass/merged_basecalled.fastq -o tmp/test.sam 
```
Let's analyze it a bit:
- `-t 8`: indicates the number of threads on which the alignment will be executed
- `-a`: to select BAM output instead on paf
- `-x map-ont`: preset of configuration for ONT data
- `-o tmp/test.sam `: output file. If not specified the outpus is stdout

### samtools

To visualize/analyze the output I loaded the samtools (module load samtools) and then runned
```bash
samtools flagstat 10G_aligned.bam # quick resume of some charachteristics
samtools index 10G_aligned.bam #to create index .bam.bai file
samtools tview 10G_aligned.bam #for a rough visualization
```

## Observations
> [!NOTE]  
> `samtools` is the software for managing and using BAM/CRAM files. For this project the most important feature is the possibility to do a **merge**
- I was suggested to use CRAM files instead of normal BAM since they are smaller, but they will need the genome reference to be anlyzed
> [!CAUTION] 
> The alignment requires a lot of memory, so it's important that the script specify this.
