#!/bin/bash
cat << "EOF"

 ____                        ___                                                     
/\  _`\   __                /\_ \    __                                              
\ \ \L\ \/\_\  _____      __\//\ \  /\_\    ___      __                              
 \ \ ,__/\/\ \/\ '__`\  /'__`\\ \ \ \/\ \ /' _ `\  /'__`\                            
  \ \ \/  \ \ \ \ \L\ \/\  __/ \_\ \_\ \ \/\ \/\ \/\  __/                            
   \ \_\   \ \_\ \ ,__/\ \____\/\____\\ \_\ \_\ \_\ \____\                           
    \/_/    \/_/\ \ \/  \/____/\/____/ \/_/\/_/\/_/\/____/                           
                 \ \_\                                                               
                  \/_/                                                               
 __                                        ____                        __            
/\ \                                      /\  _`\                     /\ \           
\ \ \        ___     ___      __          \ \ \L\ \     __     __     \_\ \    ____  
 \ \ \  __  / __`\ /' _ `\  /'_ `\  _______\ \ ,  /   /'__`\ /'__`\   /'_` \  /',__\ 
  \ \ \L\ \/\ \L\ \/\ \/\ \/\ \L\ \/\______\\ \ \\ \ /\  __//\ \L\.\_/\ \L\ \/\__, `\
   \ \____/\ \____/\ \_\ \_\ \____ \/______/ \ \_\ \_\ \____\ \__/.\_\ \___,_\/\____/
    \/___/  \/___/  \/_/\/_/\/___L\ \         \/_/\/ /\/____/\/__/\/_/\/__,_ /\/___/ 
                              /\____/                                                
                              \_/__/                                                 

EOF



module load samtools

# Color for bash echo
RED="\033[0;31m"
GREEN="\033[0;32m"
CYAN="\033[0;36m"
RESET="\033[0m"

# Input parameters are the json and what node I am on the list 
json_file=$1
my_index=$2

# Read from config.json file 
threads=$(jq -r --argjson my_index "$my_index" '.Alignment.threads[$my_index]' "$json_file")
input_file=$(jq -r '.Alignment.input_file' "$json_file") 
output_file=$(jq -r '.Alignment.output_file' "$json_file") 
reference_genome=$(jq -r '.Alignment.reference_genome' "$json_file")

node_queue=$(jq -r --argjson my_index "$my_index" '.Resources.nodes_queue[$my_index]' "$json_file")
node_name=$(jq -r --argjson my_index "$my_index" '.Resources.nodes_list[$my_index]' "$json_file")

# Brief output for checking everything it's correct
# Update config.json file
if  [ "$node_name" == "" ]; then
  echo -e "${RED}|||Update node_name from $node_name to $SLURM_NODELIST|||${RESET}"
  node_name=$SLURM_NODELIST
fi
echo -e "${RED}I am this node_name: $node_name${RESET}, and for Slurm: $SLURM_NODELIST"
echo -e "${RED}-----------------------${RESET}"
echo "Threads to be used: $threads"
echo "Reference genome: $reference_genome"
echo "Node queue: $node_queue"
echo "Input : $input_file"
echo "Output : $output_file"
echo -e "${RED}-----------------------${RESET}"

# Load virtualenv for Python:
# None for now

#Maye trimming and filtering

# Start the alignment (fastq to sam )
minimap2 -t $threads -ax map-ont $reference_genome  $input_file -o $output_file
echo "Alignment completed"

#Create unsorted bam file
bam_file_not_sorted="${output_file%.sam}_not_sorted.bam"
echo "Converting to $bam_file_not_sorted"

samtools view -Sb $output_file > $bam_file_not_sorted
echo "Conversion ended"

#Sort the bam file
bam_file="${output_file%.sam}.bam"
samtools sort $bam_file_not_sorted -o $bam_file
echo "Sorting completed"

#index the bam file
samtools index $bam_file
echo "Indexing completed"

#compute alignment statistic
samtools flagstat $bam_file

#QC report
#module load R
#module load java
#unset DISPLAY
#qualimap bamqc -bam batch_293436308571279569530322515273083434315.bam -outdir qualimap

#Variant calling/Mark duplicates
