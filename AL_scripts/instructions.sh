#!/bin/bash

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
logs_dir=$(jq -r '.Alignment.logs_dir' "$json_file") 
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
echo "Logs Directory: $logs_dir"
echo "Input Directory: $input_file"
echo "Output Directory: $output_file"
echo -e "${RED}-----------------------${RESET}"

# Load virtualenv for Python:
# None for now

# Start the alignment
minimap2 -t $threads -a -x map-ont $reference_genome  $input_file -o $output_file