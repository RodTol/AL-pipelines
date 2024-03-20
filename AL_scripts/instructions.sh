#!/bin/bash

# Color for bash echo
RED="\033[0;31m"
GREEN="\033[0;32m"
CYAN="\033[0;36m"
RESET="\033[0m"

# Input parameters are the json and what node I am on the list 
json_file=$1
my_index=$2

# Read from config.json file (necessary)
host_index=$(jq -r '.Resources.index_host' "$json_file")
node_name=$(jq -r --argjson my_index "$my_index" '.Resources.nodes_list[$my_index]' "$json_file")
node_queue=$(jq -r --argjson my_index "$my_index" '.Resources.nodes_queue[$my_index]' "$json_file")

gpus_settings=$(jq -r --argjson my_index "$my_index" '.Resources.gpus[$my_index]' "$json_file") #debug

# Brief output for checking everything it's correct
echo -e "${RED}I am this node_name: $node_name${RESET}, and for Slurm: $SLURM_NODELIST"
echo $CUDA_VISIBLE_DEVICES
echo -e "${RED}GPUs selected: $gpus_settings${RESET}"
echo -e "${RED}-----------------------${RESET}"
echo "Node queue: $node_queue"
echo -e "${RED}-----------------------${RESET}"
ip a