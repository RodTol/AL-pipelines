#!/bin/bash

######################## BOT INFO ############################
BOT_TOKEN="$BC_TOKEN_BOT"
CHAT_ID="-4270864261"

# Function to send a message to Telegram
send_message() {
 local message="$1"
 curl -s -X POST "https://api.telegram.org/bot$BOT_TOKEN/sendMessage" \
 -d "chat_id=$CHAT_ID" \
 -d "text=$message"
}

# Function to send a message to Telegram with code block formatting
send_formatted_message() {
    local message=$1
    local formatted_message="\`\`\` $message \`\`\`"

    response=$(curl -s -X POST "https://api.telegram.org/bot$BOT_TOKEN/sendMessage" \
    -d "chat_id=$CHAT_ID" \
    -d "text=$formatted_message" \
    -d "parse_mode=MarkdownV2")
    # Extract message ID from response and save it to a variable
    message_id=$(echo "$response" | jq -r '.result.message_id')    
}

export TEMP_RUN_NAME=$(jq -r '.General.run_name' < $1)
export SLURM_OUTPUT=$(jq -r '.Slurm.output' < $1)
message='cat $SLURM_OUTPUT  | awk '/Indexing completed/,0''
# Check if BC_TOKEN_BOT is defined
if [ -n "$BC_TOKEN_BOT" ]; then
    # Send a "Hello World" message to the Telegram bot
    send_message "Jenkins performed alignment build $TEMP_RUN_NAME succesfully. Here is a quick report:" > /dev/null
    send_formatted_message $message
    #send qualimap report
else
    # Print an error message if BC_TOKEN_BOT is not found
    echo "BC_TOKEN_BOT not found "
fi