#!/bin/bash

# This script will generate a plotfile and post that to slack.
# It takes the input arguments time, step number and plotfile name.
# It reads in the slack OAuth token and channel name from the file slack_vars. This file
# should look like
#   export slack_token="xorb-xxx-xxx"
#   export slack_channel="xxxxx"

# slack details
source .slack_vars

webhook_url="https://slack.com/api/files.upload"

export PATH=$MEMBERWORK/ast106/MAESTROeX/Exec/SCIENCE/rotating_star/yt-conda/bin:$PATH
export PYTHONPATH=/sw/xk6/xalt/0.7.5/site:/sw/xk6/xalt/0.7.5/libexec

# process script arguments
t_new=$1
istep=$2
plotfile=$3

plotvar="Hnuc"

# make plot and store plotname. For some reason yt is segfaulting on me at the end,
# so I'm just going to redirect that error
plotname=$(python make_plot.py $plotfile $plotvar 2> tmperr)
rm -f tmperr

text="Output plotfile \`$plotfile\` at timestep \`n=$istep\`, time \`t=$t_new\`. Plotting variable \`$plotvar\`."

echo "Sending figure $plotname to slack"

if [ -z "$plotname" ]; then
    webhook_url="https://slack.com/api/chat.postMessage"
    text="Output plotfile \`$plotfile\` at timestep \`n=$istep\`, time \`t=$t_new\` but there was an error generating plotfile :("
    json="{\"channel\": \"${slack_channel}\", \"text\": \"${text}\"}"
    curl -X POST -H "Authorization: Bearer $slack_token" -H 'Content-type: application/json' --data "${json}" $webhook_url
else
    curl -F file=@${plotname} -F "initial_comment=${text}" -F "channels=${slack_channel}" -H "Authorization: Bearer $slack_token" $webhook_url

    # delete plotfile
    rm $plotname
fi
