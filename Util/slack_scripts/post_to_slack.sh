#!/bin/bash

# This script will generate a plotfile and post that to slack.
# It takes the input arguments time, step number and plotfile name.
# It reads in the slack OAuth token and channel name from the file slack_vars. This file
# should look like
#   export slack_token="xorb-xxx-xxx"
#   export slack_channel="xxxxx"

# slack details
source slack_vars

echo $slack_token

webhook_url="https://slack.com/api/files.upload"

# process script arguments
t_new=$1
istep=$2
plotfile=$3

plotvar="magvel"

# make plot and store plotname. For some reason yt is segfaulting on me at the end,
# so I'm just going to redirect that error to null.
plotname=$(python make_plot.py $plotfile $plotvar 2> :)

text="Output plotfile \`$plotfile\` at timestep \`n=$istep\`, time \`t=$t_new\`. Plotting variable \`$plotvar\`."

curl -F file=@${plotname} -F "initial_comment=${text}" -F "channels=${slack_channel}" -H "Authorization: Bearer $slack_token" $webhook_url

 # delete plotfile
rm $plotname
