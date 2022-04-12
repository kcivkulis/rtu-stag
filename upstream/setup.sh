#!/bin/bash
#
# clear out the old installation folders
cd ../..
rm -rf stag-mwc

# first thing's first clone stag if it's not here already
git clone https://github.com/ctmrbio/stag-mwc

# Clone Bracken
git clone https://github.com/jenniferlu717/Bracken
