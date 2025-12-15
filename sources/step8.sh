#!/bin/bash

# Prepare our model database
cat ../data/hmm/*.hmm > ../data/st.hmm
hmmpress ../data/st.hmm


