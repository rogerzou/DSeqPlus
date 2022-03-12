#!/usr/bin/env bash

cd /mnt/c/Users/rzou4/Downloads/

# Download DISCOVER-Seq (Wienert et al) replicate 1 (negative control)
prefetch SRR8553804 -O /mnt/c/Users/rzou4/Downloads/SRR8553804
cd SRR8553804/SRR8553804
fasterq-dump SRR8553804.sra

# Download DISCOVER-Seq (Wienert et al) replicate 2 (negative control)
prefetch SRR8553806 -O /mnt/c/Users/rzou4/Downloads/SRR8553806
cd SRR8553806/SRR8553806
fasterq-dump SRR8553806.sra
