#!/bin/bash

cat crap.txt | cut -d" " -f1,2,3 > test12
cat crap.txt | cut -d" " -f1,2,4 > test13
cat crap.txt | cut -d" " -f1,3,4 > test23
