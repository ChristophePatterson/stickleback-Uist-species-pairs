#!/bin/bash

source ~/.bashrc
conda activate LiftOver-env

wkdir=/gpfs01/home/mbzcp2/data/sticklebacks/genomes
## Files to workwith
ls -1 $wkdir/Prior_gasAcu-results/*v4.bed

## EcoPeaks
liftOver $wkdir/Prior_gasAcu-results/Roberts-et-al-2021-Sensitive-EcoPeaks-v4.bed $wkdir/v4_to_v5.chain.txt $wkdir/Prior_gasAcu-results/Roberts-et-al-2021-Sensitive-EcoPeaks-v5.bed $wkdir/Prior_gasAcu-results/Roberts-et-al-2021-Sensitive-EcoPeaks-v5_unLifted.bed
liftOver $wkdir/Prior_gasAcu-results/Roberts-et-al-2021-Specific-EcoPeaks-v4.bed $wkdir/v4_to_v5.chain.txt $wkdir/Prior_gasAcu-results/Roberts-et-al-2021-Specific-EcoPeaks-v5.bed  $wkdir/Prior_gasAcu-results/Roberts-et-al-2021-Specific-EcoPeaks-v5_unLifted.bed

## NE EcoPeaks
liftOver $wkdir/Prior_gasAcu-results/Roberts-et-al-2021-Sensitive-NE-EcoPeaks-v4.bed $wkdir/v4_to_v5.chain.txt $wkdir/Prior_gasAcu-results/Roberts-et-al-2021-Sensitive-NE-EcoPeaks-v5.bed $wkdir/Prior_gasAcu-results/Roberts-et-al-2021-Sensitive-NE-EcoPeaks-v5_unLifted.bed
liftOver $wkdir/Prior_gasAcu-results/Roberts-et-al-2021-Specific-NE-EcoPeaks-v4.bed $wkdir/v4_to_v5.chain.txt $wkdir/Prior_gasAcu-results/Roberts-et-al-2021-Specific-NE-EcoPeaks-v5.bed $wkdir/Prior_gasAcu-results/Roberts-et-al-2021-Specific-NE-EcoPeaks-v5_unLifted.bed

## TempoPeaks
liftOver $wkdir/Prior_gasAcu-results/Roberts-et-al-2021-Sensitive-TempoPeaks-v4.bed $wkdir/v4_to_v5.chain.txt $wkdir/Prior_gasAcu-results/Roberts-et-al-2021-Sensitive-TempoPeaks-v5.bed $wkdir/Prior_gasAcu-results/Roberts-et-al-2021-Sensitive-TempoPeaks-v5_unLifted.bed
liftOver $wkdir/Prior_gasAcu-results/Roberts-et-al-2021-Specific-TempoPeaks-v4.bed $wkdir/v4_to_v5.chain.txt $wkdir/Prior_gasAcu-results/Roberts-et-al-2021-Specific-TempoPeaks-v5.bed $wkdir/Prior_gasAcu-results/Roberts-et-al-2021-Specific-TempoPeaks-v5_unLifted.bed

## Jones et al 
liftOver $wkdir/Prior_gasAcu-results/Jones-et-al-2012-CSS-02-v4.bed $wkdir/v4_to_v5.chain.txt $wkdir/Prior_gasAcu-results/Jones-et-al-2012-CSS-02-v5.bed $wkdir/Prior_gasAcu-results/Jones-et-al-2012-CSS-02-v5_unLifted.bed
liftOver $wkdir/Prior_gasAcu-results/Jones-et-al-2021-CSS-05-v4.bed $wkdir/v4_to_v5.chain.txt $wkdir/Prior_gasAcu-results/Jones-et-al-2021-CSS-05-v5.bed $wkdir/Prior_gasAcu-results/Jones-et-al-2021-CSS-05-v5_unLifted.bed

## Inversions
liftOver $wkdir/Prior_gasAcu-results/Jones-et-al-2021_inv_locations_v1.bed $wkdir/v1_withChrUn_to_v5.chain.txt $wkdir/Prior_gasAcu-results/Jones-et-al-2021_inv_locations_v5.bed $wkdir/Prior_gasAcu-results/Jones-et-al-2021_inv_locations_v5_unLifted.bed

