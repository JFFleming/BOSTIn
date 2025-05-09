#!/bin/bash

# Initialize an array to hold the selected options
OPTIONS=()

# Parse options
while [[ "$1" == --* ]]; do
    case "$1" in
        --blh|--s|--ch)
            OPTIONS+=("$1")
            shift
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Make sure we have exactly two arguments left
if [ $# -ne 3 ]; then
    echo "Usage: $0 [--blh] [--s] [--ch] <dna or protein> <list of fasta files> <genewise output prefix>"
    exit 1
fi

ARG1="$1"
ARG2="$2"

# Run Bostin.pl with each relevant argument
for i in $(cat "$ARG2"); do
    Bostin.pl "${OPTIONS[@]}" "$ARG1" "$i" "$i"
done

# Now process the options
for opt in "${OPTIONS[@]}"; do
    case "$opt" in
        --ch)
            grep 'nRCFV' *.RCFV.txt | sed -e s'/.RCFV.txt:nRCFV//'g -e '1i\'$'\n''Fasta File	nRCFV'$'\n' -e s'/    /\t/' > $3.nRCFV.fastaSummary.txt;
            Genewise.nRCFVReport.r $3.nRCFV.fastaSummary.txt
            ;;
        --blh)
			grep 'standardDeviation' *.LBSummary.txt | sed -e s'/.LBSummary.txt:standardDeviation//'g -e '1i\'$'\n''Fasta File	LB-Score Standard Deviation'$'\n' -e s'/    /\t/' > $3.LBScore.fastaSummary.txt;
            Genewise.BLHReport.r $3.LBScore.fastaSummary.txt
            ;;
        --s)
            if [ "$ARG1" = "protein" ]; then
				grep 'DE-Score:' *.SiteSaturation.TotalFrequencies.txt | sed -e s'/.SiteSaturation.TotalFrequencies.txt:DE-Score://'g -e '1i\'$'\n''Fasta File	DE-Score'$'\n' -e s'/    /\t/' > $3.DEScore.fastaSummary.txt;
                Genewise.DEScoreReport.r $3.DEScore.fastaSummary.txt
            elif [ "$ARG1" = "dna" ]; then
				grep 'CScore' *.SiteSaturation.TotalCScore.txt | sed -e s'/.SiteSaturation.TotalFrequencies.txt:CScore//'g -e '1i\'$'\n''Fasta File	CScore'$'\n' -e s'/    /\t/' > $3.CScore.fastaSummary.txt; 
                Genewise.CScoreReport.r $3.CScore.fastaSummary.txt
            fi
            ;;
    esac
done
