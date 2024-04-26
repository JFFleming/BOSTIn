function main(compHet, branchLength, siteSat, args) {
    console.log("Other things found on the command line:");
    args.forEach(arg => console.log(arg));

    if (compHet) {
        doCompHet(args);
    }

    if (branchLength) {
        if (siteSat) {
            if (args[0].includes('dna')) {
                console.log("HAHA FOUND YOU!\n We do a super analysis here as both need the NJ tree");
                doBranchLengthSiteSat(args);
            } else {
                doBranchLength(args);
            }
        } else {
            doBranchLength(args);
        }
    }

    if (!args[0].includes('dna') && !branchLength) {
        if (siteSat) {
            doSiteSat(args);
        }
    }
}

function doCompHet(args) {
    if (args[0].includes('dna')) {
        console.log("RCFV Reader DNA!");
        console.log(`RCFV_Reader.pl dna ${args[1]} ${args[2]}`);
        console.log(`nRCFVAnalysis_Script.R ${args[2]}.ntRCFV.txt ${args[2]}.ncsRCFV.txt`);
        console.log(`CompHetNarrativeReport.pl ${args[2]} dna > ${args[2]}.CompositionalHeterogeneity.NarrativeReport.txt`);
    } else if (args[0].includes('protein')) {
        console.log("RCFV Reader Protein!");
        console.log(`RCFV_Reader.pl protein ${args[1]} ${args[2]}`);
        console.log(`nRCFVAnalysis_Script.R ${args[2]}.ntRCFV.txt ${args[2]}.ncsRCFV.txt`);
        console.log(`CompHetNarrativeReport.pl ${args[2]} protein > ${args[2]}.CompositionalHeterogeneity.NarrativeReport.txt`);
    } else {
        console.log(`To run a Compositional Heterogeneity analysis in BOSTIn, you need to specify whether you are using Amino Acid or Nucleotide data. For amino acid data, use:

        BOSTIn.pl -C protein <AlignmentFile> <Prefix for Output files>

        For nucleotide data:

        BOSTIn.pl -C dna <AlignmentFile> <Prefix for Output files>`);
    }
}

function doBranchLength(args) {
    if (args[0].includes('dna')) {
        console.log("LB-Score DNA!");
        console.log(`CombinedNJ_LBScore.R ${args[1]} ${args[2]}.tre ${args[2]}.LBi-scores ${args[2]}.LBSummary.txt DNA`);
        console.log(`LBNarrativeReport.pl ${args[2]} > ${args[2]}.BranchHeterogeneity.NarrativeReport.txt`);
    } else if (args[0].includes('protein')) {
        console.log("LB-Score AA!");
        console.log(`CombinedNJ_LBScore.R ${args[1]} ${args[2]}.tre ${args[2]}.LBi-scores ${args[2]}.LBSummary.txt AA`);
        console.log(`LBNarrativeReport.pl ${args[2]} > ${args[2]}.BranchHeterogeneity.NarrativeReport.txt`);
    } else {
        console.log(`To run a Branch Length analysis in BOSTIn, you need to specify whether you are using Amino Acid or Nucleotide data. For amino acid data, use:
        BOSTIn.pl -C protein <AlignmentFile> <Prefix for Output files>

        For nucleotide data:

        BOSTIn.pl -C dna <AlignmentFile> <Prefix for Output files>`);
    }
}

function doSiteSat(args) {
    if (args[0].includes('dna')) {
        console.log("Site Saturation DNA Here. Unfortunately, it has not yet been implemented. It will be coming soon!");
        console.log("then we append the narrative to narrative results.");
        console.log("then we evaluate it with an R Script.");
    } else if (args[0].includes('protein')) {
        console.log("Site Saturation AA!");
        console.log(`DEScoreCalculator.pl ${args[1]} ${args[2]}`);
        console.log(`DEScoreAnalysis_Script.r ${args[2]}.SiteSaturation.TaxaFrequencies.txt`);
        console.log(`DEScore.NarrativeReport.pl ${args[2]} > ${args[2]}.SiteSaturation.NarrativeReport.txt`);
    } else {
        console.log(`To run a site saturation analysis in BOSTIn, you need to specify whether you are using Amino Acid or Nucleotide data. For amino acid data, use:

        BOSTIn.pl -s protein <AlignmentFile> <Prefix for Output files>

        For nucleotide data:

        BOSTIn.pl -s dna <AlignmentFile> <Prefix for Output files>`);
    }
}

function doBranchLengthSiteSat(args) {
    console.log("combination Site Sat Branch Length analysis for DNA to save NJ Tree construction time has not yet been implemented.");
}

// Example usage
main(true, true, false, ["dna", "alignment.fasta", "output"]);
