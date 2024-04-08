# BOSTIn
BOSTIn (Broad Overview of Sequence & Topology Incongruence) is a new, user-friendly phylogenetic artifact identification package.
![Screen_Shot_2024-03-22_at_4](https://github.com/JFFleming/BOSTIn/assets/84772413/352a827b-1935-4722-96c5-8bcf96a93263)

# Installing BOSTIn
# Prequisites
BOSTIn requires Perl and R. 
Many distributions already come bundled with Perl, but if you don't have access to it (which you can find out by typing "which perl" into your command line terminal), you can install it here:
https://www.perl.org/get.html
To install R, a download link can be found here:
https://cran.r-project.org/

In addition, BOSTIn requires one Perl module (BIO::SeqIO), and one R package (Phangorn). Once you have installed Perl and R, you can install those packages with these two commands:

```
cpanm install BIO::SeqIO
Rscript -e 'install.packages("phangorn", repos="https://cloud.r-project.org")'
```

# Making BOSTIn work
Once you have the prerequisites installed, then BOSTIn will work provided you have all of the BOSTIn scripts in the same directory as the Fasta file. However, moving them around your system every time you want to do an analysis is a pain, so this next step turns BOSTIn into an executable you can use anywhere on your computer!
First, place the BOSTIn directory in the place you want to keep it. Use
```
pwd
```
to find the complete address of that directory.

In your home directory (~/), you will find a file called either .bashrc or .bash_profile. Add this line to the bottom of that file.
```
export PATH=$PATH:<pwd address>/BOSTIn_v1_0
```
This makes sure that the next time you open your terminal, it registers the BOSTIn directory as a place to look for executables.
Finally, go inside the BOSTIn_v1_0 directory, with all the scripts, and run the following command:
```
chmod u+x+r *
```
This elevates these files to make sure they can do what you need them to do. If you are on a Mac, you might also need to run this
```
xattr -d com.apple.quarantine *
```
Which removes the "quarantine" status from the software package - basically telling your mac it isn't harmful.

Then, you should be ready to go! If you're having trouble installing BOSTIn, please don't hesitate to get in touch through the Github!

# Using BOSTIn
# The BOSTIn command
The BOSTIn command is formatted as follows:

```
Bostin.pl <options> <datatype> <input fasta file> <output file prefix>
```
For Example:
```
Bostin.pl --blh --ch --s protein test.fa test
```
runs analyses to detect branch length heterogeneity, compositional heterogeneity and site saturation tests on the protein dataset of test.fa, and all the output files will start with test.

Meanwhile,
```
Bostin.pl --ch dna nucleotide.fa hello
```
runs just an analysis to detect compositional heterogeneity on the dna dataset nucleotide.fa, and all the output files will start with hello

We will break down each step.

# Options
If you want to use more options, just add them one after the other before you insert the datatype. The order of options does not matter.

--blh : this option assesses the dataset for Branch Length Heterogeneity using the LB-Score, as per Struck (2014) (see references)
--ch : this option assesses the dataset for Compositional Heterogeneity using nRCFV, as per Fleming & Struck (2023) (see references)
--s : this option assesses the dataset for Site Saturation using the DE-Score, which is an experimental new methodology which will have a preprint soon! I describe the concept in more detail in the DE-Score file in this GitHub. Currently, the --s option is not enabled for nucleotide datasets, but our intention is to implement the C-Value in the coming weeks.

# Datatype
Datatype has two options:
protein
dna
protein is for amino acid datasets, while dna is for nucleotide datasets.

# Results from BOSTIn
# The Results Files
BOSTIn produces a huge number of results files, and some of these will be of more general use than others. For each option:

# Branch Length Heterogeneity

\<prefix\>.BranchHeterogeneity.NarrativeReport.txt  - This is the narrative report explaining your Branch Length Heterogeneity results (see narrative reports section below)

\<prefix\>.LBScore.Taxa.redflags.txt - This is a list of the taxa that BOSTIn identifies as "red flags" in your dataset based on branch length heterogeneity - ones that should be considered for removal or alteration.

\<prefix\>.LBScore.Taxa.yellowflags.txt	- This is a list of the taxa that BOSTIn identifies as "yellow flags" in your dataset based on branch length heterogeneity - ones that could cause problems.

\<prefix\>.LBSummary.histogram.pdf - This is a histogram of the taxon-specific LB-Scores, so that you can see what the overal distribution of branch length heterogeneity looks like in your dataset.

\<prefix\>.LBSummary.txt - This is a summary of the whole dataset LB-Scores - the mean and quartile values.

\<prefix\>.LBi-scores - This is a list of the LB-Scores for each taxa

# Compositional Heterogeneity
\<prefix\>.CompositionalHeterogeneity.NarrativeReport.txt  - This is the narrative report explaining your CompositionalHeterogeneity results (see narrative reports section below)

\<prefix\>.Frequencies.txt - This is a list of the basic frequencies of each amino acid or nucleotide across each taxa in the dataset

\<prefix\>.RCFV.txt - This shows the nRCFV and RCFV of the dataset


\<prefix\>.ncsRCFV.txt - This is a simple list of the normalised csRCFV score for each character in your dataset.

\<prefix\>.ncsRCFV.CompHet.Character.redflags.txt - This is a list of the characters that BOSTIn identifies as "red flags" in your dataset based on compositional heterogeneity - ones that should be considered for removal or alteration.

\<prefix\>.ncsRCFV.CompHet.Character.yellowflags.txt - This is a list of the characters that BOSTIn identifies as "yellow flags" in your dataset based on compositional heterogeneity - ones that could cause problems.

\<prefix\>.ncsRCFV.histogram.pdf - This is a histogram of the normalised csRCFV scores, so that you can see what the overal distribution of compositional heterogeneity looks like across characters in your dataset.

\<prefix\>.ncsRCFV.summary.txt - This shows a list of summary statistics for the normalised csRCFV scores. The mean and median, alongside the standard deviation and median absolute deviation, help us understand the distribution of compositional heterogeneity across a dataset.


\<prefix\>.ntRCFV.txt - This is a simple list of the normalised tRCFV score for each taxa in your dataset.

\<prefix\>.ntRCFV.CompHet.Taxa.redflags.txt - This is a list of the taxa that BOSTIn identifies as "red flags" in your dataset based on compositional heterogeneity - ones that should be considered for removal or alteration.

\<prefix\>.ntRCFV.CompHet.Taxa.yellowflags.txt - This is a list of the taxa that BOSTIn identifies as "yellow flags" in your dataset based on compositional heterogeneity - ones that could cause problems.

\<prefix\>.ntRCFV.histogram.pdf - This is a histogram of the normalised tRCFV scores, so that you can see what the overal distribution of compositional heterogeneity looks like across taxa in your dataset.

\<prefix\>.ntRCFV.summary.txt - This shows a list of summary statistics for the normalised tRCFV scores. The mean and median, alongside the standard deviation and median absolute deviation, help us understand the distribution of compositional heterogeneity across a dataset.


# Site Saturation

\<prefix\>.SiteSaturation.NarrativeReport.txt - This is the narrative report explaining your Site Saturation results (see narrative reports section below)

\<prefix\>.SiteSaturation.TotalFrequencies.txt - This is a file containing the DE-score of your entire dataset.

\<prefix\>.SiteSaturation.TaxaFrequencies.SiteSat.Taxa.redflags.txt - This is a list of the taxa that BOSTIn identifies as "red flags" in your dataset based on site saturation - ones that should be considered for removal or alteration.

\<prefix\>.SiteSaturation.TaxaFrequencies.SiteSat.Taxa.yellowflags.txt - This is a list of the taxa that BOSTIn identifies as "yellow flags" in your dataset based on site saturation - ones that could cause problems.

\<prefix\>.SiteSaturation.TaxaFrequencies.histogram.pdf - This is a histogram of the taxa-specific DE-Scores scores, so that you can see what the overal distribution of site saturation looks like across taxa in your dataset.

\<prefix\>.SiteSaturation.TaxaFrequencies.summary.txt This shows a list of summary statistics for the taxa-specific DE-scores. The mean and median, alongside the standard deviation and median absolute deviation, help us understand the distribution of site saturation across a dataset.

\<prefix\>.SiteSaturation.TaxaFrequencies.txt - This is a simple list of the DE-scores for each taxa in your dataset.


# The Narrative Report
Unique to BOSTIn is the Narrative Report function. This is intended as a way to explain not just which sequences are good or bad, but also why! No more hunting through opaque documentation, instead, BOSTIn customises answers for your dataset, and tells you why it has made the decisions it has made.
BOSTIn implements a number of simple summary statistic-based tests to identify red and yellow flagged taxa or characters, depending on the metric. What that means will differ from dataset to dataset. Red flags generally require attention, while yellow flags might be ignorable, especially if a taxa is particularly important to the dataset.
The Narrative report will also suggest potential next steps, such as data recoding, and will explain why that might be useful (and potential risks, too), so you can be informed about how you approach your data.
It's the part we're most excited about, but also the part that needs the most fine-tuning, so if you have any suggestions for extra details you'd like to see, do get in touch through the GitHub!


# References
Struck, T.H., 2014. TreSpEx–-detection of misleading signal in phylogenetic reconstructions based on tree information. Evolutionary Bioinformatics, 10, pp.EBO-S14239.

Fleming, J.F. and Struck, T.H., 2023. nRCFV: a new, dataset-size-independent metric to quantify compositional heterogeneity in nucleotide and amino acid datasets. BMC bioinformatics, 24(1), p.145.
