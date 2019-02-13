# Population Genetic Structure Analysis

This exercise uses data from two host-specialist species of ladybirds from Japan that are morphologically very similar: N=Henosepilachna niponica, a specialist on thistle (Cirsium sp, Asteraceae) and Y= H. yasutomii, a specialist on blue cohosh (Caulophyllum robustum, Berberiaceae). These data are generated from samples taken where the species’ ranges overlap, one separate population of each and one of each where they occur together.

ddRAD data were generated and here we will use a subset of the full data set (2,000 loci in 20 individuals, 5 from each population) to investigate if there is hybridisation in sympatry.

We will evaluate population genetic structure using two methods:
1.	Running a Bayesian cluster algorithm implemented in the program STRUCTURE.
2.	Running a Maximum-likelihood model based on genotype likelihoods in the program NGSADMIX.


For the analysis we will connect to Sheffield University’s high performance computing server, known as iceberg, through a [web browser](https://www.sheffield.ac.uk/wrgrid/using/access). The software we will use is installed on the server, and we will interact with it through the command line interface. Command-lines to be typed are indicated throughout the protocol are coded in Markdown, e.g. to navigate to your home directory type:

```markdown
cd
```

Log into iceberg. In order to access the data & software for the analysis you need to log in to a worker node: 

```markdown
qsh 
```

You will be logged in to your home directory, `/home/USERNAME/` where `USERNAME` represents your University of Sheffield or guest account username, e.g. cs4abxx. We will navigate through directories using the command `‘cd’` (change directory). To check which directory you are in you can use the command to print the working directory type:

```markdown
pwd
```

## 1. Run STRUCTURE on Computing Cluster as array job

STRUCTURE can be run with population genomic data, although the runs can take time especially if 

a) models with large K are tested

b) many individuals were sampled

c) many markers are used

d) many repeats are attempted

Make a new folder in your home directory in which to work and navigate into it.

type:

```markdown
mkdir Structure
cd Structure
```

To copy the data and analysis scripts from the shared directory, type:

```markdown
cp –r /usr/local/extras/Genomics/workshops/February2019/Structure/* ./
```

In this part you will start some dummy runs (with very few generations) to learn how to use a job submission script. Runs with sufficient number of generations (>100k would take unfortunately about 24 hours). The structure input file has been created from a vcf file of 2000 loci using PGDSpider. The runs with real results will be in the folder 'runs'.

Look at the `mainparams` file first - this contains the necessary specifications of the STRUCTURE model run (see also manual for more details).

```markdown
less mainparams
```

- How many individuals are used?
- How many markers?
- What are burn-in and run lengths?

Press `‘q’` to return to the command prompt.

The `.sh` file is a small script used to submit a job to the cluster. Submitted jobs will run remotely and leave your worker node available for other tasks, they will not necessarily start immediately. Parallel jobs or an array of jobs have to be submitted to the queue, this is how we will run STRUCTURE. To look at the submission job file:

```markdown
less TGK2.sh
```


```markdown
#!/bin/bash

#Settings for the Sun Grid Engine

#$ -l h_rt=7:59:59
# = maximum length of job

#$ -l rmem=2G
# maximum memory available for each core

#$ -P molecol
#$ -q molecol.q
# = job is referred to NBAF node 'molecol' which has 64 cores. If this option is ‘commented out’ with a double hash, the job will be assigned to an alternative node

#$ -t 1-3
# number of jobs run in array: 3 - these will be run roughly at the same time

#$-N TGK2
#=name of the job

sleep "$SGE_TASK_ID"0
#=we specified in the mainparams file that seed is generated randomly based on system clock
#to avoid that all jobs have the same seed, the will start after each other with a short 'sleep' period in between

structure -K 2 -i TG_variants_str_2000 -o /home/$USER/Structure/K2_$SGE_TASK_ID > /home/$USER/Structure/runseqK2_$SGE_TASK_ID
#the structure command line
#using the command line one can override some settings of the mainparams file, for example we specified different in/out files 
#in the command line version of STRUCTURE every K needs to be run separately (in our case we start with K = 2)

```

Press `‘q’` to return to the command prompt. To submit the job use `qsub`:

```markdown
qsub TGK2.sh
Qstat 
```

Check for the status of all your jobs (qw = waiting, r = running). 'QRLOGIN' is the status of the worker node that you are working on. `Qstat` provides the time & date when job was submitted. If the job is not showing up anymore it is finished (or has failed). You should see three lines for job 'TGK2' one for each run within the array when the job is running. Successful jobs will create a `_f` file with results and a `runseqK` file that contains the usual screen output of STRUCTURE and can be useful to track the progress of a run. `TGK2.oXXX` and `TGK2.eXXXX` files are other output or error files although they should remain empty.

We can use the first script as a template to prepare .sh files for other values of K using the sed ‘find and replace command’ – finding ‘K2’ (or ‘K 2’) and replacing it with ‘K3’ (or ‘K 3’).

```markdown
sed 's/K2/K3/g' TGK2.sh > TGK3.sh
sed –i 's/K 2/K 3/g' TGK3.sh

sed 's/K2/K4/g' TGK2.sh > TGK4.sh
sed –i 's/K 2/K 4/g' TGK4.sh
```

Open the `.sh` files to check what changes have been made before continuing, then submit each `.sh` file to the queue:

```markdown
qsub TGK3.sh
qsub TGK4.sh
```

It will take a few mins until the runs are finished, once they have check that the `_f-files` have been created and have a look at them.

**Evaluating the results from the longer runs we prepared earlier ...**

```markdown
cd runs
grep 'Estimated Ln Prob of Data' K*f
```

These are the Ln Probabilities of each independent run for each value of K. We will use the output with the highest likelihood for each of the values of K tested to generate bar plots to visualise the genetic structure across the four populations. These are `K2_1_f`, `K3_4_f`, `K4_2_f`.

We will visualise the STRUCTURE results using the software DISTRUCT. We will extract the cluster membership proportions for each individuals from the three results files listed above.

```markdown
grep -A 20 Label K2_1_f | sed "1d" > ../K2.outfile
grep -A 20 Label K3_4_f | sed "1d" > ../K3.outfile
grep -A 20 Label K4_2_f | sed "1d" > ../K4.outfile
```

We also need several parameter files:

a) `.perm` file which specifies the colour for each genetic cluster

b) `.qfile` which is needed to run Distruct and usually has the population averages. This is critical if you want to draw population averages, since we want to draw each individual we can use a dummy (e.g. `TG.qfile`)

c) `.group` file which has the labels for the four different populations

d) a `drawparams` which contains the names of each parameter file and also plot characteristics 

Navigate back into the Structure directory and have a look at the drawparams file:

```markdown
cd ..
less drawparams
```

For plotting two clusters there are specific settings to point the software to the `K2.outfile`, the `K2TG.perm` file, the output to save as `K2TG.ps` and the number of clusters being plotted as `K 2`. Press `‘q’` to return to the command prompt, and to generate the plot:

```markdown
distruct
```

type: `'ls'` to check that the `K2TG.ps` file has been created. Then go on to create two further plots for K=3 and K=4. To generate plots for three and four clusters we will modify the `drawparams` file.

```markdown
nano drawparams
```

Change every K2 to K3 and change the K-parameter on the line that starts `#define K`, save as the same name and re-run distruct. Repeat this for K=4.

You should now have three plots. If you are using an interactive session (with qsh) you can open this in the gv viewer:

```markdown
gv K2TG.ps
```

Alternatively you can convert these files to pdf and email them to yourself to view:

```markdown
ps2pdf K2TG.ps
echo "Text body" | mail -s "Subject: K2 Structure plot" -a K2TG.pdf your@email
```

## 2. Run NGSADMIX using genotype likelihoods

NGSadmix uses genotype likelihood information rather than absolute genotype calls. Genotype likelihood information is often output with the genotype in a vcf file as either ‘GL’ = Genotype Likelihood, or ‘PL’ = phred-scaled genotype likelihoods. The input file for NGSAdmix should be in beagle format. For this exercise that has been generated using vcftools with the following command:
`vcftools  --vcf variants_2000.vcf --out TG2000 --BEAGLE-PL --chr pseudoscaff_000010`

This generates a file of genotype likelihood values for each SNP for each sample, named `TG2000.BEAGLE.PL`. The command lines to run NGSadmix for K=2, K=3 and K=4 are in the script ‘NGSa.sh’. To visualise this:

```markdown
less NGSa.sh
```

The settings are similar to those for the structure submission scripts, although for the purpose of this exercise we will run each value of K only once. Press `‘q’` to return to the command prompt, then submit the job to run:

```markdown
qsub NGSa.sh
```

The job should complete quite quickly after it has started to run and we will use the .qopt files, plus the supplied pop.info file, for visualising the results. The R commands for generating the plots are in the file `plots.R`, which you can view in the usual way with `‘less’` before opening R.

Start R:

```markdown
R
```

Within R:

```markdown
source("plots.R")
q()
```

You can email the pdf file to yourself to open it 

```markdown
echo "Text body" | mail -s "Subject: NGSAdmix plot" -a Rplots.pdf your@email
```

Alternative if you are using an interactive (qsh) session: 

```markdown
echo gv Rplots.pdf
```

**Software used**

[STRUCTURE](https://web.stanford.edu/group/pritchardlab/structure.html)

[DISTRUCT](https://web.stanford.edu/group/rosenberglab/distruct.html)

[NGSADMIX](http://popgen.dk/software/index.php/NgsAdmix)
