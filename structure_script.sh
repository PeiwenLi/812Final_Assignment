#First you have to run STRUCTURE in the cluster
#Make a directory
mkdir STRUCTURE
cd STRUCTURE


#make an extraparams file headers have no "#" while parameters do
nano extraparams
PROGRAM OPTIONS

#define NOADMIX     0 // (B) Use no admixture model (0=admixture model, 1=no-ad
#define LINKAGE     0 // (B) Use the linkage model model
#define USEPOPINFO  0 // (B) Use prior population information to pre-assign ind to clusters
#define LOCPRIOR    0 //(B)  Use location information to improve weak data

#define FREQSCORR   1 // (B) allele frequencies are correlated among pops
#define ONEFST      0 // (B) assume same value of Fst for all subpopulations.

#define INFERALPHA  1 // (B) Infer ALPHA (the admixture parameter)
#define POPALPHAS   0 // (B) Individual alpha for each population
#define ALPHA     1.0 // (d) Dirichlet parameter for degree of admixture
(this is the initial value if INFERALPHA==1).

#define INFERLAMBDA 0 // (B) Infer LAMBDA (the allele frequencies parameter)
#define POPSPECIFICLAMBDA 0 //(B) infer a separate lambda for each pop(only if INFERLAMBDA=1).
#define LAMBDA    1.0 // (d) Dirichlet parameter for allele frequencies




PRIORS

#define FPRIORMEAN 0.01 // (d) Prior mean and SD of Fst for pops.
#define FPRIORSD   0.05  // (d) The prior is a Gamma distribution with these pa$

#define UNIFPRIORALPHA 1 // (B) use a uniform prior for alpha;
otherwise gamma prior
#define ALPHAMAX     10.0 // (d) max value of alpha if uniform prior
#define ALPHAPRIORA   1.0 // (only if UNIFPRIORALPHA==0): alpha has a gamma prior with mean A*B, and
#define ALPHAPRIORB   2.0 // variance A*B^2.


#define LOG10RMIN     -4.0   //(d) Log10 of minimum allowed value of r under li$
#define LOG10RMAX      1.0   //(d) Log10 of maximum allowed value of r
#define LOG10RPROPSD   0.1   //(d) standard deviation of log r in update
#define LOG10RSTART   -2.0   //(d) initial value of log10 r

USING PRIOR POPULATION INFO (USEPOPINFO)

#define GENSBACK    0  //(int) For use when inferring whether an individual is an immigrant, or has an immigrant ancestor in the past GENSBACK generations.  eg, ifGENSBACK==2, it tests for immigrant ancestry back to grandparents.
#define MIGRPRIOR 0.01 //(d) prior prob that an individual is a migrant (used only when USEPOPINFO==1).  This should be small, eg 0.01 or 0.1.
#define PFROMPOPFLAGONLY 0 // (B) only use individuals with POPFLAG=1 to update. This is to enable use of a reference set of individuals for clustering additional "test" individuals.

LOCPRIOR MODEL FOR USING LOCATION INFORMATION

#define LOCISPOP      1    //(B) use POPDATA for location information
#define LOCPRIORINIT  1.0  //(d) initial value for r, the location prior
#define MAXLOCPRIOR  20.0  //(d) max allowed value for r

OUTPUT OPTIONS

#define PRINTNET     1 // (B) Print the "net nucleotide distance" to screen dur$
#define PRINTLAMBDA  1 // (B) Print current value(s) of lambda to screen
#define PRINTQSUM    1 // (B) Print summary of current population membership to$

#define SITEBYSITE   0  // (B) whether or not to print site by site results. (Linkage model only) This is a large file!
#define PRINTQHAT    0  // (B) Q-hat printed to a separate file.  Turn this on before using STRAT.
#define UPDATEFREQ   1000  // (int) frequency of printing update on the screen. Set automatically if this is 0.
#define PRINTLIKES   0  // (B) print current likelihood to screen every rep
#define INTERMEDSAVE 0  // (int) number of saves to file during run

#define ECHODATA     1  // (B) Print some of data file to screen to check that the data entry is correct. (NEXT 3 ARE FOR COLLECTING DISTRIBUTION OF Q:)
#define ANCESTDIST   0  // (B) collect data about the distribution of ancestry coefficients (Q) for each individual
#define NUMBOXES   1000 // (int) the distribution of Q values is stored as a histogram with this number of boxes.
#define ANCESTPINT 0.90 // (d) the size of the displayed probability interval on Q (values between 0.0--1.0)



MISCELLANEOUS

#define COMPUTEPROB 1     // (B) Estimate the probability of the Data under the model.  This is used when choosing the best number of subpopulations.
#define ADMBURNIN  500    // (int) [only relevant for linkage model]: Initial period of burnin with admixture model (see$
#define ALPHAPROPSD 0.025 // (d) SD of proposal for updating alpha
#define STARTATPOPINFO 0  // Use given populations as the initial condition for population origins.  (Need POPDATA==1).  It is assumed that the PopData in the input file are between 1 and k where k<=MAXPOPS.
#define RANDOMIZE      1  // (B) use new random seed for each run
#define SEED        2245  // (int) seed value for random number generator (must set RANDOMIZE=0)
#define METROFREQ    10   // (int) Frequency of using Metropolis step to update Q under admixture model (ie use the metr. move eve i steps).  If this is set to 0, it is never used. (Proposal for each q^(i) sampled from prior.  The goal is to improve mixing for small alpha.)
#define REPORTHITRATE 0 //   (B) report hit rate if using METROFREQ

___________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

#Now generate a mainparams_kx for each value of K
nano mainparams_k2


Basic Program Parameters

#define MAXPOPS    2	  // (int) number of populations assumed
#define BURNIN    100000   // (int) length of burnin period
#define NUMREPS   1000000   // (int) number of MCMC reps after burnin

Input/Output files

#define INFILE   StoneCleb_contact.str   // (str) name of input data file
#define OUTFILE  out_StoneCleb_contact.str_1M_K2  //(str) name of output data file

Data file format

#define NUMINDS    20    // (int) number of diploid individuals in data file
#define NUMLOCI    939    // (int) number of loci in data file
#define PLOIDY       2    // (int) ploidy of data
#define MISSING     -9    // (int) value given to missing genotype data
#define ONEROWPERIND 0    // (B) store data for individuals in a single line


#define LABEL     1     // (B) Input file contains individual labels
#define POPDATA   1     // (B) Input file contains a population identifier
#define POPFLAG   0     // (B) Input file contains a flag which says
                              whether to use popinfo when USEPOPINFO==1
#define LOCDATA   0     // (B) Input file contains a location identifier

#define PHENOTYPE 0     // (B) Input file contains phenotype information
#define EXTRACOLS 0     // (int) Number of additional columns of data
                             before the genotype data start.

#define MARKERNAMES	 0  // (B) data file contains row of marker names
#define RECESSIVEALLELES 0  // (B) data file contains dominant markers (eg AFLPs)
                            // and a row to indicate which alleles are recessive
#define MAPDISTANCES     0  // (B) data file contains row of map distances
                            // between loci


Advanced data file options

#define PHASED           0 // (B) Data are in correct phase (relevant for linkage model only)
#define PHASEINFO        0 // (B) the data for each individual contains a line
                                  indicating phase (linkage model)
#define MARKOVPHASE	 0 // (B) the phase info follows a Markov model.
#define NOTAMBIGUOUS  -999 // (int) for use in some analyses of polyploid data

_____________________________________________________________________________________________________________________________________________________________

#Now write a short sbatch script 

#SBATCH -J alignment
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH --mem=6g
#SBATCH -t 14-0:0:0
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=15pl16@queensu.ca
#SBATCH -o alignment_%j.o
structure -m mainparams_k2 -e extraparams -K 2
structure -m mainparams_k3 -e extraparams -K 3
structure -m mainparams_k4 -e extraparams -K 4

#Before you submit it go back and make a replicate of your whole STRUCTURE directory for the number of replicates you are running so ... STRUCTURE_rep1,  STRUCTURE_rep2, STRUCTURE_rep3, etc
#Now cd into each one and bsub each script
#Each of mine seem to finish in a day or two
______________________________________________________________________________________________________________________________________________________________
                                                  Y
                                                  ^
                                                |..|
                                               \_ | _/
                                                  X
                                                  X
                                              _/\ ^/\_
                                                  ^
                                                  ^
                                                  \_/
___________________________________________________________________________________________________________________________________________________________


#alright now time for r we'll be using pophelper http://royfrancis.github.io/pophelper/ 
#You should be able to copy and paste this into a R script and run it (check inputs and k values)

setwd ("/home/nick/PhD/ClineAnalysis/STRUCTURE/StoneCleb/")
getwd()


# install dependencies and devtools
#install.packages(c("Cairo","ggplot2","gridExtra","gtable","tidyr","devtools"),dependencies=T)
#install.packages("Cairo")

# install pophelper package from GitHub
devtools::install_github('royfrancis/pophelper')

# load library for use
library(pophelper)
library(stringr)
_______________________________________________________________________________________________________________________________
#I have run structure with 5 replicates now is time to run CLUMPP
#Import the files, I found it easiest to do this by calling the path and listing the files then 

sfiles <- list.files(path = "/home/nick/PhD/ClineAnalysis/STRUCTURE/StoneCleb/results_StoneCleb/")
sfiles
#note do not include the _q files as these screw things up downstream.
sfiles<- file.path("/home/nick/PhD/ClineAnalysis/STRUCTURE/StoneCleb/results_StoneCleb/", c("out_StoneCleb_contact_1M_K2_f", "out_StoneCleb_contact_1M_K3_f", "out_StoneCleb_contact_1M_K4_f", "out_StoneCleb_contact_run2_1M_K2_f", "out_StoneCleb_contact_run2_1M_K3_f", "out_StoneCleb_contact_run2_1M_K4_f", "out_StoneCleb_contact_run3_1M_K2_f", "out_StoneCleb_contact_run3_1M_K3_f", "out_StoneCleb_contact_run3_1M_K4_f", "out_StoneCleb_contact_run4_1M_K2_f", "out_StoneCleb_contact_run4_1M_K3_f", "out_StoneCleb_contact_run4_1M_K4_f"))
#check
sfiles
# now turn this list of files into a 'qlist' object
slist <- readQ(files=sfiles)
#Let's get some info on these runs
#tabulateQ() produces a table of runs with various parameters
tr1<- tabulateQ(slist, writetable=TRUE)
#which is further simplified by summariseQ()
sr1<-summariseQ(tr1, writetable=TRUE)

evannoMethodStructure(summariseQ(tabulateQ(slist)))

________________________________________________________________________________________________________________________________

#now run CLUMPP, qlist is the input files, Prefix is the what you want to append to the outputs, 
#Paramode is the clumpp algorithm used typacally the best is GREEDY_OPTION = 2, 
#paramrep is the number of times the runs will be run and resampled
clumppExport(qlist=slist, prefix=, parammode=2, paramrep=5000, useexe=T)
#check if they look similar for each run
aligned_k4 <- readQ("./pop_K4/pop_K4-combined-aligned.txt")
plotQ(aligned_k4,imgoutput="join")
aligned_k3 <- readQ("./pop_K3/pop_K3-combined-aligned.txt")
plotQ(aligned_k3,imgoutput="join")
aligned_k2 <- readQ("./pop_K2/pop_K2-combined-aligned.txt")
plotQ(aligned_k2,imgoutput="join")

collectClumppOutput(prefix="pop", filetype="aligned", runsdir=NA, newdir=NA, quiet=FALSE)
collectClumppOutput(prefix="pop", filetype="merged", runsdir=NA, newdir=NA, quiet=FALSE)



_______________________________________________________________________________________________________________________
#so at this point I have already run structure and clumpp
#to generate the clumpp inputs I had to copy and paste the Proportion of membership of each pre-defined population in each of the x clusters then add put them together in the same file. 
#should be formated like this:
#1_NC1545:      0.9996  0.0003  0.0002          1
#1_NC1552:      0.9997  0.0002  0.0001          1
#1_NC1553:      0.9997  0.0002  0.0001          1
#1_NC1511:      0.9996  0.0002  0.0002          1
#1_NC1520:      0.9993  0.0003  0.0003          1
#1_NC1527:      0.9997  0.0002  0.0001          1

#Run it through CLUMPP using clumppExport()


#use the pop output to generate a rough plot
K2 <- readQ("./pop_K2-combined-merged.txt")
plotQ(K2, showindlab=TRUE, sortind="all", clustercol=c("navyblue","gold2","darkgreen"))
plotQMultiline(K2, showindlab=TRUE, clustercol=c("navyblue","darkgreen"))

K4 <- readQ("./Pcru36_popk4.popq")
plotQ(aligned4)

#*********************be careful using the merged runs unless they are in complete agreement*****************************



