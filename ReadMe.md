# Genetic Methods in Aquaculture
## Session: Data analysis
#### Date: January
#### Material below has been made by Dorte Bekkevold and Jakob Hemmer-Hansen and adapted by Belen Jimenez-Mena

## Schedule of the lecture:

This Data Analysis session will consist on the following two parts, each with a short introduction to the case study and a practical part:
- Part 1: Brief introduction to R
- Part 2: Parentage analysis
- Part 3: Population genetics analysis

The session is thought as a dynamic session - so discussions and questions throughout the practicals are encouraged. All the data can be found at DTULearn. 

## Practical 1: Introduction to R
Throughout this session, you require the use of the programming language R. To install R and Rstudio, see [here](installing_R) for a short manual on instructions how to set up R and Rstudio in your laptop.

### The objectives of this part are:
* To get R installed in our computer and use an interface called Rstudio
* To open a genetic file in R
* To install and load a package

### Dataset and tutorial
As mentioned, we are going to use the programming language R for the analysis and visualization. R is a programming language suited for statistical computing that has been developed by the scientific community and it is widely used for data analysis. 

OK. First, set up a folder in your computer where you will save the R script and the results generated in each of the exercises.
Now, we will open R, save the script into our folder. In R, appart from coding our own functions and programs, we can download packages and use the functions that have been developed and saved there. This is one of the advantages of using R, its active user community that helps us not re-invent the wheel.
So now we will load a few packages in R that contain the functions that we need. For that, write and run the following:
```
install.packages("<the package's name>")
```
R will then download the package from Internet. To load the package into your current session, you run:
```
library("<the package's name>")
```
or
```
require("<the package's name>")
```

## Practical 2: Parentage analysis 
The objective of this part is to identify the sire and dam of a list of offsprings in your aquaculture line. You need to evaluate the relatedness of the offsprings and assess its viability for future crosses. To help you do this, you can use a software called COLONY (https://www.zsl.org/science/software/colony) installed in your laptop.

### The objectives of this part are:
* To recognize the different codes that a genetic file has in relation to the parentage analysis
* To choose the best parameter criteria for the parentage analysis
* To interpret the results of a parentage analysis
* To assign individuals back to their most probable sire and dam
* To evaluate and assess how best to optimize and set-up an aquaculture line in terms of relatedness


### Dataset and tutorial
We will use a dataset consisting on 1 family material fish from mixed-family tansk (which represent families of released salmon fry currently on spawning sites or in Arctic waters, file: "YourSalmongenotypes") and a baseline of all the parents that were used as broodstock to produce the crosses in 2015 (files: 40MUMSinput and 59DADSinput). The genotypes come from a set of 16 microsatellite markers and their error type rate information is in the file "16MarkerTypeErrorRate".

Before loading the datasets into COLONY, open the input files and explore the following questions (in pairs):
- What does each row and column correspond to? 
- How is the code for each marker? And for missing data?

Once you are ready, we will go together through the loading of the files and options into COLONY. Once all files are loaded, please Run the software (1 Run will be enough for this practical but you are very welcome to try a larger number after the session).

To explore the objective of this analysis, you can look at the following Results tab ("View Results" option):
- Best (ML) Family Cluster, Best (ML) Paternity Assignment and Best (ML) Maternity Assignment: to identify the most probable sire and dam
- Best (ML) Sibship Assignment: to explore number of fullsibs and half sibs
Explore the rest of the Results' tabs if you need when discussing the questions.

Questions to discuss in pairs:
- Which is the most probable sire and dam for each offspring? Is there only one probable dam for each individual? What about sire?
- What is the proportion of half- and full- sibs for these crossings?
- What is the effective population size (Ne) of the offspring samples? What does it mean such a low Ne?

## Practical 3: Population genetics analysis
We will use a genetic marker dataset of salmon from 3 different locations: 
The individuals used for this study represent the baseline that we will be using for our analysis today. The data corresponds to the .gen file. The .gen format is a widely-used format in population genetics. The data can be found at DTULearn. 

### The objectives:
The objectives of this part are:
* To recognize a file in the genetic format and identify its parts
* To define what a population is
* To apply different genetic tools in a real dataset
* To recognize different populations in a sample of different individuals, and quantify their differences
* To assign individuals back to their population of origin

### Dataset and tutorial
#### Data exploration
As a first step, it is always a good idea to explore how a genetic file looks like. Please, open the "baseline_data.gen" file and spend a few minutes exploring how the format is.
Pay special attention to:
* How are the alleles coded?
* Is there missing data? Is there a lot of it?
Discuss with your partner.

#### Packages needed for today's exercises
There are thousands of helpful R packages for you to use. For the analysis at this lecture, we will be using the following packages. If you work in the cluster, you don't have to install the packages, as they are already installed (you only need to load them using require(), but if you had to install them again, you just uncomment the install.packages() commands). If you work from your laptop, you'd need to install them if you don't have them yet.
```
require(adegenet)
require(diveRsity)
require(hierfstat)
require(rapport)
require(graph4lg)
library(dartR)
require(StAMPP)
require(HardyWeinberg)
require(genetics)
require(pegas)
```

#### Loading the dataset

There are many ways to load a dataset in R. Because the file type is a gen file, we'll use one of the packages we just downloaded and loaded. To do this, we will introduce the path where you copied the dataset to your folder.
```
GP <- read.genepop("Baseline_file.gen",ncode = 3)
```

For the following analysis, we need to prepare the dataset a little bit. This means, for example, that we need to add the population name for each individual of the dataset and make it a factor. This is not really important right now (please, ask us if you want to know exactly what the following commands are doing). Otherwise, just type:
```
# include the pop map
GP_pop <- read.table("pop_map.txt",stringsAsFactors = F)
GP_pop <- as.data.frame(GP_pop)
GP_pop_fac <- as.factor(GP_pop$V2)
table(GP_pop_fac)
```

We will also convert it to a special format to analyse genomic datasets, a "genind" objects. As a first step we just have to type:
```
# convert into genepop file
gp <- genepop_to_genind("Baseline_file.gen", n.loci=96,allele.digit.coding = 3)
gp$pop <- GP_pop_fac
```
This will convert the dataset into the format we want and store it in an object called "data". You can use whatever name you want to store the data (with certain rules, e.g. no spaces in between words).


#### First look

First, let's look how the data we just imported into R looks like. This is the first step when analysing a genomic dataset. Genomic datasets are quite large and therefore it is a bit difficult ("humanly impossible") to check each data row one by one, by hand. That's why we write scripts and use functions, so we can automate the checks and the analysis. It reduces errors (we are humans!) and assures reproducibility.

By typing the name where we store the dataset (see earlier commands),
```
GP
```
it will show us several lines of information, to summarize all the dataset. Basically we will be able to see a summary of what is inside "data". 
```
> GP
/// GENIND OBJECT /////////

 // 487 individuals; 96 loci; 192 alleles; size: 457.8 Kb

 // Basic content
   @tab:  487 x 192 matrix of allele counts
   @loc.n.all: number of alleles per locus (range: 2-2)
   @loc.fac: locus factor for the 192 columns of @tab
   @all.names: list of allele names for each locus
   @ploidy: ploidy of each individual  (range: 2-2)
   @type:  codom
   @call: read.genepop(file = "Baseline_file.gen", ncode = 3)

 // Optional content
   @pop: population of each individual (group size range: 62-215)
```
What does all this mean? The first row is telling us that the data file is a genlight object, and inside we can find 487 rows of genotypes (so 487 individuals genotyped for 96 loci, bi-allelic).
Please spend a few minutes familiarizing yourself with the information in each entry of the dataset. You can do that by using the command
```
head(GP@loc.fac)
```
The dataset we are working on consists on 487 individuals sampled from this 4 populations found in cod. 

### Calculating allele frequencies
We will again change the format of the file to calculate allele frequencies per population. You do the conversion here:
```
# convert to hierfstat
GP_hf <- genind2hierfstat(GP,pop=GP_pop_fac)
head(GP_hf[1:3,1:3])
```

And now, use the command below to calculate the allele frequencies per SNP, per population:
```
GP_stats <- basic.stats(GP_hf,diploid=TRUE,digits=4)
str(GP_stats$pop.freq)
```

Explore the dataset a little bit. You can save the table into your laptop if it is easier:
```
# Explore the dataset
#GP_stats$pop.freq
GP_stats$pop.freq[1:2]
write.table(GP_stats$pop.freq, "GP_stats")
```

Can you spot some differences already between the four populations?
For example, have a look at the following SNPs, which are called "Rhod_1_1", "HbBeta1_1", "Gm394_0364". You can access each of them by the following command:
```
GP_stats$pop.freq$Rhod_1_1
GP_stats$pop.freq$HbBeta1_1
GP_stats$pop.freq$Gm394_0364
```
Is there any population that has some of these locus fixed? Are there differences between the populations? Discuss with your partner.


#### Filtering our dataset

As a second step in a population genetics/genomic analysis, we would need to filter the SNP-data. There are many ways of filtering a dataset, and it all depends what we are interested in. One thing to remember is that whatever we choose for filtering steps, we would need to report all the steps, so other scientists can replicate our analysis and understand why the results are the way they are. Some of the parameters that genomicists filter their data on is in the % of missing data, or minor allele frequencies. For the sake of simplicity, we won't be doing any filtering process for our exercises today, for two reasons: (1) the dataset has already been filtered (see the original article for details), and (2) for simplicity (one can spend lots of time in filtering - and one should!). Of course you are very welcome to try different filters at home.

#### Hardy-Weinberg equilibrium
Next, let’s determine if our populations/loci are in Hardy-Weinberg equilibrium. Using the DartR package we will compute the χ2 statistic over the entire dataset and compute two P-values, one analytical and one derived from permutations:
First we will change the format again (to a genind object) and change a few things into the object in order the package to work.
```
gp_gi <- gi2gl(gp, verbose = NULL)
# And we need to add a metrics entry into the dataset, and convert it into a dataframe for it to work
genli_vcf<- gl.recalc.metrics(gp_gi)
genli_vcf$other$loc.metrics <- as.data.frame(genli_vcf$other$loc.metrics)
```
And now we can finally run the HWE test per locus. You can save the results into your folder if you want.
```
gg <- gl.report.hwe(genli_vcf, p = 0.05, subset = "each",bonf = TRUE)
#write.table(adeg_genind_noNA_hw, "HW_each.txt")
```

#### Principal Component Analysis (PCA)

We will perform a principal component analysis (PCA) where we will visualize the population structure at an individual level. For the principal component analysis (PCA), we will use some functions that handles genomic objects very "quickly". 
```
adeg_pca <- scaleGen(gp, NA.method="mean",scale=F)
pca.adeg_pca <- dudi.pca(adeg_pca, scale=F, nf = 10,scannf = F)
s.class(pca.adeg_pca$li, fac=gp@pop, col=funky(15),cpoint=1)
s.class(pca.adeg_pca$li,fac=gp@pop,  xax=1, yax=2,
        col=transp(funky(12),.9),
        axesel=F, cstar=0,
        cpoint=2)
```

We can save the plot with this:
```
jpeg(filename = "PCA_cod.JPEG",
     width = 600, height = 600, units="px", pointsize=12, quality=300)
s.class(pca.adeg_pca$li,fac=gp@pop,  xax=1, yax=2,
        col=transp(funky(12),.9),
        axesel=F, cstar=0,
        cpoint=2)
dev.off()
```
How do you see the distribution of populations in the plot: 
* Do baseline samples represent different populations? 
* Are all populations clustered together in the same area, or some are further appart? 
* How many clusters do you see in the baseline file? 
Discuss in pairs.

Which loci seem to be driving the differenciation in the PCA? For looking into this, we can make a loading analysis. Select the first 5 axes.
```
deg_genind <- gi2gl(gp, parallel = TRUE)
pca_genlight <- glPca(deg_genind)
```
And now the loading plot:
```
loci_axis <- loadingplot(pca_genlight,axis = 1, lab.jitter = 1)
```
If we want to save the plot into your folder, you can simply run:
```
jpeg(filename = "cod_loadingPlot.JPEG",
     width = 600, height = 600, units="px", pointsize=12, quality=300)
loadingplot(pca_genlight,axis = 1, lab.jitter = 1)
dev.off()
```
At a first glance in your opinion, is there any special loci that seem to be driving the differentiation? 

#### Degree of differentiation (Fst)

We will now have a look at the measure of genetic differentiation between the populations. We will measure this by calculating the Fst. For this, we'll make use of the package dartR. 
First, what would you expect to see? Remember that the Fst is a measure of how different a population is from another one.
When you are ready to check your expectations, you can type the following command lines to estimate Fst:
```
fst_test <- stamppFst(deg_genind, nboots = 10000, percent = 95)
fst_test$Fsts
fst_test$Pvalues
```
Are the results statistically significant? What genetic processes could have made that these populations are different from each other?

#### Population structure analysis 

For this part, we will use another software, STRUCTURE, which is a popular software to look at population structure within populations. "Its uses include inferring the presence of distinct populations, assigning individuals to populations, studying hybrid zones, identifying migrants and admixed individuals, and estimating population allele frequencies in situations where many individuals are migrants or admixed. It can be applied to most of the commonly-used genetic markers, including SNPS" (see more information here in the website: https://web.stanford.edu/group/pritchardlab/structure.html)
