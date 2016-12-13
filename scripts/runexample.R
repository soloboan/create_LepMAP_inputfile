# Creating a LepMAP input file by generating full-sib families from pedigree information

#------------------------------------------------------------------------------------------------#
# The createLepMapinpute.R/createLepMapinpute_win3264.R scripts requires the following parameters and information to run sucessfully

#1. pedfaminfo=""     +++++ The fam file of PLINK with pedigree information or pedigree file with the column (IID, IID, Sire, Dam) 
#2. genoplinkped=""   +++++ Linkage file format genotype data (eg. PLINK ped file, 6 columns of fam info and 7-end of genotypes coded in allele)
#3. outname=""        +++++ Output filename
#-------------------------------------------------------------------------------------------------#

# running the script
source('scripts/createLepMapinpute_win3264.R')
createLepMapinpute(pedfaminfo='Ex2_10K.fam',genoplinkped='Ex2_10K.ped',
                   outname='Ex2LepMap')

#directly running on the command line in linux
./createLepMapinpute.R Examplefiles/Ex2_10K.fam Examplefiles/Ex2_10K.ped Examplefiles/Ex2LepMap

#If the above gives the error ‘R interpreter parse error’ then please find where your R is located on the server eg. My R is loca$
opt/bin/Rscript createLepMapinpute.R Examplefiles/Ex2_10K.fam Examplefiles/Ex2_10K.ped Examplefiles/Ex2LepMap
