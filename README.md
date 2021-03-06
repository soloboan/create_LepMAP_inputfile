# create an input file format for LepMAP
A script for creating the input file format (Linkage format with ONLY full sib families) required by LepMap program.

Note:: The script can handle multi-generation data with overlapping generations.

#### There are two script, one for running directly in R (both windows and linux) and the other for running it as a command line script  
        - createLepMapinpute_win3264.R     --- running directly in R  
        - createLepMapinpute.R             --- running on linux as command line script  

#### The example file can be untarred and used. tar -xzvf Examples.tar.gz

#### The script requires 4 Arguments  
       - 1. pedfaminfo     == The fam file from PLINK. This can be generated with –-make-bed and should contain the pedigree information        - 2. genoplinkped   == The genotype data in PLINK ped format, can be generated with –-recode  
       - 3. minfamsize     == The minum family size 
       - 4. outname        == The output file name (only prefix)

### how to run the “createLepMapinpute_win3264.R ” file in R
open R either in windows or linux and then source the script and use the function to generate the LepMap file as below  
         
        - source('createLepMapinpute_win3264.R')  
        - createLepMapinpute(pedfaminfo='Ex2_10K.fam',genoplinkped='Ex2_10K.ped',minfamsize=7,outname='Ex2LepMap')  

### how to run the “createLepMapinpute_win3264.R ” directly on the command line in linux – normal script execution  

        - ./createLepMapinpute.R Ex2_10K.fam Ex2_10K.ped 7 Ex2LepMap  

If the above gives the error ‘R interpreter parse error’ 
then please find where your R is located on the server eg. My R is located at “opt/bin/Rscript”, thus I run the script as  

       - opt/bin/Rscript createLepMapinpute.R Ex2_10K.fam Ex2_10K.ped 7 Ex2LepMap

The script generates two output files  
         
       - 'prefix.oldids'   -- This is the recoded Ids and corresponding oldids (eg. Ex2LepMap.oldids)  
       - 'prefix.linkage'  -- The LepMap file format (eg. Ex2LepMap. linkage)  

Note that if the genotype file and the total number of animal is large, then it take a while to do this in R  

Benchmarking :: On a window 64bit - 3.60GHz – intel i7, 16GB ram and an SSD harddisk, it takes about 20mins to run dataset of 3000 animal and 10,000 markers. But you can always split dataset into chunks (eg. Chromosome) and supply to the file to the script.  

### Example files
 
     - Files ==> “Ex2_10K.fam + Ex2_10K.ped “ -
 A 10,000 marker panel with 2,500 genotype animals comparing two (2) generations, 
the mating ratio was 1:4 in each generation. Thus the last generation offspring are have their grandparent genotyped.  
After running the script 3050 animal are generated with genotypes.  

### Enjoy LepMaping 
