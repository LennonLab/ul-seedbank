################################################################################
#                                                                              #
# University Lake Respiration                                                  #
#                                                                              #
################################################################################
#                                                                              #
# Written by: Mario Muscarella                                                 #
#                                                                              #
# Last update: 2015/12/03                                                      #
#                                                                              #
################################################################################
#                                                                              #
# Notes: This code provides the initial analysis for the University Lake       #
#        Bacterial respiration analysis. The code imports all of the raw       #
#        PreSens data and uses linear regression to calculate respiration      #
#        rate during the first 5 quality hours of the experiments. At the      #
#        moment the user has to define the start time for the calcuation       #
#                                                                              #
# Dependencies:                                                                #
#         1. PreSensRespiration.R                                              #
#         2. PreSensInteractiveRegression.R                                    #
#                                                                              #
# Issues: It seems inefficient to do each analysis by itself                   #
#                                                                              #
# Recent Changes:                                                              #
#                                                                              #
# Future Changes (To-Do List):                                                 #
#         1.                                                                   #
#                                                                              #
################################################################################

# Initial Setup
rm(list=ls())
getwd()
setwd("~/GitHub/UniversityLake/analyses")

# Load PreSens Package
# Inport the function from source file
source("../bin/PreSensInteractiveRegression.r") # Use to pick the time windows
source("../bin/PreSensRespiration.R")

# UL Respiration Analyses
# UL001
input  <-  "../data/Respiration/Output/20130419_UL001_Oxygen.txt"
output <-  "../data/Respiration/Rates/20130419_UL001_BR.txt"
in.name <- c(rep("Filtered", 4), rep("Empty", 2),
             rep("Empty", 6),
             rep("Unfiltered", 4), "Empty", "Blank",
             rep("Empty", 5), "Blank")
PreSens.Respiration2(infile = input, outfile = output, start = 1,
                                 end = 6, name.in = in.name)

# UL002
input  <-  "../data/Respiration/Output/20130425_UL002_Oxygen.txt"
output <-  "../data/Respiration/Rates/20130425_UL002_BR.txt"
in.name <- c(rep("Filtered", 4), rep("Empty", 2),
             rep("Empty", 6),
             rep("Unfiltered", 4), rep("Empty", 2),
             rep("Empty", 4), rep("Blank", 2))
PreSens.Respiration2(infile = input, outfile = output, start = 1.5,
                                 end = 6.5, name.in = in.name)

# UL003
input  <-  "../data/Respiration/Output/20130502_UL003_Oxygen.txt"
output <-  "../data/Respiration/Rates/20130502_UL003_BR.txt"
in.name <- c(rep("Filtered", 4), rep("Empty", 2),
             rep("Empty", 6),
             rep("Unfiltered", 4), "Empty", "Blank",
             rep("Empty", 5), "Blank")
PreSens.Respiration2(infile = input, outfile = output, start = 3.5,
                                 end = 8.5, name.in = in.name)

# UL004
input  <-  "../data/Respiration/Output/20130509_UL004_Oxygen.txt"
output <-  "../data/Respiration/Rates/20130509_UL004_BR.txt"
in.name <- c(rep("Filtered", 4), rep("Empty", 2),
             rep("Empty", 6),
             rep("Unfiltered", 4), "Empty", "Blank",
             rep("Empty", 5), "Blank")
PreSens.Respiration2(infile = input, outfile = output, start = 2,
                                 end = 7, name.in = in.name)

# UL005
input  <-  "../data/Respiration/Output/20130517_UL005_Oxygen.txt"
output <-  "../data/Respiration/Rates/20130517_UL005_BR.txt"
in.name <- c(rep("Filtered", 4), rep("Empty", 2),
             rep("Unfiltered", 4), rep("Empty", 2),
             rep("Empty", 4), rep("Blank",2),
             rep("Empty", 6))
PreSens.Respiration2(infile = input, outfile = output, start = 5,
                     end = 10, name.in = in.name)

# UL007
input  <-  "../data/Respiration/Output/20130529_UL007_Oxygen.txt"
output <-  "../data/Respiration/Rates/20130529_UL007_BR.txt"
in.name <- c(rep("Unfiltered", 4), rep("Empty", 2),
             rep("Filtered", 4), rep("Empty", 2),
             rep("Blank",2), rep("Empty", 4),
             rep("Empty", 6))
# These data are odd and thus I'm skipping

# UL008
input  <-  "../data/Respiration/Output/20130607_UL008_Oxygen.txt"
output <-  "../data/Respiration/Rates/20130607_UL008_BR.txt"
in.name <- c(rep("Unfiltered", 4), rep("Empty", 2),
             rep("Filtered", 4), rep("Empty", 2),
             rep("Blank",2), rep("Empty", 4),
             rep("Empty", 6))
PreSens.Respiration2(infile = input, outfile = output, start = 1.5,
                     end = 6.5, name.in = in.name)

# UL009
input  <-  "../data/Respiration/Output/20130614_UL009_Oxygen.txt"
output <-  "../data/Respiration/Rates/20130614_UL009_BR.txt"
in.name <- c(rep("Unfiltered", 4), rep("Empty", 2),
             rep("Filtered", 4), rep("Empty", 2),
             rep("Blank",2), rep("Empty", 4),
             rep("Empty", 6))
PreSens.Respiration2(infile = input, outfile = output, start = 4,
                     end = 9, name.in = in.name)

# UL010
input  <-  "../data/Respiration/Output/20130621_UL010_Oxygen.txt"
output <-  "../data/Respiration/Rates/20130621_UL010_BR.txt"
in.name <- c(rep("Unfiltered", 4), rep("Empty", 2),
             rep("Filtered", 4), rep("Empty", 2),
             rep("Blank",2), rep("Empty", 4),
             rep("Empty", 6))
PreSens.Respiration2(infile = input, outfile = output, start = 5,
                     end = 10, name.in = in.name)

# UL011
input  <-  "../data/Respiration/Output/20130628_UL011_Oxygen.txt"
output <-  "../data/Respiration/Rates/20130628_UL011_BR.txt"
in.name <- c(rep("Unfiltered", 4), rep("Empty", 2),
             rep("Filtered", 4), rep("Empty", 2),
             rep("Blank",2), rep("Empty", 4),
             rep("Empty", 6))

# These data are odd and thus I'm skipping

# UL012
input  <-  "../data/Respiration/Output/20130621_UL012_Oxygen.txt"
output <-  "../data/Respiration/Rates/20130621_UL012_BR.txt"
in.name <- c(rep("Filtered", 4), rep("Empty", 2),
             rep("Unfiltered", 4), rep("Empty", 2),
             rep("Blank",2), rep("Empty", 4),
             rep("Empty", 6))
PreSens.Respiration2(infile = input, outfile = output, start = 5,
                     end = 10, name.in = in.name)


# Test Code Using Interactive Form of Analysis
PreSens.Respiration(input, "../data/Respiration/Rates/test.txt")
