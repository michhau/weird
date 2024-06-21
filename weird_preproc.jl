######################################################
###        PREPROCESSING SCRIPT FOR WEIRD          ###
###            author: Michi Haugeneder            ###
######################################################
#=
To be run on local machine with access to the IR-data
stored as .NetCDF4-files. Preparation step.
'Raw' IR input data must be cut before!

- read in single .NetCDF4-files
- perform preprocessing
- store output as NetCDF4
=#
using NCDatasets
#define 'importdir' as directory containing source code, include modules
importdir = @__DIR__
include(joinpath(importdir, "modules", "preproc.jl"))
import .prep
#path to file containing arguments
argsfile = joinpath(@__DIR__, "args", "preproc.txt")

######################################################
###                ARGUMENT PARSING                ###
######################################################
#Parse the arguments from file
(sourcefile, fps, gaussparam, backlookingavgtime, lp_threshold, outfolder) =
prep.argparse(argsfile)

######################################################
###             PREPROCESSING AND SAVING           ###
######################################################
prep.preprocessing(sourcefile, gaussparam, fps,
backlookingavgtime, lp_threshold, outfolder)
