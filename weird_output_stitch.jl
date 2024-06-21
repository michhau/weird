######################################################
###        STICHING WINDFIELD OUTPUT FILES         ###
###            author: Michi Haugeneder            ###
######################################################
#=
Stiching the single NetCDF4 output files to one big file
files are assumed to be named $prefix$nrstartframe[to]$nrendframe[.nc],
for example 'wind_irdata_per3_fr351to352.nc'
=#
using NCDatasets

#define 'importdir' as directory containing source code, include modules
importdir = @__DIR__
include(joinpath(importdir, "modules", "stitch.jl"))
import .sti
#path to file containing arguments
argsfile = joinpath(@__DIR__, "args", "stitch.txt")

######################################################
###                ARGUMENT PARSING                ###
######################################################
#Parse the arguments from file
(directory, prefix, defllvl) = sti.argparse(argsfile)

#read directory and look for windfield NetCDF4-files
windncfiles = sti.windncfilesindir(directory, prefix)

######################################################
###                  FILE STITCHING                ###
######################################################
#create (if not already existing) output folder
outfolder = string(joinpath(directory, "stitched"),"/")
mkpath(outfolder)

#extract start frame from file names
(startfr, endfr, windncsorted) = sti.sortstartendandfilevec(
    windncfiles)

#checkforgaps
nogaps = sti.checkgaps(startfr, endfr)

windtot = sti.concentatewindfields(directory, windncsorted)

target = string(outfolder, prefix, startfr[1], "to", endfr[end], ".nc")

sti.savewindasnetcdf(windtot, "wind", target, defllvl)
