######################################################
###           WEIRD WIND FIELD ESTIMATION          ###
###            author: Michi Haugeneder            ###
######################################################
#=
Obtaining a 2D-windfield from the IR-screen data using WEIRD
=#

using NCDatasets, Dates
#define 'importdir' as directory containing source code, include modules
importdir = @__DIR__
include(joinpath(importdir, "modules", "proc.jl"))
import .proc
#path to file containing arguments
argsfile = joinpath(@__DIR__, "args", "proc.txt")

######################################################
###                ARGUMENT PARSING                ###
######################################################
#Parse the arguments from file
(sourcefolder, windowsize, maxshift, nini, alpha, NaNthresh) =
proc.argparse(argsfile)

######################################################
###                  PREPARATION                   ###
######################################################
#read directory and look for NetCDF4-files
irncfiles = proc.irncfilesindir(sourcefolder)
if length(irncfiles)>1
    @warn("Processing not implemented for more than one NetCDF4-File. Exiting.")
    exit()
end

#note: adapt the following lines for the computer
#you are using (splitting pairs of frames to different tasks)
println("####################################")
println("####### GENERAL INFORMATION ########")
@show(gethostname())
@show(Threads.nthreads())
taskid = parse(Int64, ENV["SLURM_ARRAY_TASK_ID"])
@show taskid
@show Dates.now()
println("####################################")
println()
println()
println("####################################")
println("############## WEIRD ###############")

######################################################
###               WINDFIELD ESTIMATION             ###
######################################################
outfolder = string(joinpath(sourcefolder, "WEIRD_output"),"/")
mkpath(outfolder)

irfile = NCDataset(joinpath(sourcefolder, irncfiles[1]), "r")
irinput = irfile["irdata"]
#starting at frame 351 to have a complete backlooking average
framestart = 250+taskid
frameend = framestart+1
irin = irinput[:,:,framestart:frameend]
close(irfile)

@show framestart
@show frameend

#WEIRD processing
(windfield, maxqual) = proc.createwindfield(irin, NaNthresh,
windowsize, collect(1:size(irin, 3)), maxshift, nini, alpha)

#save output
proc.savewindasnetcdf(windfield, "wind", string(outfolder, "wind_", irncfiles[1][1:end-3], "_fr", framestart, "to", frameend, ".nc"), 8)
