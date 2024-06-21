######################################################
###           MODULE FOR STICHING OUTPUT           ###
###            author: Michi Haugeneder            ###
######################################################
module sti

using ProgressMeter, NCDatasets
export argparse, windncfilesindir, sortstartendandfilevec, checkgaps,
concentatewindfields

"Parse the arguments from 'file' and return as tuple"
function argparse(file::String)
    #read lines of arguments containing file
    argsraw = readlines(file)
    #assign VARIABLES
    directory = argsraw[4]
    filestam = argsraw[7]
    defllvl = parse(Int64, argsraw[10])
    return (directory, filestam, defllvl)
end

"Return vector with all wind .nc (NetCDF4-files) in 'folder'
containing 'prefix'"
function windncfilesindir(folder::String, prefix::String)::Vector{String}
    files = readdir(folder)
    re = Regex("$(prefix).*nc")
    return files[occursin.(re, files)]
end

"Sort 'startvec', 'endvec', and 'filevec'"
function sortstartendandfilevec(filevec::Vector)
    startframevec = zeros(Int64, length(filevec))
    endframevec = zeros(Int64, length(filevec))
    for i in eachindex(filevec)
        tmp = match(r"\d{1,4}to\d{1,4}.nc", filevec[i]).match
        startframevec[i] = parse(Int64,match(r"\d{1,4}to", tmp).match[1:end-2])
        endframevec[i] = parse(Int64, match(r"\d{1,4}\.nc", tmp).match[1:end-3])
    end
    #get permutation for sorting in increasing order
    p = sortperm(startframevec)
    return (startframevec[p], endframevec[p], filevec[p])
end

"Check for no gaps in frame data"
function checkgaps(startsorted::Vector, endsorted::Vector)::Bool
    gapcheckpassed = true
    for igapcheck in eachindex(startsorted)
        if igapcheck != 1 && endsorted[igapcheck-1] != startsorted[igapcheck]
        @warn("There is a gap!", startsorted[igapcheck], igapcheck)
        gapcheckpassed = false
        end
    end
    if gapcheckpassed
        @info("No frame gap in the data. Complete!", startsorted[1], endsorted[end])
    end
    return gapcheckpassed
end
"Concentate all dual frame files"
function concentatewindfields(directory::String, windsorted::Vector)::Array
    innerwindfile = NCDataset(joinpath(directory, windsorted[1]), "r")
    winddat = innerwindfile["wind"]
    windtot = winddat[:,:,:]
    close(innerwindfile)
    @showprogress for idx in 2:size(windsorted,1)
        innerwindfile = NCDataset(joinpath(directory, windsorted[idx]), "r")
        winddat = innerwindfile["wind"]
        windtmp = winddat[:,:,:]
        close(innerwindfile)
        windtot = cat(windtot, windtmp, dims=3)
    end
    return windtot
end

"Saving windfield-sniplet as NetCDF4"
function savewindasnetcdf(data::Array, name::String, target::String, deflatelvl::Int64=5)
    @info("Saving windfield-sniplet output to NetCDF4-file")
    ds = NCDataset(target, "c")
    defDim(ds, "row", size(data, 1))
    defDim(ds, "col", size(data, 2))
    defDim(ds, "frame", size(data, 3))
    defDim(ds, "u, w", size(data, 4))
    defVar(ds, name, data, ("row", "col", "frame", "u, w"); shuffle=true, deflatelevel=deflatelvl)
    close(ds)
end
end #module