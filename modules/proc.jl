######################################################
###           MODULE FOR WEIRD PROCESSING          ###
###            author: Michi Haugeneder            ###
######################################################

module proc

using PyCall, Statistics, DelimitedFiles, Dates, Random, FFTW, Dates, Profile, ProgressMeter, NCDatasets

export argparse, irncfilesindir, createwindfield, savewindasnetcdf

"Parse the arguments from 'file' and return as tuple"
function argparse(file::String)
    #read lines of arguments containing file
    argsraw = readlines(file)

    #assign VARIABLES
    sf = argsraw[3] #source folder
    wndwsz = parse(Int64, argsraw[6])
    rowleft = parse(Int64, argsraw[9])
    rowright = parse(Int64, argsraw[10])
    colup = parse(Int64, argsraw[11])
    colbottom = parse(Int64, argsraw[12])
    maxshift = [rowleft, rowright, colup, colbottom]
    nini = parse(Int64, argsraw[15]) #n for WEIRD-step1
    alpha = parse(Float64, argsraw[18]) #threshold for WEIRD-step2
    NaNthresh = parse(Float64, argsraw[21])
    return (sf, wndwsz, maxshift, nini, alpha, NaNthresh)
end

"Return vector with all irdata .nc (NetCDF4-files) in 'folder'"
function irncfilesindir(folder::String)::Vector{String}
    files = readdir(folder)
    re = r".*preproc\.nc"
    return files[occursin.(re, files)]
end

"Calculate fraction of NaNs in a given array"
function nanratio(data::Array)::Number
    nrnans = count(x->isnan(x), data)
    nrtot = count(x->true, data)
    return nrnans/nrtot
end

"Return boolean 2D-Array, if window with size 'windowsize' and array-index as lefttop
has a NaN-ratio below 'nanthresh'"
function possibletargetwindows(data::Array, windowsize::Int64, nanthresh::Number)::BitArray
    resarray = BitArray(undef, size(data, 1)-windowsize+1, size(data, 2)-windowsize+1)
    for irow in 1:(size(data, 1)-windowsize+1)
        for jcol in 1:(size(data, 2)-windowsize+1)
            nantmp = nanratio(data[irow:irow+windowsize-1, jcol:jcol+windowsize-1])
            resarray[irow, jcol] = nantmp <= nanthresh
        end
    end
    return resarray
end

"Return rows and cols (coords) for left top of reference window taking edge space
for shifting into account; search to bottom"
function limitsforsearch(windowsize::Int,
    shiftvector::Vector, rowsizedata::Int, colsizedata::Int)
    #define rowmin
    if shiftvector[1] >= 0
        rowmin = 1
    else
        rowmin = 1 - shiftvector[1]
    end
    #define rowmax
    rowmax = rowsizedata - windowsize + 1
    #define colmin
    if shiftvector[3] >= 0
        colmin = 1
    else
        colmin = 1 - shiftvector[3]
    end
    #define colmax
    if shiftvector[4] >=0
        colmax = colsizedata - (windowsize + shiftvector[4]) + 1
    else
        colmax = colsizedata - windowsize + 1
    end
    #quality check: all values should be >0 and xxmax > xxmin! if false, then bad
    qualitycheck = 0 < rowmin < rowmax && 0 < colmin < colmax
    if !qualitycheck
        return [0 0 0 0]
        println("ERROR: Windowsize is too small!")
    else
       # println(rowmin, rowmax, colmin, colmax)
        return [rowmin, rowmax, colmin, colmax]
    end
end

"For the WEIRD algorithm create an array with a random order of positions to visit
col1:row, col2:col"
function createrandomvisitorder(wndwsize::Int64)
    temppos = zeros(Int, wndwsize*wndwsize, 2)
    pos = similar(temppos)
    for i in 0:wndwsize-1
        for j in 0:wndwsize-1
            temppos[i*wndwsize+j+1,1] = j
            temppos[i*wndwsize+j+1,2] = i
        end
    end
    perm = shuffle(collect(1:size(temppos, 1)))
    for i in 1:size(temppos, 1)
        pos[i,:] = temppos[perm[i],:]
    end
    return pos
end

"Precalculations for the WEIRD-algorithm. Return array containing maxshift1, maxshift2 per row and col"
function precalc(sizeIR1, reflefttoprow::Vector, reflefttopcol::Vector,
    wndwsize, maxshiftin1, maxshiftin2)::Array
    ######################      PRECALCULATIONS     ########################
    #determine distance to bottom as maxshift-bottom if necessary;
    #else keep original maxshift
    retref = zeros(Int, size(reflefttoprow, 1), size(reflefttopcol,1), 2)
    for rj in 1:size(reflefttopcol,1)
        for ri in 1:size(reflefttoprow,1)
            cond1 = sizeIR1- reflefttoprow[ri] + 1 - wndwsize
            if cond1 < maxshiftin2
                retref[ri, rj, 2] = cond1
                #set vertical shift to be symmetric (rowshift_pos = -rowshift_neg)
                retref[ri, rj, 1] = - retref[ri, rj, 2]
            else
                retref[ri, rj, 1] = maxshiftin1
                retref[ri, rj, 2] = maxshiftin2
            end
        end
    end
    return retref
end

"WEIRD-algorithm step 1"
function WEIRDstep1!(accdiffs::Array, nrreps::Array, IRdatain::Array, maxshiftinput::Vector, nrtargetwndwsrow::Integer,
    nrtargetwndwscol::Integer, reflefttopin::Vector, nrsamplesin::Int64, samplelocationsin)
    #######################        STEP 1        ###########################
    #create Array with number of inital differences in step1
    nrinit = zeros(Int64, maxshiftinput[2]-maxshiftinput[1]+1,
    maxshiftinput[4]-maxshiftinput[3]+1)

    #loop over all target windows
    #cols
    for i4 in maxshiftinput[3]:maxshiftinput[4]
        #rows
        for j4 in maxshiftinput[1]:maxshiftinput[2]
            #loop over sample locations
            nrsummands=0
            k4=1
            while nrsummands<=nrsamplesin
                refrow = reflefttopin[1]+samplelocationsin[k4,1]
                refcol = reflefttopin[2]+samplelocationsin[k4,2]
                toadd = abs(IRdatain[refrow,refcol,1] - IRdatain[refrow+j4,refcol+i4,2])
                if !isnan(toadd)
                    accdiffs[j4-maxshiftinput[1]+1, i4-maxshiftinput[3]+1] += toadd
                    nrreps[j4-maxshiftinput[1]+1, i4-maxshiftinput[3]+1] += 1
                    nrsummands += 1
                end
                k4 += 1                
            end
            nrinit[j4-maxshiftinput[1]+1, i4-maxshiftinput[3]+1] = nrsummands
        end
    end
    return nrinit
end

"WEIRD-algorithm step2"
function WEIRDstep2(accumulateddiffs::Array, wndwsize::Int64, nrsamples::Int64, alpha)::Float64
    #find minimum of sample point differences
    TSmin = minimum(accumulateddiffs)
    #estimate minimum if all points would be evaluated
    TMmin = TSmin*(wndwsize^2/nrsamples)
    #calculate critical value with the help of alpha
    TMcr = alpha*TMmin
    return TMcr
end

"WEIRD-algorithm step 3, operate on accumulateddiffs & nrrepitations"
function WEIRDstep3!(accumulateddiffs::Array, nrrepitations::Array, nrinit::Array,
    IRdata::Array, maxshift::Vector, reflefttop::Vector, nrsamples::Int64, samplelocations, TMcr)
    szsamplelocs = size(samplelocations, 1) # nr of pixel per targetwindow
    nrremainingwindows = (maxshift[4]-maxshift[3]+1)*(maxshift[2]-maxshift[1]+1) #nr targetwindows not exceeding TMr
    nrinit .+= 1 #next step after WEIRDstep1
    #array containing false, if corresponding target window reached threshold
    posarray = BitArray(undef, maxshift[2]-maxshift[1]+1, maxshift[4]-maxshift[3]+1)
    posarray .= true
    A = CartesianIndices(posarray)
    while nrremainingwindows>1
        for a in A[posarray]
            colshift = last(Tuple(a)) + maxshift[3] - 1
            rowshift = first(Tuple(a)) + maxshift[1] - 1
            refrow = reflefttop[1] + samplelocations[nrinit[a], 1]
            refcol = reflefttop[2] + samplelocations[nrinit[a], 2]
            toadd = abs(IRdata[refrow, refcol, 1] - IRdata[refrow+rowshift, refcol+colshift, 2])
            nrinit[a] += 1
            if !isnan(toadd)
                accumulateddiffs[a] += toadd
                nrrepitations[a] += 1
            else
                if nrinit[a] <= szsamplelocs
                    nanbreak = false
                else
                    nanbreak = true
                end
                while !nanbreak
                    refrow = reflefttop[1] + samplelocations[nrinit[a], 1]
                    refcol = reflefttop[2] + samplelocations[nrinit[a], 2]
                    toadd = abs(IRdata[refrow, refcol, 1] - IRdata[refrow+rowshift, refcol+colshift, 2])
                    nrinit[a] += 1
                    if !isnan(toadd)
                        accumulateddiffs[a] += toadd
                        nrrepitations[a] += 1
                        nanbreak = true
                    end
                    #stop WEIRDstep3! after this repitition for all target windows,
                    #since one window is completely done (too many NaN-differences)
                    if nrinit[a] > szsamplelocs
                        nanbreak = true
                    end
                end
            end
            if accumulateddiffs[a] > TMcr
                posarray[a] = false
                nrremainingwindows -= 1
            end
            #stop WEIRDstep3! after this repitition for all target windows,
            #since one window is completely done (too many NaN-differences)
            if nrinit[a] > szsamplelocs
                nrremainingwindows = -2
                @info("reached limit for one corr-search", rowshift, colshift)
            end
        end #for eachindex
        if nrremainingwindows == 1 #unique solution using nr. repitions
            nrrepitations[posarray] .+= 1 #increase repitations to define clear solution
            nrremainingwindows = -2 #set parameter to break while-loop
        end
    end #while
end #function

"Postprocessing for WEIRD-algorithm"
function WEIRDpostproc(accumulateddiffs::Array, nrrepitations::Array,
    maxshiftin::Vector)
    #for case 'below' snow surface
    if maxshiftin[1] == maxshiftin[2] == -999
        return NaN, NaN, NaN
    end
    #######################       POST-PROCESSING       #######################
    #if multiple matching windows take one with lowest accumulateddiff
    maxreps = maximum(nrrepitations)
    #factor for quality of match (how clear is the match?)
    if maxreps == 0
        matchquality = -1.5
    else
        matchquality = float(maxreps)/mean(nrrepitations)
    end
    posmaximum = findall(x->x==maxreps, nrrepitations)
    #if multiple maxima: find minimum accumulateddiff for windspeed and -1 for matchquality
    if size(posmaximum, 1) > 1
            matchquality = - size(posmaximum, 1)
            idxmultmatch = findall(x->x==maxreps, nrrepitations)
            minimumacc = minimum(accumulateddiffs[idxmultmatch])
            posmaximum = findall(x->x==minimumacc, accumulateddiffs)
            #if still >1 maxima
            if size(posmaximum, 1) > 1
                matchquality = -0.5
            end
    end
    return posmaximum[1][1] + maxshiftin[1] - 1, posmaximum[1][2] + maxshiftin[3] - 1, matchquality
end

"Print quality control output to screen"
function qualcontr(windfield::Array, maxqual::Array)
    println("######### QUALITY CONTROL ##########")
    #nr points windspeed==0
    windis0 = count(x->x==0, windfield)
    #total nur notnan windspeed
    windnotnan = count(x->!isnan(x), windfield)
    #nr points greater 0 (single maximum)
    nrptsg0 = count(x->x>0, maxqual)
    #nr points smaller 0 (multiple maxima)
    nrptss0 = count(x->x<0, maxqual)
    #all points !=0
    allpts = nrptsg0 + nrptss0
    #multiple maxima
    doubmax = count(x->x==-2, maxqual)
    tripmax = count(x->x==-3, maxqual)
    quatrmax = count(x->x==-4, maxqual)
    maxmax = -minimum(maxqual)
    notunique = count(x->x==-0.5, maxqual)
    println("ratio windspeed=0: ", round(windis0/windnotnan*100.0, digits=2), "%")
    println("unique maxima: ", nrptsg0, " -> ", round(nrptsg0*100/(allpts), digits=2), "%")
    println("multiple maxima: ", nrptss0, " -> ", round(nrptss0*100/(allpts), digits=2), "%")
    println("double maxima: ", doubmax, " -> ", round(doubmax*100/(allpts), digits=2), "%")
    println("triple maxima: ", tripmax, " -> ", round(tripmax*100/(allpts), digits=2), "%")
    println("quatruple maxima: ", quatrmax, " -> ", round(quatrmax*100/(allpts), digits=2), "%")
    println("maximal: ", maxmax, " times maximum")
    println("not unique: ", notunique, " -> ", round(notunique*100/(allpts), digits=2), "%")
    println("####################################")
end

"Perform WEIRD algorithm; Flag 'qualitycheck' returns additional values -> slower"
function WEIRD(IRdata::Array, samplelocations::Array, reflefttop::Vector,
    wndwsize::Int64, maxshiftin::Vector, nrtargetrow::Integer, nrtargetcol::Integer,
    accdiffs::Array, nrreps::Array, nrsamples::Int64, alpha::Float64)
    #check if 'below' snow surface
    if maxshiftin[2] < 0
        return WEIRDpostproc(accdiffs, nrreps, [-999, -999, -999, -999])
    end
    #reset from previous step
    fill!(accdiffs, 0.0)
    fill!(nrreps, 0)
    #step1
    nrinit = WEIRDstep1!(accdiffs, nrreps, IRdata,
    maxshiftin, nrtargetrow, nrtargetcol, reflefttop, nrsamples, samplelocations)
    #step2
    TMcr = WEIRDstep2(accdiffs, wndwsize, nrsamples, alpha)
    #step3
    WEIRDstep3!(accdiffs, nrreps, nrinit, IRdata, maxshiftin, reflefttop, nrsamples, samplelocations, TMcr)
    return WEIRDpostproc(accdiffs, nrreps, maxshiftin)
end

"Unify functions to create a windfield estimation"
function createwindfield(IRdata::Array, NaNthresh::Number, windowsize::Int,
    pics::Vector, maxshiftvector::Vector, nrsamples::Int, alpha::Float64)
    #correlation search: define indices for reference window left top (min/max rows/cols)
    limitssearch = limitsforsearch(windowsize, maxshiftvector,
    size(IRdata, 1), size(IRdata, 2))
    #create vectors containing all allowed rows/cols for correlation search
    rowsearch = collect(limitssearch[1]:limitssearch[2])
    colsearch = collect(limitssearch[3]:limitssearch[4])

    samplelocations = createrandomvisitorder(windowsize)

    #initialize array later containing the windfield with same size as IRdata
    windfield = fill(NaN, size(IRdata, 1), size(IRdata, 2), size(pics, 1) - 1, 2)

    #initialize array with nr. of maxima from WEIRD-method per position
    maxqual = zeros(Float64, size(IRdata, 1), size(IRdata, 2), size(pics, 1) - 1)

    #precalculation to obtain maximum allowed vertical (row) shift (at the bottom) per column
    maxshiftarray = precalc(size(IRdata, 1), rowsearch ,
    colsearch, windowsize, maxshiftvector[1], maxshiftvector[2])

    #create BitArray for NaN-threshold-check (e.g. snowsurface or between screens)
    NaNthreshcheck = possibletargetwindows(IRdata[:,:,1], windowsize, NaNthresh)

    #initialize arrays with maximum dimensions for each thread
    maxrowdim = (2*maximum(maxshiftarray[:,:,2])+1)
    maxcoldim = (maxshiftvector[4]-maxshiftvector[3]+1)
    accdiffs = zeros(Float64, maxrowdim, maxcoldim)
    nrreps = zeros(Int64, maxrowdim, maxcoldim)

    #loop over all allowed left top positions for the correlation search
    println("Filling windfield-array...")
    #@showprogress for i2 in 1:size(colsearch, 1)
    for i2 in 1:size(colsearch, 1)
        if mod(i2,50) == 0
            println(string(round(Int, i2*100/size(colsearch, 1)), "%"))
        end
        nrtargetwndwscol = (maxshiftvector[4]-maxshiftvector[3]+1)
        #print progress
        #println(i2, "/", size(colsearch, 1))
        for j2 in 1:size(rowsearch, 1)
            nrtargetwndwsrow = (2*maxshiftarray[j2,i2,2]+1)
            #check for NaN-threshold (e.g. snowsurface, between screens)
            #if one target window is not possible -> skip search for this interrogation window
            any(x->!x, NaNthreshcheck[(maxshiftarray[j2,i2,1]:maxshiftarray[j2,i2,2]).+rowsearch[j2], (maxshiftvector[3]:maxshiftvector[4]).+colsearch[i2]]) && continue
            #filling the windfield array
            for l2 in 1:size(pics, 1)-1
                (windfield[j2+limitssearch[1]-1,i2+limitssearch[3]-1,l2,1],
                windfield[j2+limitssearch[1]-1,i2+limitssearch[3]-1,l2,2],
                maxqual[j2+limitssearch[1]-1,i2+limitssearch[3]-1,l2]) =
                WEIRD(IRdata[:,:,pics[l2]:pics[l2+1]],samplelocations,
                [rowsearch[j2],colsearch[i2]],windowsize,
                [maxshiftarray[j2, i2, 1], maxshiftarray[j2, i2, 2], maxshiftvector[3], maxshiftvector[4]],
                nrtargetwndwsrow, nrtargetwndwscol,
                accdiffs[:,:], nrreps[:,:], nrsamples, alpha)
            end
        end
    end
    #printing some quality control output
    qualcontr(windfield, maxqual)
    println("Done")
    return windfield, maxqual
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