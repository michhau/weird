######################################################
###          MODULE FOR WEIRD PREPROCESSING        ###
###            author: Michi Haugeneder            ###
######################################################
module prep

using ReadWriteDlm2, MAT, DelimitedFiles, PyCall,
ImageFiltering, ProgressMeter, NCDatasets, LsqFit,
FFTW, Statistics, Dates

export argparse, preprocessing

######################################################
###                       I/O                      ###
######################################################
"Parse the arguments from 'file' and return as tuple"
function argparse(file::String)
    #read lines of arguments containing file
    argsraw = readlines(file)
    #assign VARIABLES
    sf = argsraw[4] #source file
    fps = parse(Float64, argsraw[7])
    gaussrow = parse(Int64, argsraw[10])
    gausscol = parse(Int64, argsraw[11])
    gp = [gaussrow, gausscol]   #gaussparameters
    hpt = parse(Float64, argsraw[14]) #high-pass
    lpf = parse(Float64, argsraw[17]) #low-pass
    outfldr = argsraw[20]
    #variables for task splitting
    #firstframe = parse(Int64, argsraw[24])
    #lastframe = parse(Int64, argsraw[27])
    #nrtasks = parse(Int64, argsraw[30])
    return (sf, fps, gp, hpt, lpf, outfldr)
end

"Load variable from NetCDF4 file"
function loadfromNetCDF4(file::String, varname::String)
    netfile = NCDataset(file, "r")
    dat = netfile[varname]
    dattot = dat[:,:,:]
    close(netfile)
    return dattot
end

"Saving 3D IRdata-array as NetCDF4"
function saveirasnetcdf(data::Array, name::String, target::String, deflatelvl::Int64=5)
    @info("Saving IRdata output to NetCDF4-file")
    ds = NCDataset(target, "c")
    defDim(ds, "row", size(data, 1))
    defDim(ds, "col", size(data, 2))
    defDim(ds, "frame", size(data, 3))
    defVar(ds, name, data, ("row", "col", "frame"); shuffle=true, deflatelevel=deflatelvl)
    close(ds)
end

######################################################
###                     FILTERS                    ###
######################################################
"Only back-looking moving average."
function backlookingmovavg(X::Vector, numofele::Int)::Vector
    len = size(X, 1)
    #create vector with vec_isnan[i]=1 if X[i]=NaN, 0 otherwise
    vec_isnan = isnan.(X)
    if numofele >= len
        println("#avg elements >= size vector! Returning original vector!")
        return X
    else
        Y = fill(NaN, size(X, 1))
        summed = 0
        curr_nans = 0
        n = 1
        while n <= numofele    #first part until numofele elements available
            #if not NaN
            if !vec_isnan[n]
                summed += X[n]
                Y[n] = summed/(n - curr_nans)
            #if NaN
            else
                curr_nans += 1
                if n>1
                    Y[n] = Y[n-1]
                else
                    Y[n] = NaN
                end
            end
            n += 1
        end
        while n <= len    #second part until end (length of input vector)
            curr_nans += vec_isnan[n]
            curr_nans -= vec_isnan[n-numofele]
            if !vec_isnan[n]
                summed += X[n]
            end
            if !vec_isnan[n-numofele]
                summed -= X[n-numofele]
            end
            if numofele>curr_nans
                Y[n] = summed/(numofele-curr_nans)
            else
                Y[n] = NaN
            end
            #println(curr_nans)
            n += 1
        end
    end
    return Y
end

"Perform back-looking time-averaging for all points as due to 'Inagaki et al. (2013)'
'nrelements' refers to third dimension of 'IRdata' (time-dimension)"
function subtractbacklookingtimeaverage(IRdatain::Array, nrelements::Int)::Array
    IRdatareturn = similar(IRdatain)
    @info("Performing back-looking time averaging...")
    @sync for i in 1:size(IRdatain, 2)
        Threads.@spawn for j in 1:size(IRdatain, 1)
            IRdatareturn[j,i,:] = IRdatain[j,i,:] - backlookingmovavg(IRdatain[j,i,:], nrelements)
        end
    end
    println("Done")
    return IRdatareturn
end

"Set the pixels showing the snow surface (see condition) to NaN
and return row nr. for each column"
function setsnowsurfacetonan(irdata::Array, thresh::Float64=1.5)
    @info("Setting snow-surface to NaN")
    irout = irdata
    snowsurf = zeros(Int64, size(irdata, 2))
    #only do it for the lowest 30% assuming top is row=1
    lim = round(Int, size(irdata, 1)*0.7)
    @sync for icol in 1:size(irdata, 2)
        Threads.@spawn for jrow in lim:size(irdata, 1)
            if any(x->x<thresh,irdata[jrow,icol,:])
                irout[jrow,icol,:] .= NaN
            end
        end
    end
    for jcol in 1:size(irdata, 2)
        for irow in lim:size(irdata, 1)
            if isnan(irdata[irow, jcol, round(Int, (1+size(irdata, 3))/2)])
                snowsurf[jcol] = irow
                break
            end
        end
    end
    return irout
end

"Apply picture-wise spatial filtering using Gauss-kernel with sigmarow and sigmacol"
function spatialgaussfilter(irdata::Array, sigmrow::Int, sigmcol::Int):Array
    @info("Picture-wise spatial Gauss-filtering with sigma_row=", sigmrow,
    " and sigma_col=", sigmcol)
    out = similar(irdata)
    for i in 1:size(irdata, 3)
        out[:,:,i] = imfilter(irdata[:,:,i], Kernel.gaussian((sigmrow, sigmcol)), NA())
    end
    return out
end

######################################################
###             FOURIER TRANSFORMATION             ###
######################################################
"After Stull p.307ff. Detrend data with linear fit"
function detrend(invec::Vector)::Vector
    model(x, p) = p[1] * x .+ p[2]
    p0 = [0.0001, 3.0]
    fit = curve_fit(model, collect(1:size(invec, 1)), invec, p0)
    param = fit.param
    return invec - (param[1]*collect(1:size(invec, 1)) .+ param[2])
end

"After Stull p.310 Apply bell tapering to smoothen edges."
function belltaper(invec::Vector)::Vector
    #check equation 8.4.3 on p.310 in Stull
    outvec = similar(invec)
    N = size(invec, 1)
    lim1 = round(Int, 0.1*N)
    lim2 = round(Int, 0.9*N)
    for i in eachindex(invec)
        if (i <= lim1) || (i >= lim2)
            outvec[i] = invec[i] * (sin(5*pi*i/N))^2
        else
            outvec[i] = invec[i]
        end
    end
    return outvec
end

"Fourier-trafo with the necessary folding according to stull"
function fourierplus(invec::Vector)
    #compare JuliaFFT-doc (https://juliamath.github.io/AbstractFFTs.jl/stable/api/#Public-Interface-1)
    #with Stull page 303 (8.4.1b) to see the division by N
    FTvec = fft(invec)/length(invec) #(8.4.1b with n->n+1 (Stull n=0 -> here n=1))

    #folding back the second half according to Stull p. 313
    #resulting in dexcrete spectral intensity (or energy) E_A(n)
    nyfreq = ceil(Int, length(FTvec)/2.0)  #Nyquist frequency
    EAvec = zeros(Float64, nyfreq)
    EAvec[1] = abs.(FTvec[1])^2

    if mod(length(FTvec),2) == 0 #even
        EAvec[2:end-1] = 2 * abs.(FTvec[2:nyfreq-1]).^2
        EAvec[end] = abs.(FTvec[nyfreq])^2
    else
        EAvec[2:end] = 2 * abs.(FTvec[2:nyfreq]).^2
    end
    return FTvec*length(invec), EAvec
end

"Converting to spectral energy density (c.f. (8.6.2b) in Stull)"
function spectralenergydensity(FTin::Vector, duration)
    freq = zeros(Float64, length(FTin))
    for i in 2:length(freq)
        freq[i] = (i-1)/duration
    end
    Sout = similar(FTin)
    Sout[1] = FTin[1]
    Sout[2:end] = FTin[2:end]*duration
    return Sout, freq
end

"For lowpass-filter, below a certain freq 'flim' return idx of next bigger entry in freq-vector"
function findcutofffreq(freqin::Vector, flim)::Int
    #find index idxflim of limit freq >= flim in freqin
    idxflim = size(freqin, 1)
    for i in 1:size(freqin, 1)
        if freqin[i] >= flim
            idxflim = i
            break
        end
    end
    #println("Desired LPF cutoff-frequency: ", flim, " Hz,")
    #println("actual LPF cutoff-frequency: ", freqin[idxflim], " Hz")
    return idxflim
end

"Apply a low-pass filter and return the Fourier-backtransformed data"
function lowpass(ftrawin::Vector, idxlim::Int)::Vector
    lpft = ftrawin
    lpft[idxlim:end-idxlim+2] .= 0 #setting desired freqs to 0

    #backtrafo
    lpift = ifft(lpft)
    lpiftreal = zeros(Float64, length(lpift))
    for i in eachindex(lpift)
        if imag(lpift[i]) > exp10(-13)
            println("ERROR: Return of inverse FT is not REAL. Please check.")
            println("imaginary part: ", imag(lpift[i]))
            println("EXIT. Returning 0!")
            return [0.0]
        end
    end
    lpiftreal = real.(lpift)
    return lpiftreal
end

"Inverse bell-tapering after Stull."
function ibelltaper(invec::Vector)::Vector
    #check equation 8.4.3 on p.310 in Stull
    outvec = similar(invec)
    N = size(invec, 1)
    lim1 = round(Int, 0.1*N)
    lim2 = round(Int, 0.9*N)
    for i in eachindex(invec)
        if (i <= lim1) || (i >= lim2)
            outvec[i] = invec[i] / (sin(5*pi*i/N))^2
        else
            outvec[i] = invec[i]
        end
    end
    return outvec
end

"Apply low-pass filter to input vector. All-in-one-function."
function lowpassallinone(invect::Vector, dur, fco)::Vector
    #apply bell taper and detrending
    vec1 = belltaper(detrend(invect))

    #apply Fourier-Trafo
    (FTraw, FTvec) = fourierplus(vec1)
    (temp,freq) = spectralenergydensity(FTvec, dur)

    #apply inverse Fourier-Trafo
    idxfco = findcutofffreq(freq, fco)

    #back-trafo
    lpift = lowpass(FTraw, idxfco)
    #apply inverse bell tapering
    lpf = ibelltaper(lpift)
    return lpf
end

"Apply 'lowpassallinone' for a matrix"
function fouriercutoffmatrix(inmatrix::Array, dur, fco)
    out = similar(inmatrix)
    IRdetrend = similar(inmatrix)
    @info("Lowpass filtering for array...")
    p = Progress(size(inmatrix, 2))
    update!(p, 0)
    jj=Threads.Atomic{Int}(0)
    l=Threads.SpinLock()
    Threads.@threads for i in 1:size(inmatrix, 2)
      #  println(i, "/", size(inmatrix, 2))
        for j in 1:size(inmatrix, 1)
            if !isnan(maximum(inmatrix[j,i,:]))
                IRdetrend[j,i,:] = detrend(inmatrix[j,i,:])
                out[j,i,:] = lowpassallinone(inmatrix[j,i,:], dur, fco)
            else
                out[j,i,:] .= NaN
                IRdetrend[j,i,:] .= NaN
            end
        end
        Threads.atomic_add!(jj,1)
        Threads.lock(l)
        update!(p,jj[])
        Threads.unlock(l)
    end
    return out, IRdetrend
end

######################################################
###              UNIFY ALL FUNCTIONS               ###
######################################################
"Performing the actual preprocessing and saving the output"
function preprocessing(sourcefile::String, gaussparam::Vector,
    fps::Float64, backlookingavgtime::Float64, lp_threshold::Float64,
    outfolder::String)
    @info("Preprocessing IR-data for WEIRD-windfield-estimation")
    #create outputfolder
    mkpath(outfolder)
    #irin = fromdirloadMATfile(joinpath(sourcefolder, fl[1]))
    irin = loadfromNetCDF4(sourcefile, "irdata")
    lenfile = size(irin, 3)
    #setting the snow to NaN and spatial Gauss filtering
    irtmp = setsnowsurfacetonan(irin)
    irtmp = spatialgaussfilter(irtmp, gaussparam[1], gaussparam[2])
    #high-pass filtering = subtract backlooking time average
    irtmp2 = subtractbacklookingtimeaverage(irtmp,round(Int, backlookingavgtime*fps))
    #low-pass filtering using Fourier transformation
    (irout, ) = fouriercutoffmatrix(irtmp2, round(Int, fps*size(irtmp2, 3)), lp_threshold)
    saveirasnetcdf(irout, "irdata", string(outfolder, sourcefile[end-6:end-3], "_preproc.nc"), 8)
end
end #module