# time-corr.jl -- Compute time correlation functions of mouse data
#
# Copyright (C) 2022 by the AUTHORS
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see
# <https://www.gnu.org/licenses/>.

include("./code/rdmouse.jl")
include("./code/corr.jl")

function process_mouse(mouse,Δt;layers=settings.layerlist,rotations=1)
    print("* $mouse\n")
    tcorrdir=settings.time_corr_dir*"/rotations=$(rotations)/"
    mkpath(tcorrdir)
    for layer in layers
        try
            print("** Reading layer $(layer)\n")
            mpos,msig=mouse_pos_sig(mouse,layer)
            Nneuro=size(mpos,1)
            print("   Read $Nneuro neurons, $(size(msig,2)) times\n")
            print("*** Computing self timecorr\n")
            t,c=tcorr_self(msig,Δt)
            if count(!isnan,c)>0
                cfile=tcorrdir*"mouse_$(mouse)_layer_$(layer)_self.dat"
                writedlm(cfile,[t c/c[1]])
            end
            for W in settings.Wlist
                print("*** Computing block timecorr, W=$W, rotations=$rotations\n")
                t,c=tcorr_box(mpos,msig,W,Δt,rotations=rotations)
                if count(!isnan,c)>0
                    cfile=tcorrdir*"mouse_$(mouse)_layer_$(layer)_W$W.dat"
                    writedlm(cfile,[t c/c[1]])
                end
            end
            print("*** Computing global timecorr\n")
            tsig=sum(msig,dims=1)[1,:]
            c=time_correlation(tsig,connected=true,normalized=true)
            t=collect(0:size(c,1)-1)
            if count(!isnan,c)>0
                cfile=tcorrdir*"mouse_$(mouse)_layer_$(layer)_global.dat"
                writedlm(cfile,[Δt*t c])
            end
        catch err
            if isa(err,MouseNotFound)
                print("** WARN: Layer $layer not processed, missing file\n")
            else
                rethrow(err)
            end
        end
    end
end

function findtau(mouse,layerlist,Wlist,tcorrdir)
    tt=[]
    for layer in layerlist
        cfile=tcorrdir*"mouse_$(mouse)_layer_$(layer)_self.dat"
        if !isfile(cfile) continue end
        C=readdlm(cfile)
        tauW=correlation_time_spectral(C[:,2],C[2,1]-C[1,1])
        tauth=tau_threshold(C[:,1],C[:,2]/C[1,2])
        push!(tt,[mouse layer 0 tauW tauth])
        for W in Wlist
            cfile=tcorrdir*"mouse_$(mouse)_layer_$(layer)_W$W.dat"
            if !isfile(cfile) continue end
            C=readdlm(cfile)
            tauW=correlation_time_spectral(C[:,2],C[2,1]-C[1,1])
            tauth=tau_threshold(C[:,1],C[:,2]/C[1,2])
            push!(tt,[mouse layer W tauW tauth])
        end
    end
    return vcat(tt...)  # Reallocates the whole thing, not good for large arrays
end

function do_time_corr()
    rotations=9

    # Compute correlation functions
    for i in eachindex(settings.mouse_list)
        mouse=settings.mouse_list[i]
        Δt=1. /settings.sfreq_list[i]
        process_mouse(mouse,Δt,rotations=rotations)
    end

    # Compute correlation time

    tautab=[  "#mouse" "layer" "W" "tau_HH"  "tau_threshold"]
    tcorrdir=settings.time_corr_dir*"/rotations=$(rotations)/"
    for mouse in settings.mouse_list
        print("Doing tau, $mouse\n")
        tautab=vcat(tautab,findtau(mouse,settings.layerlist,settings.Wlist,tcorrdir))
    end
    writedlm(tcorrdir*"tauW.dat",tautab)
end
