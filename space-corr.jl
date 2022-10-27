# space-corr.jl -- Compute space correlation functions of mouse data
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

function process_mouse(mouse;threshold=0.,filterth=0,layers=settings.layerlist,rotations=1)
    print("* $mouse w/ threshold $threshold, filtering = $filterth\n")
    corrdir=settings.space_corr_dir*"/athresh=$(threshold)_fthresh=$(filterth)_rotations=$(rotations)"
    mkpath(corrdir)
    for layer in layers
        try
            print("** Reading layer $(layer)\n")
            mpos,msig=mouse_pos_sig(mouse,layer)
            Nneuro_before=size(mpos,1)
            print("   Read $Nneuro_before neurons, $(size(msig,2)) times\n")
            if filterth!=0
                print("** Filtering with $filterth\n")
                if filterth<=0
                    mpos,msig=filter_inactive(mpos,msig)
                else
                    mpos,msig=filter_inactive(mpos,msig,threshold,filterth)
                end
            end
            Nneuro=size(mpos,1)
            print("   After filtering: $Nneuro neurons, $(size(msig,2)) times\n")
            for W in settings.Wlist
                print("*** Computing corrs in box, W=$W\n")
                r,c,cd,g=three_corrs_box(mpos,msig,W,settings.Δr,normalize=true,threshold=threshold,rotations=rotations)
                if count(!isnan,c)>0
                    cfile=corrdir*"/Cr_$(mouse)_layer_$(layer)_W$W.dat"
                    writedlm(cfile,[r c])
                else
                    print("** WARN: Layer $layer W=$W: C(r) failed\n")
                end
                if count(!isnan,cd)>0
                    cfile=corrdir*"/Crd_$(mouse)_layer_$(layer)_W$W.dat"
                    writedlm(cfile,[r cd])
                else
                    print("** WARN: Layer $layer W=$W: Cd(r) failed\n")
                end
                if count(!isnan,g)>0
                    cfile=corrdir*"/Gr_$(mouse)_layer_$(layer)_W$W.dat"
                    writedlm(cfile,[r g])
                else
                    print("** WARN: Layer $layer W=$W: G(r) failed\n")
                end
            end
            print("*** Computing global corrs(r)\n")
            r,c,cd,g=three_corrs(mpos,msig,settings.Δr,normalize=true,threshold=threshold)
            cfile=corrdir*"/Cr_$(mouse)_layer_$(layer)_global.dat"
            writedlm(cfile,[r c])
            cfile=corrdir*"/Crd_$(mouse)_layer_$(layer)_global.dat"
            writedlm(cfile,[r cd])
            cfile=corrdir*"/Gr_$(mouse)_layer_$(layer)_global.dat"
            writedlm(cfile,[r g])
        catch err
            if isa(err,MouseNotFound)
                print("** WARN: Layer $layer not processed, missing file\n")
            else
                rethrow(err)
            end
        end
    end
end

function findr0(mouse,layerlist,Wlist,pref,corrdir)
    tt=[]
    for layer in layerlist
        for W in Wlist
            cfile=corrdir*"/$(pref)_$(mouse)_layer_$(layer)_W$W.dat"
            if !isfile(cfile) continue end
            C=readdlm(cfile)
            # ad-hoc corrections for mouse2 and 3
            if (mouse=="mouse2" && layer==1) || (mouse=="mouse3" && (layer==9 || layer==10))
                    C=C[1:end .!=2,:]
            end
            r0=correlation_length_r0(C[:,1],C[:,2])
            push!(tt,[mouse layer W r0])
        end
    end
    if length(tt)==0 return [missing missing missing missing] end
    tt=vcat(tt...)    # Reallocates the whole thing, not good for large arrays
    for W in Wlist
        mask=tt[:,3].==W
        if count(mask)<1 continue end
        rm=Statistics.mean(tt[mask,4])
        tt=vcat(tt,[mouse "avg" W rm])
    end
    return tt
end

function findr0Gr(mouse,layerlist,Wlist,corrdir)
    tt=[]
    for layer in layerlist
        for W in Wlist
            cfile=corrdir*"/Gr_$(mouse)_layer_$(layer)_W$W.dat"
            if !isfile(cfile) continue end
            C=readdlm(cfile)
            C[1,2]=1
            r0=correlation_length_r0(C[:,1],C[:,2].-1)
            push!(tt,[mouse layer W r0])
        end
    end
    if length(tt)==0 return [missing missing missing missing] end
    tt=vcat(tt...)
    for W in Wlist
        mask=tt[:,3].==W
        if count(mask)<1 continue end
        rm=Statistics.mean(tt[mask,4])
        tt=vcat(tt,[mouse "avg" W rm])
    end
    return tt
end

function do_space_corr()
    threshold=0.
    filterth=0
    rotations=9
    # Compute correlations

    for mouse in settings.mouse_list
        process_mouse(mouse,threshold=threshold,filterth=filterth,rotations=rotations)
    end

    # Compute r0

    r0Crtab=["#mouse" "layer" "W" "r0" ];
    r0Crdtab=["#mouse" "layer" "W" "r0" ];
    r0Grtab=["#mouse" "layer" "W" "r0" ];
    corrdir=settings.space_corr_dir*"/athresh=$(threshold)_fthresh=$(filterth)_rotations=$(rotations)"
    for mouse in settings.mouse_list
        print("Computing r0, $mouse, $corrdir\n")
        r0Crtab=vcat(r0Crtab,findr0(mouse,settings.layerlist,settings.Wlist,"Cr",corrdir))
        r0Crdtab=vcat(r0Crdtab,findr0(mouse,settings.layerlist,settings.Wlist,"Crd",corrdir))
        r0Grtab=vcat(r0Grtab,findr0Gr(mouse,settings.layerlist,settings.Wlist,corrdir))
    end
    writedlm(corrdir*"/r0W_Cr.dat",r0Crtab)
    writedlm(corrdir*"/r0W_Crd.dat",r0Crdtab)
    writedlm(corrdir*"/r0W_Gr.dat",r0Grtab)
end
