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
include("./code/scorr.jl")

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

function do_space_corr()
    for mouse in settings.mouse_list
        process_mouse(mouse,threshold=0.,rotations=9)
    end
end
