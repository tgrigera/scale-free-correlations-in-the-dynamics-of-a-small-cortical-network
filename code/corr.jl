# corr.jl -- Functions for correlation functions
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

"""
    boxing(pos,W)

Return an array of masks that select a subset of all (2-d) positions
lying in boxes of side `W`.
"""
function boxing(pos,W)
    boxes=[]
    L=maximum(pos)
    if W<L
        r=0:W:L-W
        for x in r,y in r
            mask= x.<=pos[:,1].<=x+W .&& y.<=pos[:,2].<=y+W
            push!(boxes,mask)
        end
    else
        mask=trues(size(pos,1))
        push!(boxes,mask)
    end
    return boxes
end

"""
    boxing_3d(pos,W)

Return an array of masks that select a subset of all (3-d) positions
lying in boxes of side `W`.
"""
function boxing_3d(pos,W)
    boxes=[]
    Lx=maximum(pos[:,1])
    Ly=maximum(pos[:,2])
    Lz=maximum(pos[:,3])
    L=min(Lx,Ly,Lz)
    r=0
    if W<L
        r=0:W:L-W
    end
    for x in r,y in r,z in r
        mask= x.<=pos[:,1].<=x+W .&& y.<=pos[:,2].<=y+W .&& z.<=pos[:,3].<=z+W
        push!(boxes,mask)
    end
    return boxes
end

function GR(binning::DistanceBinning,signal)
    GRden=Float64.(length.(binning))
    Ntot=size(signal,1)
    GRden./=Ntot*Ntot

    GR=zeros(Float64,length(GRden))
    ntimes=size(signal,2)
    for time=1:ntimes
        Nact=count( x->x>0,signal[:,time])
        if Nact==0 continue end
        grt=zeros(Int32,length(GR))
        for ib in eachindex(binning)
            for (i,j) in binning[ib]
                grt[ib]+=signal[i,time]*signal[j,time]
            end
        end
        GR .+= grt / (Nact*Nact)
    end
    GR/=ntimes
    GR./=GRden
    return vcat(0,collect(range(binning))),vcat(0,GR)
end    

function three_corrs(pos,signal,Δr;rmax=nothing,normalize=false,threshold=0.,return_pair_count=false)
    bins=distance_binning(pos,Δr,rmax=rmax)
    r,Cs=space_correlation(bins,signal,connected=true,normalized=normalize)
    dsignal=csig_to_dsig(signal,threshold=threshold)
    r,Cds=space_correlation(bins,dsignal,connected=true,normalized=normalize)
    r,G=GR(bins,dsignal)
    npairs = return_pair_count ? vcat(0,length.(bins.pairs)) : nothing
    return return_pair_count ? (r,Cs,Cds,G,npairs) : (r,Cs,Cds,G)
end

function three_corrs_box(mpos,signal,W,Δr;normalize=false,threshold=0.,
                         rotations=1,return_pair_count=false)
    is3d= size(mpos,2)==3
    rmax=is3d ?  sqrt(3)*W : sqrt(2)*W

    Crnbins = Int(ceil(rmax/Δr))+1
    Cr=zeros(Float64,Crnbins)
    nsCr=zero(Cr)
    Crd=zero(Cr)
    nsCrd=zero(Cr)
    Gr=zero(Cr)
    nsGr=zero(Cr)
    r=zero(Cr)

    if is3d rotations=1 end
    npairs=0
    for irot=0:rotations-1
        θ=irot*2π/rotations
        pos= θ > 0 ? rotate_pos(mpos,θ)  : mpos
        bxs= is3d ? boxing_3d(pos,W) : boxing(pos,W)
        for box in bxs
            tseries=sum(signal[box,:],dims=1)
            if any(tseries.>0)
                ret=three_corrs(pos[box,:],signal[box,:],Δr;rmax=rmax,normalize=normalize,threshold=threshold,
                                     return_pair_count)
                if return_pair_count r,c,cd,g,npairs=ret else r,c,cd,g=ret end
                mask=.!(isnan.(c))
                Cr[mask].+=c[mask]
                nsCr[mask].+=1
                mask=.!(isnan.(cd))
                Crd[mask].+=cd[mask]
                nsCrd[mask].+=1
                mask=.!(isnan.(g))
                Gr[mask].+=g[mask]
                nsGr[mask].+=1
            end
        end
    end

    Cr./=nsCr
    Crd./=nsCrd
    Gr./=nsGr
    return return_pair_count ? (r,Cr,Crd,Gr,npairs) : (r,Cr,Crd,Gr)
end

function tcorr_self(signal,Δt=1;normalize=false)
    corr=zeros(Float64,size(signal,2)÷2)
    nparts=0
    for i=1:size(signal,1)
        if any(signal[i,:].>0)
            c=time_correlation(signal[i,:],connected=true,normalized=normalize)
            corr .+= c
            nparts+=1
        end
    end
    return Δt*collect(0:size(corr,1)-1),corr/nparts
end

function tcorr_box(mpos,signal,W,Δt=1;rotations=1)
    is3d= size(mpos,2)==3
    rmax=is3d ?  sqrt(3)*W : sqrt(2)*W

    someneg=any(signal.<0) 
    corr=zeros(Float64,size(signal,2)÷2)
    c=zero(corr)
    tseries=zeros(Float64,size(signal,2))
    nbxs=0
    for irot=0:rotations-1
        θ=irot*2π/rotations
        pos= θ > 0 ? rotate_pos(mpos,θ)  : mpos
        bxs= is3d ? boxing_3d(pos,W) : boxing(pos,W)
        for box in bxs
            tseries=sum(signal[box,:],dims=1)[1,:]
            if someneg || any(tseries.>0)
                c=time_correlation(tseries,connected=true,normalized=false)
                corr .+= c
                nbxs+=1
            end
        end
    end
    return Δt*collect(0:size(corr,1)-1),corr/nbxs
end

function tau_threshold(t,C,threshold=0.1)
    it=findfirst(x->x<threshold,C/C[1])
    if !isnan(C[it-1]) return t[it-1] - (C[it-1]/C[1]-threshold)*(t[it]-t[it-1])/((C[it]-C[it-1])/C[1])
    else return t[it]
    end
end    
