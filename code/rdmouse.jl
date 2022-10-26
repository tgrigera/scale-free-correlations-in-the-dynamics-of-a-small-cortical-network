# rdmouse.jl -- reading mouse experimental data
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

include("./settings.jl")

using DelimitedFiles
import Statistics

struct MouseNotFound <: Exception
    file::String
end

Base.showerror(io::IO, e::MouseNotFound) = print(io, "mouse not found: file ", e.file)

function mouse_files(mouse,layer)
    mdir=settings.edata_dir*"/"*mouse*"/"
    return mdir*mouse*"_layer$(layer)_positions.dat",mdir*mouse*"_layer$(layer)_ca.dat",
    mdir*mouse*"_layer$(layer)_positions.bin",mdir*mouse*"_layer$(layer)_ca.bin"
end

function mouse_read_bin_matrix(fname)
    mat=zeros(Float64,2,2)
    nf=size(mat,1)
    f=open(fname*"_size","r")
    nf=read(f,typeof(nf))
    nc=read(f,typeof(nf))
    close(f)
    f=open(fname,"r")
    mat=zeros(Float64,nf,nc)
    read!(f,mat)
    close(f)
    return mat
end

""" 
    mouse_pos_sig(mouse,layer;raw=false)

Read positions and activity of requested `mouse` and `layer`.  If
available, read from binary file, if not, use ASCII file and write
binary file.  `raw` means do not shift origin of positions.
"""
function mouse_pos_sig(mouse,layer;raw=false)
    posf,sigf,posbf,sigbf=mouse_files(mouse,layer)
    if !isfile(posf)
	throw(MouseNotFound(posf))
    elseif !isfile(sigf)
	throw(MouseNotFound(sigf))
    end
    mpos = isfile(posbf) ? mouse_read_bin_matrix(posbf) : readdlm(posf)
    msig = isfile(sigbf) ? mouse_read_bin_matrix(sigbf) : readdlm(sigf)
    if !isfile(posbf)
        open(posbf*"_size","w") do f write(f,size(mpos)...) end
        open(posbf,"w") do f write(f,mpos); end
    end
    if !isfile(sigbf)
        open(sigbf*"_size","w") do f write(f,size(msig)...) end
        open(sigbf,"w") do f write(f,msig) end
    end
    if raw return mpos,msig
    else
        xmin,ymin=minimum(mpos[:,1]),minimum(mpos[:,2])
        mpos.-=[xmin ymin]
    end
    return mpos,msig
end

"""
    mouse_pos_sig_3D(mouse;raw=false)

Read all layers of `mouse` as 3-d data.  `raw` means do not shift
origin of positions.
"""
function mouse_pos_sig_3D(mouse;raw=false)
    p3=[]
    s3=[]
    for l in 1:11
        try
            pos,sig=mouse_pos_sig(mouse,l,raw=true)
            z=35*(l-1)*ones(size(pos,1))
            pos=hcat(pos,z)
            p3= l>1 ? vcat(p3,pos) : pos
            s3= l>1 ? vcat(s3,sig) : sig
        catch err
            if !isa(err,MouseNotFound)
                rethrow(err)
            end
        end
    end
    if raw return p3,s3
    else
        xmin,ymin=minimum(p3[:,1]),minimum(p3[:,2])
        p3.-=[xmin ymin 0]
    end
    return p3,s3
end

"""
     rotate_pos(mpos,θ)

Rotate all (2-d) postions counterclockwise by the given angle,
with respect to the center of the region.
"""
function rotate_pos(mpos,θ)
    rmat=[ cos(θ)  -sin(θ) ; sin(θ) cos(θ) ]
    xmin,xmax=extrema(mpos[:,1])
    ymin,ymax=extrema(mpos[:,2])
    rc = [(xmin+xmax)/2 , (ymin+ymax)/2]
    rpos=zero(mpos)
    for i=1:size(mpos,1)  rpos[i,:]=rmat*(mpos[i,:]-rc) + rc end
    return rpos
end

"""
    csig_to_dsig(csig;threshold=0.)

Transform continuous signal `csig` to a point process where
inactive=-1, active=0.  Active sites are those where `csig>threshold`.
"""
function csig_to_dsig(csig;threshold=0.)
    dsig=zeros(Int8,size(csig))         # set all inactive
    dsig[csig.>threshold].=1            # turn on if continous sig above thershold
    return dsig
end

"""
    filter_inactive(mpos,msig)

Return new `(mpos, msig)` tuple where sites with no activity over the
whole time series are filtered out.
"""
function filter_inactive(mpos,msig)

    amask=reshape(sum(msig,dims=2).>0,size(msig,1))
    return mpos[amask,:],msig[amask,:]
end

"This will filter on the discrete signal"
function filter_inactive(mpos,msig,threshold,firing_threshold)
    dsig=csig_to_dsig(msig,threshold=threshold)
    amask=reshape(sum(dsig,dims=2).>=firing_threshold,size(dsig,1))
    return mpos[amask,:],msig[amask,:]
end
