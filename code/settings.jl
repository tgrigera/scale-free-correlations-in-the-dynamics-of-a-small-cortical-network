# settings.jl -- load defaults and directory locations
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

push!(LOAD_PATH,homedir()*"/software/BioStatPhys.jl")

using BioStatPhys

Base.@kwdef struct Settings
    Δr::Float64 = 5
    Wlist = 100:100:1000
    layerlist = 1:11

    mouse_list = [ "mouse1", "mouse2", "mouse3", "mouse4a", "mouse4b", "mouse4c", "mouse5", "mouse6", "mouse7"]
    sfreq_list = [ 3, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5]

    edata_dir = "./Stringer-data"
    space_corr_dir = "./space_corr"
end

settings = Settings()
