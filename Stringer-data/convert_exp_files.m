% This routine converts .mat experimental files to .dat
% Files are deposited at:
% https://doi.org/10.25378/janelia.6845348.v4][DOI 10.25378/janelia.6845348.v4

clear all
close all
clc


Files = dir('spont_M1*.mat');
mouse_labels = ["mouse1" "mouse2" "mouse3" "mouse4a" "mouse4b" "mouse4c" "mouse5" "mouse6" "mouse7"] ;

% file_number: 1, mouse1, spont_M150824_MP019_2016-04-05.mat
% file_number: 2, mouse2, spont_M160825_MP027_2016-12-12.mat
% file_number: 3, mouse3, spont_M160907_MP028_2016-09-26.mat
% file_number: 4, mouse4a, spont_M161025_MP030_2016-11-20.mat
% file_number: 5, mouse4b, spont_M161025_MP030_2017-06-16.mat
% file_number: 6, mouse4c, spont_M161025_MP030_2017-06-23.mat
% file_number: 7, mouse5, spont_M170714_MP032_2017-08-04.mat
% file_number: 8, mouse6, spont_M170717_MP033_2017-08-18.mat
% file_number: 9, mouse7, spont_M170717_MP034_2017-08-25.mat

% Select file
file_number = 9;

clear mouse calc_sig x y z layer

mouse = load(eval(sprintf('Files(%d).name', file_number)));
fprintf(1, '\n Loaded %s, %s\n', mouse_labels(file_number), Files(file_number).name)

[total_n,total_time] = size(mouse.Fsp);
fprintf(1, '\nneurons: %d\ntime: %d\n', total_n, total_time);

calc_sig = double(mouse.Fsp(1:total_n,1:total_time));
x = mouse.med(1:total_n,1);
y = mouse.med(1:total_n,2);
z = mouse.med(1:total_n,3);

layer = unique(z);

fprintf(1, '\n\n%d layers\n\n', length(layer));


for k=1:length(layer)
    clear z_layer zi nn_z xx yy positions
    
    z_layer = layer(k);
    zi = find(z==z_layer);
    nn_z = length(zi);
    
    calc_sig_z = calc_sig(zi,:);
    xx = x(zi);
    yy = y(zi);
    
    positions = [xx yy];
        
    dlmwrite(sprintf('%s_layer%d_positions.dat', mouse_labels(file_number), k), positions, 'delimiter', '\t')
    dlmwrite(sprintf('%s_layer%d_ca.dat', mouse_labels(file_number), k), calc_sig_z, 'delimiter', '\t','precision', 4)
           
end


