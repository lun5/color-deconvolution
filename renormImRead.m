function renorm = renormImRead()
for k=1:6
    list_methods = {'Luong', 'Macenko', 'Reinhard', 'Khan', 'Vahadane', 'Vahadane Fast'};
    im_dir = strcat('C:\Users\tam128\Documents\MATLAB\Norm\15 Targets Tiles 512\', list_methods{k});
    source_dir = 'C:\Users\tam128\Documents\MATLAB\Norm\Tiles_512\';
    renorm_dir = 'C:\Users\tam128\Documents\MATLAB\Norm\15 Targets Renorm Tiles 512\';
    imlist = dir(fullfile(im_dir,'*.tif')); imlist = {imlist.name}';
    
    for m=1:length(imlist)
        imname = imlist{m};
        norm_source = imread(fullfile(im_dir,imname));
        filename = strsplit(imname, '-');
        sourcename =  filename{2};
        source = imread(fullfile(source_dir,sourcename));
        if k==1
            im2 = NormLuong(norm_source, source);
            printfile(strcat(renorm_dir, 'Luong\'), im2, '', imname);
        elseif k==2
            im2 = NormMacenko(norm_source, source);
            printfile(strcat(renorm_dir, 'Macenko\'), im2, '', imname);
        elseif k==3
            im2 = im2uint8(NormReinhard(norm_source, source));
            printfile(strcat(renorm_dir, 'Reinhard\'), im2, '', imname);
        elseif k==4
            im2 = NormSCDLeeds(norm_source, source);
            printfile(strcat(renorm_dir, 'Khan\'), im2, '', imname);
        elseif k==5
            im2 = SNMFnorm(norm_source, source);
            printfile(strcat(renorm_dir, 'Vahadane\'), im2, '', imname);
        elseif k==6
            im2 = Demo_colornorm(norm_source, source);
            printfile(strcat(renorm_dir, 'Vahadane Fast\'), im2, '', imname);
        end
    end
end