function printf = printfile(direc, figure, targetname, imname)

    targetfile = strrep(targetname, '.tif', '-');
    %filename = strcat(direc, targetfile, imname);
    filename = fullfile(direc, [targetfile, imname]);
    try
        print(figure, '-dtiff', filename);
    catch
        imwrite(figure, filename);
    end
    
end