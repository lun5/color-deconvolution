function vahadane_stains = vahadaneStain(str)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate stain vectors for hematoxylin and eosin stain                 %
%                                                                         %
% Based on code written by Abhishek Vahadane and Tingying Peng and Shadi  %
% Albarqouni and Maximilian Baust and Katja Steiger and Anna Melissa      %
% Schlitter and Amit Sethi and Irene Esposito and Nassir Navab,           %
% Booktitle:{IEEE International Symposium on Biomedical Imaging},         %
% Date-Modified: {2015-01-31 17:49:35 +0000},                             %
% Title:{Structure-Preserved Color Normalization for Histological Images},%
% Year = {2015}}                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf;
im_dir = 'C:\Users\tam128\Documents\MATLAB\Norm\tissue_images_tp10-867-1\';
imlist = dir(fullfile(im_dir,'*.tif'));
imlist = {imlist.name}';
vahadane_stains = cell(length(imlist), 1);

%% Calculate stain vectors
nstains=2;    % number of stains
lambda=0.1;   % default value sparsity regularization parameter
total_hem = zeros(1, 3);  total_eos = zeros(1,3);

for i=1:length(imlist)
    imname = imlist{i};
    source_im = imread(fullfile(im_dir,imname));
    
    if str=='fast'
        % Fast stain separation (V=WH)
        display('Fast stain separation is running....')
        tic;
        [Wi, Hi,Hiv,stains,Hso_Rmax]=Faststainsep(source_im,nstains,lambda);
        time=toc
    elseif str=='slow'
        % Slow stain separation
        display('Slow/direct stain separation is running....')
        tic;
        [Wi, Hi,Hiv,stains,Hso_Rmax]=stainsep(source_im,nstains,lambda);
        time=toc
    end
    
%% Visuals (for 2 stains)
%     figure(1);
%     set(figure(1), 'Position', [200, 200, 700, 450]);
%     subplot_tight(1,3,1,0.01);imshow(source_im);xlabel('Source')
%     subplot_tight(1,3,2,0.01);imshow(stains{1});xlabel('stain1')
%     subplot_tight(1,3,3,0.01);imshow(stains{2});xlabel('stain2')
    
    stain_graph = figure(2);
    view(0,0);      % 2-D plot, comment for 3-D plot
    set(stain_graph, 'Position', [1000, 200, 650, 650]);
    xlabel('Red'); ylabel('Green'); zlabel('Blue');
    Wi = Wi.';
    color = exp(-Wi*Hso_Rmax);
    vahadane_stains{i} = color;
    total_hem = total_hem + color(1, :);
    total_eos = total_eos + color(2, :);
    scatter3(color(1, 1), color(1, 2), color(1, 3), 50, color(1, :), 'filled');
    hold on;
    scatter3(color(2, 1), color(2, 2), color(2, 3), 50, color(2, :), 'filled');
end

%% Calculate and plot average stain vectors 

%num_values = sum(~cellfun(@isempty,vahadane_stains),2);
avg_hem = total_hem ./length(imlist);
avg_eos = total_eos ./length(imlist);

figure(2);
%h = quiver3(0, 0, 0, avg_hem(1, 1), avg_hem(1, 2), avg_hem(1, 3));
%set(h, 'Color', avg_hem(1, :)); set(h, 'LineWidth', 1.5);
hs = scatter3(avg_hem(1, 1),avg_hem(1, 2), avg_hem(1, 3), 220, avg_hem(1, :), 'p', 'filled');
set(hs, 'MarkerEdgeColor', 'w', 'LineWidth', 1); hold on;

%e = quiver3(0, 0, 0, avg_eos(1, 1), avg_eos(1, 2), avg_eos(1, 3));
%set(e, 'Color', avg_eos(1, :)); set(e, 'LineWidth', 1.5);
es = scatter3(avg_eos(1, 1),avg_eos(1, 2), avg_eos(1, 3), 220, avg_eos(1, :), 'p', 'filled');
set(es, 'MarkerEdgeColor', 'w', 'LineWidth', 1);

filename = strcat('C:\Users\tam128\Documents\MATLAB\Norm\tissue_images_tp10-867-1\Stain vectors-tissue_images_tp10-867-1.tif');
print(stain_graph, '-dtiff', filename);
save('stain_vectors_tissue_images_tp10-867-1.mat', 'vahadane_stains', 'avg_hem', 'avg_eos');

end