function qssim_matrix = qssimMontage(list_methods)

if nanargin<1
    im_dir = 'C:\Users\tam128\Documents\MATLAB\Norm\Images';
    imlist = dir(fullfile(im_dir,'*.tif'));
    imlist = {imlist.name}';
    
    for m=1:5
        targetname = imlist{m};
        target = imread(fullfile(im_dir,targetname));
        for n=1:length(imlist)
            imname = imlist{n};
            if strcmp(targetname, imname)== 0
                source_im = imread(fullfile(im_dir,imname));
                [list_methods, ~] = normalize(source_im, target);
                num_methods = length(list_methods);
                
                for i = 1:num_methods
                    for j = 1:num_methods
                        image1 = list_methods{i}; im1 = image1;
                        image1(:,:,2:3) = [];
                        [h1, w1] = size(image1);
                        image2 = list_methods{j}; im2 = image2;
                        image2(:,:,2:3) = [];
                        [h2, w2] = size(image2);
                        
                        if h1 > h2
                            im1 = imresize(im1, [h2, w2]);
                            im2 = list_methods{j};
                        elseif h1 < h2
                            im2 = imresize(im2, [h1, w1]);
                            im1 = list_methods{i};
                        end
                        
                        [mqssim, ~] = qssim(im1, im2);
                        qssim_matrix(i, j) = mqssim;
                    end
                    if i>=3
                        %qssim_total((i-2), 1) = qssim_total((i-2), 1) + qssim_matrix(2, i);
                        qssim_total((i-2), 1) = qssim_total((i-2), 1) + qssim_matrix(1, i);     %Renorm
                    end
                end
                
                qssim_figure = figure('position', [200, 120, 1500, 800]);
                imagesc(qssim_matrix); colormap autumn; colorbar;
                textStrings = num2str(qssim_matrix(:),'%0.4f');
                textStrings = strtrim(cellstr(textStrings));
                [x,y] = meshgrid(1:10);
                hStrings = text(x(:),y(:),textStrings(:),'HorizontalAlignment','center');
                midValue = mean(get(gca,'CLim'));
                textColors = repmat(qssim_matrix(:) > midValue,1,3);
                labels = {'Source'; 'Target'; 'Luong'; 'Macenko'; 'Reinhard'; 'Khan'; 'Khan Hist'; 'Khan Leeds'; 'Vahadane'; 'Vahadane Fast'};
                set(gca,'XTick', [1:10]);
                set(gca,'XTickLabel', labels); set(gca,'xaxisLocation','top')
                set(gca,'YTick',[1:10]); set(gca,'YTickLabel', labels);
                title('QSSIM Matrix', 'FontSize', 14);
                
                printfile('C:\Users\tam128\Documents\MATLAB\Norm\QSSIM\Random target QSSIM\Renormalized target\QSSIM Renormalized\',...
                    qssim_figure, targetname, imname);
            end
        end
    end
end
end