%% sssss
figure; imshow(target_im); zoom(2);
target_rect_pink = getrect;
target_pink = imcrop(target_im, target_rect_pink);
target_pink_rgb = raw2rgb(target_pink);
imshow(target_im);zoom(2);
target_rect_purple = getrect;
target_purple = imcrop(target_im, target_rect_purple);
target_purple_rgb = raw2rgb(target_purple);
[ target_pink_oppCol, target_pink_brightness, target_pink_theta, target_pink_sat] = rgb2oppCol( target_pink_rgb, target_rotation_matrix);
[ target_purple_oppCol, target_purple_brightness, target_purple_theta, target_purple_sat] = rgb2oppCol( target_purple_rgb, target_rotation_matrix); 
%
% figure;
% subplot(3,2,1); imshow(target_pink);
% subplot(3,2,2); imshow(target_purple);
% %
% subplot(3,2,3);scatter(target_pink_oppCol(1,:),target_pink_oppCol(2,:),20,target_pink_rgb'./255,'filled');
% axis([-1 1 -1 1]); axis square
% subplot(3,2,4);scatter(target_purple_oppCol(1,:),target_purple_oppCol(2,:),20,target_purple_rgb'./255,'filled');
% axis([-1 1 -1 1]); axis square
% %
% subplot(3,2,5);circ_plot(target_pink_theta(:),'hist',[],40,true,true,'linewidth',2,'color','r');
% subplot(3,2,6);circ_plot(target_purple_theta(:),'hist',[],40,true,true,'linewidth',2,'color','r');

figure;
subplot(1,2,1); imshow(target_pink);
subplot(1,2,2); imshow(target_purple);
%
figure;
subplot(1,2,1);scatter(target_pink_oppCol(1,:),target_pink_oppCol(2,:),20,target_pink_rgb'./255,'filled');
axis([-1 1 -1 1]); axis square
subplot(1,2,2);scatter(target_purple_oppCol(1,:),target_purple_oppCol(2,:),20,target_purple_rgb'./255,'filled');
axis([-1 1 -1 1]); axis square
%
figure;
subplot(1,2,1);circ_plot(target_pink_theta(:),'hist',[],40,true,true,'linewidth',2,'color','r');
subplot(1,2,2);circ_plot(target_purple_theta(:),'hist',[],40,true,true,'linewidth',2,'color','r');
%
bincol = [0.8,0.8,0.8];
legendstr = {'pink','purple'};
A = {target_pink_brightness,target_purple_brightness};
figure; nhist(A,'legend',legendstr,'color',bincol,...
    'separate','samebins','noerror','smooth','binfactor',1,'smooth','pdf');
B = {target_pink_sat,target_purple_sat};
figure; nhist(B,'legend',legendstr,'color',bincol,...
    'separate','samebins','noerror','smooth','binfactor',1,'smooth','pdf');

%%
figure; imshow(source_im);zoom(2);
source_rect_pink = getrect;
source_pink = imcrop(source_im, source_rect_pink);
source_pink_rgb = raw2rgb(source_pink);
imshow(source_im);zoom(2);
source_rect_purple = getrect;
source_purple = imcrop(source_im, source_rect_purple);
source_purple_rgb = raw2rgb(source_purple);
[ source_pink_oppCol, source_pink_brightness, source_pink_theta, source_pink_sat] = rgb2oppCol( source_pink_rgb, source_rotation_matrix);
[ source_purple_oppCol, source_purple_brightness, source_purple_theta, source_purple_sat] = rgb2oppCol( source_purple_rgb, source_rotation_matrix); 
%
% figure;
% subplot(3,2,1);imshow(source_pink);
% subplot(3,2,2);imshow(source_purple);
% %
% subplot(3,2,3);scatter(source_pink_oppCol(1,:),source_pink_oppCol(2,:),20,source_pink_rgb'./255,'filled');
% axis([-1 1 -1 1]); axis square
% subplot(3,2,4);scatter(source_purple_oppCol(1,:),source_purple_oppCol(2,:),20,source_purple_rgb'./255,'filled');
% axis([-1 1 -1 1]); axis square
% %
% subplot(3,2,5);circ_plot(source_pink_theta(:),'hist',[],40,true,true,'linewidth',2,'color','r');
% subplot(3,2,6);circ_plot(source_purple_theta(:),'hist',[],40,true,true,'linewidth',2,'color','r');
figure;
subplot(1,2,1);imshow(source_pink);
subplot(1,2,2);imshow(source_purple);
%
figure;
subplot(1,2,1);scatter(source_pink_oppCol(1,:),source_pink_oppCol(2,:),20,source_pink_rgb'./255,'filled');
axis([-1 1 -1 1]); axis square
subplot(1,2,2);scatter(source_purple_oppCol(1,:),source_purple_oppCol(2,:),20,source_purple_rgb'./255,'filled');
axis([-1 1 -1 1]); axis square
%
figure;
subplot(1,2,1);circ_plot(source_pink_theta(:),'hist',[],40,true,true,'linewidth',2,'color','r');
subplot(1,2,2);circ_plot(source_purple_theta(:),'hist',[],40,true,true,'linewidth',2,'color','r');

A = {source_pink_brightness,source_purple_brightness};
figure; nhist(A,'legend',legendstr,'color',bincol,...
    'separate','samebins','noerror','smooth','binfactor',1,'smooth','pdf');
B = {source_pink_sat,source_purple_sat};
figure; nhist(B,'legend',legendstr,'color',bincol,...
    'separate','samebins','noerror','smooth','binfactor',1,'smooth','pdf');

%% Then you should do the same for the image after normalization
source_eq_pink = imcrop(source_eq_image, source_rect_pink);
source_eq_pink_rgb = raw2rgb(source_eq_pink);
source_eq_purple = imcrop(source_eq_image, source_rect_purple);
source_eq_purple_rgb = raw2rgb(source_eq_purple);

[ source_eq_oppCol, source_eq_brightness, source_eq_theta, source_eq_sat] = rgb2oppCol( double(source_rgb_eq_uint8), target_rotation_matrix);

[ source_eq_pink_oppCol, source_eq_pink_brightness, source_eq_pink_theta, source_eq_pink_sat] = rgb2oppCol( source_eq_pink_rgb, target_rotation_matrix);
[ source_eq_purple_oppCol, source_eq_purple_brightness, source_eq_purple_theta, source_eq_purple_sat] = rgb2oppCol( source_eq_purple_rgb, target_rotation_matrix); 
%
%% 
% figure;
% subplot(3,2,1);imshow(source_eq_pink);
% subplot(3,2,2);imshow(source_eq_purple);
% %
% subplot(3,2,3);scatter(source_eq_pink_oppCol(1,:),source_eq_pink_oppCol(2,:),20,source_eq_pink_rgb'./255,'filled');
% axis([-1 1 -1 1]); axis square
% subplot(3,2,4);scatter(source_eq_purple_oppCol(1,:),source_eq_purple_oppCol(2,:),20,source_eq_purple_rgb'./255,'filled');
% axis([-1 1 -1 1]); axis square
% %
% subplot(3,2,5);circ_plot(source_eq_pink_theta(:),'hist',[],40,true,true,'linewidth',2,'color','r');
% subplot(3,2,6); circ_plot(source_eq_purple_theta(:),'hist',[],40,true,true,'linewidth',2,'color','r');
figure;
subplot(1,2,1);imshow(source_eq_pink);
subplot(1,2,2);imshow(source_eq_purple);
%
figure;
subplot(1,2,1);scatter(source_eq_pink_oppCol(1,:),source_eq_pink_oppCol(2,:),20,source_eq_pink_rgb'./255,'filled');
axis([-1 1 -1 1]); axis square
subplot(1,2,2);scatter(source_eq_purple_oppCol(1,:),source_eq_purple_oppCol(2,:),20,source_eq_purple_rgb'./255,'filled');
axis([-1 1 -1 1]); axis square
%
figure;
subplot(1,2,1);circ_plot(source_eq_pink_theta(:),'hist',[],40,true,true,'linewidth',2,'color','r');
subplot(1,2,2); circ_plot(source_eq_purple_theta(:),'hist',[],40,true,true,'linewidth',2,'color','r');
%
A = {source_eq_pink_brightness,source_eq_purple_brightness};
figure; nhist(A,'legend',legendstr,'color',bincol,...
    'separate','samebins','noerror','smooth','binfactor',1,'smooth','pdf');
B = {source_eq_pink_sat,source_eq_purple_sat};
figure; nhist(B,'legend',legendstr,'color',bincol,...
    'separate','samebins','noerror','smooth','binfactor',1,'smooth','pdf');

%% why not bimodal?
bincol = [0.8,0.8,0.8];
legendstr = {'purple','pink','inbetween','just pink + purple','all together'};
purple = normrnd(1,1,1500,1);
pink = normrnd(5,1,4000,1);
inbetween = normrnd(3,1,1000,1);
justpp = [purple; pink];
together = [purple; pink; inbetween];
A = {purple, pink, inbetween, justpp, together};
figure; nhist(A,'legend',legendstr,'color',bincol,...
    'separate','samebins','noerror','smooth','binfactor',1,'smooth','pdf');
