%% Test routine for the Sharp class



%% Reading one image
% folder where the images are located :
folder = '/Users/awojdyla/Desktop/SHARP_2015-06-23_LBNL';
mask_name = 'IntelCalMask';
folder = '/Users/awojdyla/Documents/MATLAB/SHARP/data/SHARP_2015-01-23_LBNL';
mask_name = 'IntelCalMask';
series_nb = 3;
image_nb = 1;
[ img,  meta] = MIP.read(folder,mask_name,series_nb,image_nb);
[ imgs,  metas] = MIP.read(folder,mask_name,series_nb);
%%
% or
filename = '/Users/awojdyla/Desktop/SHARP_2015-06-23_LBNL/IntelCalMask-150623-0001-0001.png';
[img, meta] = MIP.read(filename);

% or
filename = '/Users/awojdyla/Desktop/SHARP_2015-06-23_LBNL/IntelCalMask-150623-0001-0001.png';
img  = double(imread(filename));
meta = MIP.metadata(img);

imagesc(img)
title(sprintf('image : %d/%d',meta.image_number,meta.series_nsteps))
axis image off


%% Reading a series
%folder = '/Users/awojdyla/Desktop/SHARP_2015-06-23_LBNL';
[imgs, metas] = MIP.read(folder,'IntelCalMask',1);
imagesc(imgs{7})
axis image off

%%
idx = 4;
imagesc(imgs{idx})
title(sprintf('defocus : dz = %+1.2f um',metas{idx}.series_dz*1e3))
axis image off

%% Reading metadata
meta = MIP.metadata(img);
meta.ma_arg2;

%% Removing background
[ img_wo_bg , bg ] = MIP.remove_bg(img);
imagesc(img_wo_bg)
axis image off

%% Extracting the ROI
roi_size_px = 256;
[ img_roi, x_roi, y_roi ] = MIP.ROI( img, roi_size_px);
imagesc(img_roi)
axis image off

%% Adjusting the ROI
roi_x_decenter = -9;
roi_y_decenter = -15;
roi_size_px = 256;
[ img_roi, x_roi, y_roi ] = MIP.ROI( img, roi_size_px,...
                                       roi_x_decenter,roi_y_decenter);
imagesc(img_roi)
axis image off


%% Rotate an image
angle_deg = 1.2;
img_rotated = MIP.rotate( MIP.ROI(img), angle_deg );
imagesc(img_rotated)
axis image off


%% Downscaling an image
x_size_px = 256;
y_size_px = 256;
img_ds = MIP.resize(img,x_size_px,y_size_px);
imagesc(img_ds)
axis image off


%% Upscaling an image
x_size_px = 3000;
y_size_px = 3000;
img_us = MIP.resize(img,x_size_px,y_size_px);
imagesc(img_us)
axis image off


%% Image shift
% nearest pixel
img256 = MIP.ROI(img,255);
x_shift_px = 15;
y_shift_px = -6;
img_shifted1 = MIP.circshift2(img256,x_shift_px,y_shift_px);
imagesc(img_shifted1)
axis image off
%%
% sub-pixel
img_shifted2 = MIP.spshift2(img256,x_shift_px,y_shift_px);
imagesc(img_shifted2)
axis image off



%% Through-Focus vertical shift compensation

v_offsets1 = MIP.offaxis_correction(imgs);
v_offsets2 = MIP.offaxis_correction(21, 900, 0.400e-6, 6);
idx = 4;
imagesc(MIP.ROI(imgs{idx},256,0,round(v_offsets2(idx))))
axis image off



%% Through-Focus data formatting
img_tff = MIP.format_tf(imgs,400,0,0);
idx = 1;
imagesc(img_tff{idx})
title('formatting through focus stacks')
axis image off


%% Finer Through-Focus vertical shift compensation
% identical but with a finer compensation
[bb,dz_m] = MIP.format_tf(imgs,400,0,700,'fine');
idx = mod(1-1,7)+1;
imagesc(bb{idx})
title(sprintf('delta-z = %+1.2f nm',dz_m(idx)*1e9))
axis image off


%% Measuring pitch
imgs = MIP.read('/Users/awojdyla/Documents/MATLAB/SHARP/data/SHARP_2015-06-23_LBNL','IntelCalMask',40);
im = MIP.crop2(MIP.rotate(MIP.ROI(imgs{4},200,0,0),1.2),100);
%MIP.measure_halfpitch(im(:,:),15,'ConsistencyTest')
imagesc(im)
%% Extracting Single line
im = MIP.crop2(MIP.rotate(MIP.ROI(imgs{4},200,0,0),1.2),100);
img_lines = MIP.extract_lines(im);
imagesc(im)
%%
idx = 4;
imagesc(img_lines{mod(idx-1,length(img_lines))+1})
axis image off

%% Extracting LER
idx = 3;
threshold_ct = 100*24;
[ll_nm, lr_nm] = MIP.extract_ler(img_lines{idx}, threshold_ct,1);
subplot(121)
imagesc(img_lines{idx})
axis image
subplot(122)
plot(fliplr(ll_nm),(0:(size(img_lines{idx},1)-1))*15e-3,fliplr(lr_nm),(0:(size(img_lines{idx},1)-1))*15e-3);
ylabel('height [µm]')
xlabel('position [µm]')
title(sprintf('thresold : %04.0f',threshold_ct))
xlim([0 1000*15e-3])