%% Test routine for the MIP class


%% Reading one image
% folder where the images are located :
filename = '/Users/awojdyla/Documents/MATLAB/mip/sample_data/img_40s_s03.png';
im = MIP.read(filename);

%% Read a folder
folder = '/Users/awojdyla/Documents/MATLAB/mip/sample_data/defect_tf';
img = MIP.read(folder);
imagesc(img{10})


%% Extracting a region of interest
roi_size_px = 256;
roi_x_decenter = -1;
roi_y_decenter = -6;

[ img_roi, x_roi, y_roi ] = MIP.ROI( img, roi_size_px, roi_x_decenter, roi_y_decenter);
imagesc(img_roi{1})
axis image off


%% Rotate an image
angle_deg = 0;
img_rotated = MIP.rotate( img, angle_deg );
imagesc(img_rotated{1})
axis image off


%% Downscaling an image
x_size_px = 256;
y_size_px = 256;
img_ds = MIP.resize(img,x_size_px,y_size_px);
imagesc(img_ds{1})
axis image off

%% Upscaling an image
x_size_px = 3000;
y_size_px = 3000;
img_us = MIP.resize(img,x_size_px,y_size_px);
imagesc(img_us{1})
axis image off


%% Image shift
% nearest pixel
img256 = MIP.ROI(img,255);
x_shift_px = 15;
y_shift_px = -6;
img_shifted1 = MIP.circshift2(img256,x_shift_px,y_shift_px);
imagesc(img_shifted1)
axis image off

%% sub-pixel shift
img = MIP.ROI(im,30,0,0);
dx_px = 0.4;
dy_px = 11.2;
img_shifted = MIP.spshift2(img*(1+1i).^2,dx_px,dy_px);
MIP.imagec((img_shifted));


%% Display and image export
% Here are some functions that are useful to display bipolar or complex
% data, and also to export complex data to other languages (HDF5 and KIFC)

%% display a unipolar signal with full dynamic range
% generate a gaussian
[X,Y] = meshgrid(-64:63);
Z = MIP.gaussian(X,0,10).*MIP.gaussian(Y,0,20)+1;

% display in black and white
img_bw_rgb = MIP.image_bw(Z);
image(img_bw_rgb)
title('displaying a gaussian+background'); axis image off;

%% display a bipolar signal, with blue/red hue
% generate a bipolar gaussian derivative
[X,Y] = meshgrid(-64:63);
dZ = diff(MIP.gaussian(X,0,10).*MIP.gaussian(Y,0,20))';

% display in blue-white-red colorscale
img_bwr_rgb = MIP.image_bwr(dZ);
image(img_bwr_rgb)
title('displaying a (bipolar) gaussian derivative'); axis image off;

%% Automatically pick appropriate scale for display
[X,Y] = meshgrid(-64:63);
dZ = diff(MIP.gaussian(X,0,10).*MIP.gaussian(Y,0,20))';
MIP.imagezc(dZ);


%% Display a complex image
[X,Y] = meshgrid(-64:63);
% gaussian amplitude
A = MIP.gaussian(X,0,10).*MIP.gaussian(Y,0,20);
% spiral phase
ephi = exp(1i*atan2(X,Y));
img_rgb = MIP.imagec(A.*ephi);
image(img_rgb); axis image


%% Save any image of the image above
filename = 'test';
MIP.write_png(img_rgb,filename);

%% Load KIFC format
defect = MIP.read_kifc(['sample_data' filesep 'phase_defect.kif']);
MIP.imagec(defect);


%% Feature generation
% Here is a section to demonstrate a few useful features that can be
% generated with MIP

%% Gaussian
% create a spatial scale (say, in meters)
x_m = linspace(-1e-2,1e-2,100);
[X_m,Y_m] = meshgrid(x_m,x_m);
% center of the gaussian
x_center_m = 5e-3;
y_center_m = 3e-3;
% full-width at half max in each direction
x_fwhm_m = 1e-3;
y_fwhm_m = 3e-3;
% generate the gaussian!
G = (MIP.gaussian(X_m,x_center_m,x_fwhm_m).*MIP.gaussian(Y_m,y_center_m,y_fwhm_m))';
imagesc(x_m*1e3,x_m*1e3,G.')
xlabel('x-position [mm]')
ylabel('y-position [mm]')

%% Speckle
% screen size
L_m  = 5e-6;
% pixel size
dx_m = 15e-9;
%spatial scale
x_m = 0:dx_m:L_m;

% EUV wavelength
lambda_m = 13.5e-9;

% mirror roughess
rough_m = 1e-10;

% numerical aperture
NA = 0.33/4; % .33NA wafer side, with 4x demag;

% generate speckle
E = MIP.speckle(length(x_m),2*rough_m/lambda_m);
% filter the speckle through a lens
Efilt = MIP.spatial_filter(E, L_m, lambda_m, NA, 0, 0);
imagesc(abs(Efilt).^2); axis image

% with off-axis illumination:
Efilt_off = MIP.spatial_filter(E, L_m, lambda_m, NA, 0.5, 0);
imagesc(abs(Efilt_off).^2); axis image

% with an added central obscuration of NA 0.03
Efilt_df = MIP.spatial_filter(Efilt, L_m, lambda_m, -0.03);
imagesc(abs(Efilt_df).^2); axis image

%% Binary pseudo-random array
% generate a Binary pseudo-random array (needs a prime number)
bprp = MIP.bprp(47);
imagesc(bprp); axis image
% rescale it with nearest neighbor
bprp_4x = MIP.rescale_nn(bprp,4);
imagesc(bprp_4x); axis image
% showing that the spectrum is very dense:
MIP.imagespec(bprp);

%% Non-redundant arrays
nra = MIP.nra(5,1000);
% hard-coded 6x6 NRAs:
%nra = MIP.nra6(3) 
imagesc(nra); axis image
imagesc(MIP.xcorr2(nra)); axis image


%% Bin an image
x_px = -50:49;
[X_px,Y_px] = meshgrid(x_px);
G =(MIP.gaussian(X_px,0,10).*MIP.gaussian(Y_px,0,10))';
G_bin = Sharp.bin2(G,2,2);
imagesc(G_bin); axis image


%% Shot noise
% generate a gaussian
x_m = linspace(-1e-2,1e-2,100);
[X_m,Y_m] = meshgrid(x_m,x_m);
G =(MIP.gaussian(X_m,0,1e-3).*MIP.gaussian(Y_m,0,1e-3))';

% scale the intensity 
N_phot = 100; % number of photons in the image
I_ph = MIP.shot_noise(G.^2, N_phot);
imagesc(I_ph)

%% Additive white gaussian noise
x_m = linspace(-1e-2,1e-2,100);
[X_m,Y_m] = meshgrid(x_m,x_m);
G =(MIP.gaussian(X_m,0,1e-3).*MIP.gaussian(Y_m,0,1e-3))';

% signal to-noise-ratio
snr = 1;
I_n = MIP.awg_noise(G.^2,snr);
imagesc(I_n); axis image


%% Measurements on images

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

%% Critical dimension

%% Line-width roughness

%% Fourier Ring coeffient

%% Propagation

