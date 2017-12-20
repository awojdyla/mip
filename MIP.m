classdef MIP < handle
    % MIP MATLAB IMAGE PROCESSING tool
    %
    % The software is provided "as is",
    % without warranty of any kind, express or implied
    % A Wojdyla, CXRO/LBNL
    % awojdyla@lbl.gov
    % September 2014 - December 2017
    
    
    properties (Constant)
        c   = 299792458         % speed of light [m/s]
        h   = 6.626070040*1e-34 % Planck constanct [J.s]
        eV  = 1.6021766208*1e-19% elementaty charge x 1V [J.s]
        Ag  = 6.0221408571e23   % Avogadro number [mol-1]
        deg = 180/pi            % degree/rad conversion
        sigma_fhwm = 1/2*sqrt(2*log(2))
        %data_folder = '.';
        %save_folder = '.';
    end
    
    methods(Static)
        
        function img = read(filestring)
        %READ Read a single image or folder,     
        %   MIP.read(filename)
        %   MIP.read(folder)
        %
        % See also MIP.ROI, MIP.ROTATE, MIP.RESIZE
        if exist(filestring,'file')==2
            img = double(imread(filestring));
            
        elseif exist(filestring,'dir')==7
            folder = filestring;
            subarray = @(x) {x(3:end).name};
            if folder(end)~=filesep
                folder = strcat(folder, filesep);
            end
            filenames = subarray(dir(folder));
            
            img = cell(1,length(filenames));
            for i_img=1:length(img)
                img{i_img} = double(imread(strcat(folder,filenames{i_img})));
            end
        end
        end
        
        function [ img_roi, x_roi, y_roi ] = ROI( img, roi_size_px, x_d, y_d)
        %ROI Extract the region of interest from a SHARP image
        %   [img_roi, x_roi, y_roi] = MIP.ROI(img)
        %       uses a (256px)^2 region centered around the ROI
        %   [img_roi, x_roi, y_roi] = MIP.ROI( img, roi_size_px)
        %       lets you define the size of the ROI
        %   [img_roi, x_roi, y_roi] = MIP.ROI( img, roi_size_px, x_d, y_d)
        %       lets you define an offset in the x- and y-direction
            
            if nargin ==1
                roi_size_px = max(256, min(size(img)));
            end
            if nargin <= 2
                x_d = 0;
                y_d = 0;
            elseif nargin <=3
                error('not enough arguments!')
            end
            
            if isa(img,'cell') %recursion
                img_roi = MIP.batch(img, ...
                    sprintf('MIP.ROI(x, %d, %d, %d)',roi_size_px,x_d,y_d));
                y_roi = (1:roi_size_px)+fix((size(img{1},1)-roi_size_px)/2)+y_d;
                x_roi = (1:roi_size_px)+fix((size(img{2},2)-roi_size_px)/2)+x_d;
            else
                % make sure the offset is not too large
                if abs(x_d)>(min(size(img))-roi_size_px)/2 || abs(y_d)>(min(size(img))-roi_size_px)/2
                    error('the positions offset are too large')
                end    

                y_roi = (1:roi_size_px)+fix((size(img,1)-roi_size_px)/2)+y_d;
                x_roi = (1:roi_size_px)+fix((size(img,2)-roi_size_px)/2)+x_d;
                img_roi = img(y_roi,x_roi);
            end

        end
        
        function img_rot = rotate(img, angle_deg, opt_arg)
        %ROTATE Rotate an image using interpolation
        %   img_rotate = MIP.rotate(img,angle_deg)
        %   rotates the image by the specified angle, ccw.
        %   the input can be a single image or a cell stack
        %
        %   img_rotate = MIP.rotate(img,angle_deg,'crop')
        %   does the same but crops the image to remove nan zones
        %
        % See also MIP.RESIZE, MIP.CROP2
            
            %batch mode, if the input is an image stack
            if isa(img,'cell')
                img_rot = MIP.batch(img,sprintf('MIP.rotate(x,%d)',angle_deg));
            else
                
                if isreal(img)
                    [X1,Y1] = meshgrid(linspace(-0.5,0.5,size(img,1)),...
                        linspace(-0.5,0.5,size(img,2)));
                    rot = @(theta) [cos(theta), sin(theta); -sin(theta) cos(theta)];
                    
                    mat_rot = rot(-angle_deg*pi/180);
                    
                    X2 = reshape(X1(:)*mat_rot(1,1) + Y1(:)*mat_rot(1,2),size(img,1),size(img,2));
                    Y2 = reshape(X1(:)*mat_rot(2,1) + Y1(:)*mat_rot(2,2),size(img,1),size(img,2));
                    
                    img_rot = interp2(X1,Y1,img,X2,Y2,'linear');
                else % recursive calls for amplitude and phase;
                    img_rot = abs(MIP.rotate(abs(img),angle_deg)).*...
                        exp(1i.*MIP.rotate(angle(img),angle_deg));
                end
                
                if nargin==3
                    if strcmp(opt_arg,'crop')
                        [~,m1] = min(isnan(img_rot(:,1)));
                        [~,m2] = min(isnan(img_rot(1,:)));
                        img_rot = MIP.crop2(img_rot,...
                            size(img,1)-m1,size(img,2)-m2);
                    else
                        error(sprintf('MIP.rotate > unknown argument :"%s"',opt_arg))
                    end
                end
            end
        end
        
        
        function img_resized = resize(img,x_size_px,y_size_px)
        %RESIZE  Resize an image using interpolation
        %   img_resized = MIP.resize(img, size_px)
        %   img_resized = MIP.resize(img, x_size_px, y_size_px)
            if nargin==2
                if size(img,1)==size(img,2)
                    y_size_px = x_size_px;
                else
                    error('the image should be square if you only use two arguments.')
                end
            end
            if isa(img,'cell')
                img_resized = MIP.batch(img,sprintf('MIP.resize(x,%d,%d)',x_size_px,y_size_px));
            else
                [X1,Y1] = meshgrid(linspace(-0.5,0.5,size(img,1)),linspace(-0.5,0.5,size(img,2)));
                [X2,Y2] = meshgrid(linspace(-0.5,0.5,  x_size_px),linspace(-0.5,0.5,  y_size_px));
                img_resized = interp2(X1,Y1,img',X2,Y2)';
            end
        end
        
        function img_binned = bin2(img,p,q)
        %BIN2 2-dimensionnal binning of an image
        %   img_binned = MIP.bin2(img, x_bin, y_bin)
        %
        % Example :
        %   img_binned = bin2(img, 2, 3);

        % Matt, from MathCentral
        % http://www.mathworks.com/matlabcentral/newsreader/view_thread/248189
            
            [m,n]=size(img); %M is the original matrix
            
            img=sum( reshape(img,p,[]) ,1 );
            img=reshape(img,m/p,[]).'; %Note transpose
            
            img=sum( reshape(img,q,[]) ,1);
            img_binned=reshape(img,n/q,[]).'; %Note transpose
        end
        
        function out = pad2(x, numrows, numcols)
        %PAD2  Pad a matrix with zeros, symmetrically in bothdirections
        %   out = MIP.pad2(in, newlength)
        %   assumes that the final size in both directions is the same
        %   out = MIP.pad2(in, numrows, numcols)
        %   pad with different size in both directions
        %
        %   See also MIP.CROP2

        % R Miyakawa
            
            if nargin == 2
                numcols = numrows;
            end
            [sr, sc] = size(x);
            
            if numrows < sr && numcols < sc
                out = crop2(x,numrows, numcols);
                return;
            end
            
            out = zeros(numrows, numcols);
            startr = ceil((numrows - sr)/2);
            startc = ceil((numcols - sc)/2);
            
            out(startr + 1: startr + sr, startc + 1: startc + sc) = x;
        end
        
        function out = crop2(in, lenr, lenc)
        %CROP2 Crop a matrix, symetrically in both directions
        %   out = MIP.crop2(in, length)
        %   out = MIP.crop2(in, length, width)
        %
        %   See also MIP.PAD2

        % R Miyakawa
            
            if isa(in,'cell')
                out = MIP.batch(in,sprintf('MIP.crop2(x,%d)',lenr));
            else
                
                
                [sr,sc] = size(in);
                
                if nargin<3
                    lenc = lenr;
                end
                
                if lenr > sr && lenc > sc
                    out = pad2(in, lenr, lenc);
                    return;
                end
                
                if nargin == 2
                    lenc = lenr;
                end
                
                wr = lenr/2;
                wc = lenc/2;
                out = zeros(lenr, lenc);
                out(1:lenr, 1:lenc) = ...
                    in( floor(sr/2-wr+1):floor(sr/2+wr), floor(sc/2 - wc+1):floor(sc/2+wc));
            end
        end
        
        function image_shifted = circshift2(input, x_shift_px, y_shift_px)
        %CIRCSHIFT2 Image shift using circular permutation
        %   img_shifted = MIP.circshift2(img, x_shift_px, y_shift_px)
        %   shifts the image in x- and y-direction, to the nearest pixel
        %
        % See also MIP.SPSHIFT2
            
            image_shifted = circshift(circshift(input',round(x_shift_px))',round(y_shift_px));
        end
        
        function img_shifted = spshift2( img,x_shift_px, y_shift_px )
        %SPSHIFT2 Sub-pixel image shift using spectral phase
        %   img_shifted = MIP.spshift2( img,x_shift_px, y_shift_px )
        %       shifts the image in x- and y-direction with fractional pixel
        %
        %    This function does work for complex data., and wraps around
        %
        % See also MIP.CIRCSHIFT2
            
            if size(x_shift_px,1) ~= 1 || size(x_shift_px,2) ~= 1 ...
                    || size(y_shift_px,1) ~= 1 || size(y_shift_px,2) ~= 1
                error('Sparp.spshift :: the shift must be a scalar')
            elseif isempty(x_shift_px) || isempty(x_shift_px)
                error('Sparp.spshift :: empty shift')
            end
            
            fx = MIP.fs(1:size(img,2));
            fy = MIP.fs(1:size(img,1));
            
            [Fx, Fy] = meshgrid(fx, fy);
            
            img_shifted = real(MIP.ift(exp(1i*2*pi*(Fx*x_shift_px+Fy*y_shift_px)).*MIP.ft(real(img))));
            if sum(abs(imag(img(:))))~=0
                img_in =  real(MIP.ift(exp(1i*2*pi*(Fx*x_shift_px+Fy*y_shift_px)).*MIP.ft(imag(img))));
                img_shifted = img_shifted+1i*img_in;
            end
            
        end
        
        function [ img_lowpass ] = lowpass( img, radius )
        %LOWPASS Applies a gaussian lowpass filter to an image
        %   img_lowpass = lowpass2( img, radius )
        %       performs a Gaussian lowpassing with specified radius
            
            if nargin == 1
                radius = 10;
            end
            
            [XX,YY] = meshgrid(...
                (1:size(img,1))-floor(size(img,1)/2),...
                (1:size(img,2))-floor(size(img,2)/2));
            RR = sqrt(XX.^2+YY.^2);
            GAUSSIAN_KERNEL = fftshift(fft2(ifftshift( exp(-(RR/(2*radius)).^2))));
            
            if isa(img,'cell')
                img_lowpass = cell(size(img));
                for i=1:length(img_lowpass)
                    IMG =             fftshift(fft2(ifftshift( img{i} ) ) );
                    img_temp = real(ifftshift(ifft2(fftshift( IMG.*GAUSSIAN_KERNEL))));
                    img_lowpass{i} = img_temp*norm(img{i}(:))/norm(img_temp);
                end
                
            else
                IMG =             fftshift(fft2(ifftshift( img ) ) );
                img_temp = real(ifftshift(ifft2(fftshift( IMG.*GAUSSIAN_KERNEL))));
                img_lowpass = img_temp*norm(img(:))/norm(img_temp);
            end
            
        end

        function [ img_out ] = detrend2( img_in )
        %DETREND2 Remove a 2nd order 2D polynomial trend in an image
        %   [ img_out ] = MIP.detrend2(img_in)
            
            
            [X,Y]=meshgrid(1:size(img_in,2),1:size(img_in,1));
            mask = img_in>0.5*max(img_in(:));
            
            x = X(mask);
            y = Y(mask);
            z = img_in(mask) ;
            
            C = ones(length(x),1) ;
            V = [x.*y y.^2 y C x x.^2] ;
            
            a = V \ z ;
            poly_xy = @(X,Y) a(1)*X.*Y+a(2)*Y.^2+a(3)*Y+a(4)+a(5).*X+a(6)*X.^2;
            
            %img_out = img_in./poly_xy(X,Y).*double((poly_xy(X,Y)>100));
            img_out = img_in./poly_xy(X,Y).*(1+sign(poly_xy(X,Y)-500))*0.5.*max(img_in(:));
            warning('detrending still under testing')
        end
        
        function cell_out = batch(cell_in, instruction)
        %BATCH Batch processing for cells of image
        %   cell_out = batch(cell_in, instruction)
        %   applies the 'instruction' to the cells considered
        %   the instruction is a function handle
        %                   or a string to be evaluated
        %                   (e.g. 'abs(x).^2')
        %
        % See also MIP.mat2cell, MIP.tile_cell
            
            if min(size(cell_in))~=1
                error('MIP.batch : cell arrays must be 1xN. See MIP.batch2')
            end
            cell_out = cell(size(cell_in));
            if isa(instruction,'function_handle')
                for i=1:length(cell_in)
                    cell_out{i} = instruction(cell_in{i});
                end
            else
                for i=1:length(cell_in)
                    x = cell_in{i};
                    cell_out{i} = eval(instruction);
                end
            end
        end
        
        %FIXME : fail grace typechecking
        function img_as_mat = cell2mat(img_as_cell)
        %CELL2MAT Transforms a cell-based stack of images into a 3D matrix
        %   img_as_mat = MIP.cell2mat(img_as_cell)
        %
        % See also MIP.MAT2CELL
            
            N_img = length(img_as_cell);
            img_as_mat = zeros(size(img_as_cell{1},1),...
                size(img_as_cell{1},2),...
                N_img);
            for i=1:N_img
                img_as_mat(:,:,i) = img_as_cell{i};
            end
        end
        
        %FIXME : fail grace typechecking
        function img_as_cell = mat2cell(img_as_mat)
        %MAT2CELL Transforms a 3D matrix into a cell-based stack of images
        %   img_as_cell = MIP.mat2cell(img_as_mat)
        %
        % See also MIP.CELL2MAT
            
            N_img = size(img_as_mat,3);
            img_as_cell = cell(N_img,1);
            for i=1:N_img
                img_as_cell{i} = img_as_mat(:,:,i);
            end
        end
        
        %TODO : automate on series
        function [ img_appended ] = tile(varargin)
        %TILE  Appends images horizontally
        %   img_tiled = MIP.tile(img1,img2,img3,...)
        %   creates a tiled image of all the input images
            
            if nargin == 1 % do nothing
                img_appended = varargin{1};
                
            elseif nargin == 2 % main part
                img1 = varargin{1};
                img2 = varargin{2};
                % resisizing the initial image if needed
                if size(img1,1)<size(img2,1)
                    img_temp = img1;
                    img1 = zeros(size(img2,1), size(img1,2));
                    img1(1:size(img_temp,1),1:size(img_temp,2)) = img_temp;
                end
                img_appended = zeros(size(img1,1), size(img1,2)+size(img2,2));
                img_appended(1:size(img1,1), 1:size(img1,2)) = img1;
                img_appended(1:size(img2,1), (size(img1,2)+1):end) = img2;
            else % recursive
                img_appended = append2(varargin{1},append2(varargin{2:end}));
            end
        end
        
        function tile_q = tile_cell(img_cells)
        %TILE_CELL     Transform a cell array into an image
        %   tile_q = MIP.tile_cell(img_cells)
        %
        % See also : MIP.CELL_STITCH ,MIP.TILE
        
            tile_q = img_cells{1};
            for q=2:length(img_cells)
                tile_q = MIP.tile(tile_q,img_cells{q});
            end
        end
        
        function stack_q = stack_cell(img_cells)
        %STACK_CELL
        %   tile_q = MIP.stack_cell(img_cells)
        %
        % See also : MIP.CELL_STITCH ,MIP.TILE_CELL
        
            N_r = floor(sqrt(length(img_cells)));
            img_stack = cell(N_r,N_r+1);
            for i= 1:N_r
                for j=1:(N_r+1)
                    if (i-1)*N_r+j<length(img_cells)
                        img_stack{i,j} = img_cells{(i-1)*N_r+j};
                    end
                end
                stack_q = MIP.cell_stitch(img_stack);
            end
        end
        
        %FIXME non square images, inline cells
        function img_full = cell_stitch(img_sub)
        %CELL_STITCH Stitches many 2D cells to form a complete 2D image
        %   img_full = MIP.cell_stitch(img_sub)
            
            if ~isempty(img_sub{1,1})
                N_r = size(img_sub,1);
                N_c = size(img_sub,2);
                roi_size_px = size(img_sub{1,1},1);
            else
                roi_size_px = 0;
                for i=1:size(img_sub,1)*size(img_sub,1)
                    roi_size_px = max(roi_size_px, size(img_sub{i},1));
                end
            end
            
            
            if ~(size(img_sub{1,1},1)==size(img_sub{1,1},2))
                error('MIP.cell_stitch :: sub-images must be square')
            end
            
            img_full = zeros(N_r*roi_size_px,N_c*roi_size_px,size(img_sub{1,1},3));
            
            for r = 1:N_r
                for c = 1:N_c
                    if ~isempty(img_sub{r,c})
                        img_full(((r-1)*roi_size_px+1):r*roi_size_px,...
                            ((c-1)*roi_size_px+1):c*roi_size_px,:) ...
                            = img_sub{r,c};
                    end
                end
            end
            
        end
        
        %% Image stats
        
        function img_sum = sum_img(img)
        %SUM_IMG Perfoms the sum of all images in the imag stack
        %   img_sum = sum_img(img)
        
            N_img = length(img);
            img_sum = zeros(size(img{1}));
            for i = 1:N_img
                img_sum = img_sum +img{i}/N_img;
            end
        end
        
        function rms = rms_diff(img_ref,img_comp)
        %RMS_DIFF Average pixel difference
        %   rms = rms_diff(img_ref,img_comp)
        
            rms = sum(abs(img_ref(:)-img_comp(:)).^2)...
                /sum(abs(img_ref(:)).^2);
        end
        
        function rms = nrms_diff(img_ref,img_comp)
        %NRMS_DIFF Normalized Average pixel difference   
        %   nrms = nrms_diff(img_ref,img_comp)
        %
        %   nrms = 0 for two identical image
        %   nrms = 1 for an image compared to background
        %   nrms can be larger than 1 (MIP.nrms_diff(img,-img)=2)
            rms = sqrt(sum((img_ref(:)-img_comp(:)).^2)) ...
                ./sqrt(sum( img_ref(:)             .^2));
        end
        
        %FIXME : check
        function [ array_out ] = circsum( img_in )
        %CIRCSUM Circular sum in 2D
        %   array_out = circsum(img_in)
        %
        % See also MIP.FRC
            
            if size(img_in,1) ~= size(img_in,2)
                error('the input matrix must be square')
            end
            
            N_px = size(img_in,1);
            
            [X,Y] = meshgrid(1:N_px);
            xc = floor(N_px/2)+1;
            yc = floor(N_px/2)+1;
            
            dcirc = zeros(1,floor(N_px/2));
            for i=1:floor(N_px/2)
                domain = ((X-xc).^2+(Y-yc).^2)>=(i-1)^2 & ((X-xc).^2+(Y-yc).^2)<(i)^2;
                dcirc(i) = sum(sum(domain.*(img_in)));
            end

            array_out = dcirc;
            
        end
        
        function [ freq_axis ] = fs(real_axis)
        %FS 1D or 2D Frequency scale, zero-centered in inverse spatial units
        %   freq_Hz = MIP.fs(time_s)
        %
        % See also MIP.FT, MIP.IFT
            
            fs = 1/(real_axis(2)-real_axis(1));
            Nfft = length(real_axis);
            
            df = fs/Nfft;
            freq_axis = (0:df:(fs-df)) - (fs-mod(Nfft,2)*df)/2;
        end
        
        function SIGNAL = ft(signal)
        %FT 1D or 2D zero-centered Fourier Transform
        %   MIP.FT Performs a Fourier transform using optics convention
        %   SIGNAL = MIP.ft(signal)
        %
        %   See also MIP.IFT, MIP.FS
            
            if size(signal,1) == 1 || size(signal,2) == 1
                SIGNAL = fftshift( ifft( ifftshift( signal ) ) );
            else % perform a 2D fourier Transform
                SIGNAL = fftshift( ifft2( ifftshift( signal ) ) );
            end
        end
        
        function signal = ift(SIGNAL)
        %IFT 1D or 2D zero-centered Inverse Fourier Transform
        %   MIP.IFT Performs a Inv Fourier transform using optics convention
        %   signal = MIP.ift(SIGNAL)
        %
        %   See also MIP.FT, MIP.FS
            
            if size(SIGNAL,1) == 1 || size(SIGNAL,2) == 1
                signal = fftshift( fft( ifftshift( SIGNAL ) ) );
            else %perform a 2D fourier Transform
                signal = fftshift( fft2( ifftshift( SIGNAL ) ) );
            end
        end
        
        function [ img_out ] = remove_dc( img_in )
        %REMOVE_DC Removes the DC component of an image (for better display)
        %   img_out = remove_dc(img_in)
        %
        % See also MIP.REMOVE_BG
            
            img_size = size(img_in,1);
            [X,Y] = meshgrid(1:img_size);
            mask = (X == round(1+img_size/2)) | (Y == round(1+img_size/2));
            img_out = img_in;
            img_out(mask) = 0;mean(img_in(:));
        end
        
        function s_out = sgolay(s_in)
        %SGOLAY Savitzky-Golay filtering (5pt window, 3rd order)
        %   s_out = MIP.sgolay(s_in)
        %   (implementation for missing Signal Processing toolbox)
        %
        % See also MIP.SDIFF
        
            s_out = s_in;
            s_out(3:end-2) = 1/35*(...
                -3*(s_in(1:end-4)+s_in(5:end))...
                +12*(s_in(2:end-3)+s_in(4:end-1))...
                +17* s_out(3:end-2));
        end
        
        function s_out = sdiff(s_in)
        %SDIFF Savitzky-Golay differential (5pt window, 3rd order)
        %   s_out = MIP.sdiff(s_in)
        %   (implementation for missing Signal Processing toolbox)
        %
        %   See also MIP.SGOLAY, MIP.SDIFF2
        
            s_out = s_in;
            s_out(3:end-2) = 1/12*(...
                1*(s_in(1:end-4)+s_in(5:end))...
                -8*(s_in(2:end-3)+s_in(4:end-1))...
                +0* s_out(3:end-2));
        end
        
        function s_out = sdiff2(s_in)
        %SDIFF2 Savitzky-Golay second order differential (5pt window, 3rd order)
        %   s_out = MIP.sdiff2(s_in)
        %
        %   See also MIP.SGOLAY, MIP.SDIFF
        
            s_out = s_in;
            s_out(3:end-2) = 1/7*(...
                2*(s_in(1:end-4)+s_in(5:end))...
                -1*(s_in(2:end-3)+s_in(4:end-1))...
                -2* s_out(3:end-2));
        end
        
        function [frc_array, halfbit_threshold]  = frc(img1, img2)
        %FRC Fourier Ring Coefficients between two images
        %   frc_coeffs = MIP.frc(img1, img2)
        %       computes the FRC coefficients between two images
        %   [frc_array, halfbit_threshold] = MIP.frc(img1,img2)
        %
        % See Nieuwenhuizen et al. Nature Methods 10, 6 (2013)
        %  http://doi.org/10.1038/nmeth.2448
        
            frc_array = MIP.circsum(MIP.ft(img1).*conj(MIP.ft(img2)))./...
                ((sqrt(MIP.circsum(abs(MIP.ft(img1)).^2)).*sqrt(MIP.circsum(abs(MIP.ft(img2)).^2))));
            
            n_vox = MIP.circsum(ones(size(img1)));
            halfbit_threshold = (0.2071+1.9102./sqrt(n_vox))./(1.2071+0.9102./sqrt(n_vox));
        end
        
        
        %TODO implement
        function [argout1, argout2] = cross_section(img_data,x1,y1,x2,y2)
        %CROSS_SECTION Computes the cross-section between two points
        %   [data] = MIP.cross_section(img,x1,y1,x2,y2)
        %   [px_scale,data] = MIP.cross_section(img,x1,y1,x2,y2)
            
            if x1==x2 && y1==y2
                error('the two datapoints should be distinct')
            end
            
            data = img_data;
            % angle of the cross-cut
            theta_rad = atan((y2-y1)./(x2-x1));
            N_pts = sqrt((x2-x1).^2 + (y2-y1).^2);
            pts_px = linspace(0,1,floor(N_pts+1));
            
            % intersects
            x_i_0 = pts_px*cos(theta_rad)*(x2-x1);
            y_i_0 = pts_px*sin(theta_rad)*(y2-y1);
            r_i = sqrt((x_i_0).^2+(y_i_0).^2);
            x_i = x_i_0+x1;
            y_i = y_i_0+y1;
            
            % Extracting oblique strids for linear interpolation
            ff = data(mod(floor(y_i)+size(data,1)*floor(x_i-1)-1,size(data,1)*size(data,2))+1);
            cf = data(mod(floor(y_i)+size(data,1)*floor(x_i)-1,size(data,1)*size(data,2))+1);
            fc = data(mod(floor(y_i+1)+size(data,1)*floor(x_i-1)-1,size(data,1)*size(data,2))+1);
            cc = data(mod(floor(y_i+1)+size(data,1)*floor(x_i)-1,size(data,1)*size(data,2))+1);
            
            data = (1-(y_i-floor(y_i))).*...
                ((1-(x_i-floor(x_i))).*ff...
                +(x_i-floor(x_i)).*cf)...
                +(y_i-floor(y_i)).*...
                ((1-(x_i-floor(x_i))).*fc...
                +(x_i-floor(x_i)).*cc);
            
            if nargout == 1
                argout1 = data;
            elseif nargout == 2
                argout1 = r_i;
                argout2 = data;
            end
        end
        
        %FIXME: make more general
        function half_pitch_4x_nm = measure_halfpitch(img,nm_px,opt_arg)
        %MEASURE_HALFPITCH Estimate of the Half-pitch, wafer units (4x)
        %   half_pitch_nm_4x = measure_pitch(img,nm_px)
        %
        %   half_pitch_nm_4x = measure_pitch(img)
        %       assumes .33 4x NA lens data
        %
        % See also MIP.CROSS_SECTION, MIP.EXTRACT_LER
            
            % effective pixel size (default as .33 4xNA measurements)
            if nargin <2
                nm_px = 15;
            end
            % lateral size of the image (for the cross-cut)
            roi_size_px = size(img,2);
            % frequency scale
            f_cpnm = MIP.fs((0:(roi_size_px-1))*nm_px);
            % Fourier Transform of the image
            IMG = abs(MIP.ft(img));
            
            % take a lateral cross-section to locate the first order
            s = IMG(floor(end/2+1),:);
            df_cpnm = f_cpnm(2)-f_cpnm(1);
            
            % Getting coarse estimate of the 1st order frequency
            [~,imax] = max(s.*(f_cpnm>=df_cpnm));
            f0 = f_cpnm(imax);
            
            % Refining the frequency using a 2nd order fit
            pp = polyfit(f_cpnm((imax-1):(imax+1)),s((imax-1):(imax+1)),2);
            f0_fine = -0.5*pp(2)/pp(1);
            % if the coarse estimate does not agree with the fine estimate,
            % something is wrong
            if abs(f0_fine-f0)>df_cpnm
                warning('pitch detection might be inaccurate')
            end
            
            % half-pitch in wafer units
            half_pitch_4x_nm = 1/f0_fine/2;
            
            % consistency test; if the result is not consistent, maybe the
            % region of interest hasn't been properly chosen
            if nargin>2 && strcmp(opt_arg,'ConsistencyTest')
                if size(img,2)>10
                    % iteratively redo with smaller rois
                    fp = 0.5./MIP.measure_halfpitch(img(:,1:end-2),nm_px);
                    fm = 0.5./MIP.measure_halfpitch(img(:,1:end-5),nm_px);
                    
                    if abs(f0_fine-fp)>df_cpnm || abs(f0_fine-fm)>df_cpnm
                        disp('MIP.pitch_measure :: Consistency test not passed')
                    end
                end
            end
        end
        
        function img_lines = extract_lines(img, offset_px)
        %EXTRACT_LINES Extract individual lines from an image for later analysis
        %   img_lines = MIP.extract_lines(img)
        %
        % See also MIP.EXTRACT_LER, MIP.MEASURE_HALFPITCH
            
            if nargin<2
                offset_px = 0;
            end
            
            % effective pixel size
            pxnm = 15;
            % determine the pitch
            halfpitch_4xnm = MIP.measure_halfpitch(img,pxnm);
            % set the pitch size in px
            pitch_px = halfpitch_4xnm*2/pxnm;
            % number of lines contained in the picture
            N_lines = round(size(img,2)/pitch_px);
            % find the maximum oof the center line to center the other lines
            [~,i_offset] = max(sum(img(:,floor((N_lines-1)/2*pitch_px):floor((N_lines+1)/2*pitch_px)),1));
            % if the maxiumum is after the center, make sure the first line
            % is not a cut line
            if i_offset>pitch_px/2
                i_offset = i_offset - pitch_px/2;
            end
            
            i_offset = i_offset-offset_px;
            
            % extract all the lines
            img_lines = cell(1,N_lines-1);
            for i=1:(N_lines-1)
                try
                    img_lines{i} = img(:,(floor(((1+(i-1)*pitch_px):i*pitch_px)-(i_offset)+(pitch_px)/2)));
                catch
                    img_lines{i} = [];
                end
            end
        end
        
        %FIXME make sure the line is vertical and positive
        function [ x1, x2 ] = extract_ler( img_line, threshold, interp_factor, opt_arg)
        %EXTRACT_LER Extract the left and righ edge of a line
        %   [ left_edge, right_edge ] = extract_ler( img_line, threshold, interp_factor )
        %       where img_line is an image containg only one line
        %             threshold is the threshold value
        %             interp_factor the optional interpolation factor
        %
        % This function uses linear interpolation for finding the edges
        % in a robust an efficient way
        %
        % See also MIP.EXTRACT_LINES
            
            if isa(img_line,'cell')
                error('MIP.extract_ler > please provide a valid single image');
            end
            
            
            if threshold>max(img_line(:))
                warning('MIP.extract_ler > threshold higher than the max value');
            end
            
            if threshold<min(img_line(:))
                warning('MIP.extract_ler > threshold lower than the min value');
                if min(img_line(:))>500
                    warning('you might want to check background subtraction');
                end
            end
            
            % size of th image
            x_px = size(img_line,1);
            y_px = size(img_line,2);
            % if some interpolation is needed, resize the imgae (brute force)
            if nargin<3
                interp_factor = 1;
            else
                img_line = MIP.resize(img_line,...
                    x_px*interp_factor,y_px*interp_factor);
            end
            
            % position of the edge
            x1 = zeros(1,x_px*interp_factor);
            x2 = zeros(1,x_px*interp_factor);
            for i=1:x_px*interp_factor
                % for every vertical position, find the threshold cut
                test = img_line(i,:);
                bounds = find(test>threshold);
                % here are the coarse estimate of threshold crossing
                ixl_e = min(bounds);
                ixr_e = max(bounds);
                
                % size increment
                dx = 1;
                % refine the threasold crossing estimate using
                % explicit linear interpolation
                try
                    % left edge
                    if ixl_e>1 %make sure there is a left edge
                        xl = ixl_e-(test(ixl_e)-threshold)/(test(ixl_e)-test(ixl_e-1))*dx;
                    else %otherwise, pupulate missing edge as NaNs
                        xl = NaN;
                    end
                    % right edge
                    if ixr_e<length(test)
                        xr = ixr_e-(test(ixr_e)-threshold)/(test(ixr_e+1)-test(ixr_e))*dx;
                    else
                        xr = NaN;
                    end
                catch err
                    rethrow(err)
                end
                
                % output
                x1(i) = xl;
                x2(i) = xr;
                
                if nargin == 4
                    if strcmp(opt_arg,'detrend')
                        x1 = detrend(x1);
                        x2 = detrend(x2);
                    end
                end
            end
            
            idx_l = round(x_px/2)*2;
            if nargout==0
                plot((1:size(img_line,2))*15-idx_l/2, img_line(round(end/2),:),'k.-',...
                    x1(round(end/2))*15-idx_l/2,interp1([floor(x1(idx_l)),ceil(x1(idx_l))],...
                    [img_line(end/2,floor(x1(idx_l))),img_line(end/2,ceil(x1(idx_l)))],x1(idx_l)),'rx',...
                    x2(round(end/2))*15-idx_l/2,interp1([floor(x2(idx_l)),ceil(x2(idx_l))],...
                    [img_line(end/2,floor(x2(idx_l))),img_line(end/2,ceil(x2(idx_l)))],x2(idx_l)),'rx')
                imagesc(img_line)
                axis  off
                hold on
                plot(x1,1:length(x1),'b',x2',1:length(x2),'r')
                hold off
            end
        end
        
        function [ x_d, I_d, cd_m, hp_m, NILS ] = extract_cd( x_m, I_ct, threshold_ct )
        %EXTRACT_CD extraction of CD from a measurement
        %	[ x_d, I_d, cd_m, hp_m, NILS ] = extract_cd( x_m, I_ct, threshold_ct )
        %
        % it gives you the spatial coordinates at threshold crossing (x_d) and the
        % value at this point (I_d should be threshold!).
        % It also spits out the cd in nm for each line (if the actual spatial scale is given)
        % and the (signed) NILS for each edge.
            
            if size(I_ct,1)>1 && size(I_ct,2)>1
                error('please provide 1D data')
            end
            i_l = logical(abs([diff(double(I_ct>threshold_ct)) 0]));
            i_r = logical(abs([0 diff(double(I_ct>threshold_ct))]));
            x_l = x_m(i_l);
            I_l = I_ct(i_l);
            I_r = I_ct(i_r);
            
            frac = (threshold_ct-I_l)./(I_r-I_l);
            dx_m = x_m(2)-x_m(1);
            %
            x_d = x_l+frac*dx_m;
            I_d= I_l+frac.*(I_r-I_l);
            
            
            temp = diff(x_d);
            cd_m = temp(1:2:end);
            hp_m = mean(temp);
            NILS = interp1(x_m-dx_m/2,[0 mean(cd_m)*diff(log(I_ct))/dx_m],x_d);
            
            if nargout==0
                plot(x_m,abs([0 mean(cd_m)*diff(log(I_ct))/dx_m]),x_d,abs(NILS),'o')
                
            end
        end
        
        function [ lw, xl, xr ] = extract_lw( array_line, threshold)
        %EXTRACT_LW
        %   [ lw, xl, xr ] = extract_lw( array_line, threshold)
        %
        % See also MIP.extract_ler
            
            if isa(array_line,'cell')
                error('MIP.extract_ler > please provide a valid single image');
            end
            
            if threshold>max(array_line(:))
                warning('MIP.extract_ler > threshold higher than the max value');
            end
            
            if threshold<min(array_line(:))
                warning('MIP.extract_ler > threshold lower than the min value');
                if min(array_line(:))>500
                    warning('you might want to check background subtraction');
                end
            end
            
            % for every vertical position, find the threshold cut
            bounds = find(array_line>threshold);
            % here are the coarse estimate of threshold crossing
            ixl_e = min(bounds);
            ixr_e = max(bounds);
            
            % explicit linear interpolation
            try
                % left edge
                if ixl_e>1 %make sure there is a left edge
                    xl = ixl_e-(array_line(ixl_e)-threshold)/(array_line(ixl_e)-array_line(ixl_e-1));
                else %otherwise, pupulate missing edge as NaNs
                    xl = NaN;
                end
                % right edge
                if ixr_e<length(array_line)
                    xr = ixr_e-(array_line(ixr_e)-threshold)/(array_line(ixr_e+1)-array_line(ixr_e));
                else
                    xr = NaN;
                end
            catch err
                rethrow(err)
            end
            
            lw = xr-xl;
        end
        
        % FIXME change freq definition
        function [ ephi ] = speckle( N, amp)
        % SPECKLE Generation of Low Frequency Noise (~speckle)
        %   ephi = MIP.speckle( N, amp )
        %   	generates a pure phase speckle NxN matrix with amplitude amp (in waves)
        %
        % example (generate EUV mirror speckle):
        % 
        %
        %   See also MIP.GAUSSIAN, MIP.SPATIAL_FILTER

            ephi = exp(1i*2*pi*randn(N)*amp);
            %imagesc(rlow); axis image
        end
        
        %FIXME: non-square matrices
        function Efilt = spatial_filter(E, L_m, lambda_m, NA, x_off_npc, y_off_npc)
        %SPATIAL_FILTER Spatial filtering of an image to account for limited aperture
        %   Efilt = spatial_filter(E, L_m,lambda_m, NA)
        %    where E is complex field, L_m the size of the screen 
        %    filters the object spectrum to the NA.
        %   If NA<0, then the NA is filtered from the inside (~darkfield)
        %
        %   Efilt = spatial_filter(E, L_m,lambda_m, NA, x_off_npc, y_off_npc)
        %       allows you to off-center the illumination (useful for structured illumination)
        %
        %
        % % generate speckle at 13.5nm from 1Angstrom roughness over 5um
        % E = MIP.speckle(334,2e-10/13.5e-9); 
        % % filter the speckle through a lens of NA=0.1
        % Efilt = MIP.spatial_filter(E, 5e-6, 13.5e-9, 0.1, 0, 0);
        
        % on-axis imaging
        if ~exist('x_off_npc','var')
            x_off_npc = 0;
        end
        if ~exist('y_off_npc','var')
            y_off_npc = 0;
        end
        
            % create a frequency scale
            f_cpm = MIP.fs(linspace(-L_m/2, L_m/2, size(E,1)));
            [Fx, Fy] = meshgrid(f_cpm);
            
            % find the limit imposed by the lens
            f_max = asin(NA)/lambda_m;
            
            if NA>=0 %filter the lens
                FILTER = double((Fx-f_max*x_off_npc).^2+(Fy-f_max*y_off_npc).^2<f_max.^2);
            else % inner obscuration (~dark-field)
                FILTER = double((Fx-f_max*x_off_npc).^2+(Fy-f_max*y_off_npc).^2>f_max.^2);
            end
            Efilt = MIP.ift(FILTER.*MIP.ft(E));
            
        end
        
        function I_n = shot_noise(I, N_ph)
        %SHOT_NOISE Adds multiplicative noise to an image to emulate shot noise
        %   I_n = MIP.shotnoise(I, N_ph)
        %       where I is the image amplitude, and N_ph the total number
        %       of photons accross the image
        %
        %   Shot noise follows a Poisson distribution, but when N_ph>20 it
        %   can be approximated by a gaussian.
        %
        % See also MIP.AWG_NOISE
       

            I_ph = I./sum(I(:))*N_ph;
            I_n = round(I_ph+sqrt(I_ph).*rand(size(I_ph)));
        end
        
        function I_n = awg_noise(I,snr)
        %AWG_NOISE Adds additive gaussian white noise to emulate camera noise
        %   I_n = awg_noise(I, snr)
        %       where I is the intensity of the NxN image and snr the
        %       linear signal-to-noise ratio.
        %
        % Note that the definition of signal to noise ratio is discutable,
        % since it depends on the object, really.
        % 
        % See also MIP.SHOT_NOISE
        
            noise = abs(randn(size(I)));
            I_n = I+noise.*std(I(:))/std(noise(:))*1/snr;
        end
        
        function ura = bprp(Nx_prime, Ny_prime)
        %BPRP Generate a Binary Pseudo-Radom Pattern 
        % (a.k.a modified Uniformly Redundant Array)
        %   ura = bprp(N_prime) returns a 2D-BPRP of size N_prime
        %
        %    tips: use primes(N) to get primes!
        %
        % see Fenimore & Cannon "Coded aperture imaging 
        %      with uniformly redundant arrays" dx.doi.org/10.1364/AO.17.000337
        %       dx.doi.org/10.1364/AO.28.004344
        %       dx.doi.org/10.1364/OE.22.019803
        %
        % See also MIP.RESCALE_NN
        
        if nargin==1
            Ny_prime = Nx_prime;
        end
        
        if ( ~isprime(Nx_prime) ||  ~isprime(Ny_prime) )
            error('Needs two prime numbers!');
        end
        
        % basic array
        ba = zeros(Nx_prime,Ny_prime);
        
        % K is associated with Nx_prime and M is associated with Ny_prime.
        
        % a simple method to implement the equations is to evaluate mod(x^2,r) for
        % all x from 1 to r. The resulting values give the locations (I) in Cr
        % that contains +1. All other terms in Cr are -1.
        Cr = zeros(1,Nx_prime)-1;
        cr_idx = unique(sort(mod((1:Nx_prime).^2,Nx_prime)))+1;
        Cr(cr_idx) = 1;
        
        Cs = zeros(1,Ny_prime)-1;
        cs_idx = unique(sort(mod((1:Ny_prime).^2,Ny_prime)))+1;
        Cs(cs_idx) = 1;
        
        for ix = 1:Nx_prime
            for jy = 1:Ny_prime
                if ix == 1
                    ba(ix,jy) = 0;
                elseif ( ix ~= 1 && jy == 1 )
                    ba(ix,jy) = 1;
                elseif ( Cr(ix)*Cs(jy) == 1 )
                    ba(ix,jy) = 1;
                else
                    ba(ix,jy) = 0;
                end
            end
        end
        
        %positive array
        pa = ba.*2-1;
        % b(1,1) has to be equal to 1 so that the sidelobes are flater:
        pa(1,1)=1;

        ura = pa;
        end
        
        function nra_n = nra(N, N_try)
        %NRA Generate NxN Non-Redundant Array (NRA) through random search
        % NRA have the nice property of having an autocorrlation function
        % made only of 0's and 1's (each pixel "interfere" with the others only once)
        % see dx.doi.org/10.1364/JOSAA.28.001107
        %
        %   nra_n = nra(N, N_try)
        %
        % This is done through a random search, since I am not aware
        % of a direct method (there are O(2^(NxN)) possibilities)
        % When N>5 consider using MIP.NRA6 to get quick results
        %
        % See also MIP.NRA6
            if N>4
                warning('NRA search can be quite slow')
            end
            
            if nargin==1
                N_try = 65536;
            end
            
            xcorr2 = @(x,y) MIP.xcorr2(x,y);
            
            idx = 1;
            nra = cell(1,100);
            for i = 1:N_try
                ar = round(rand(N));
                
                test =  xcorr2(ar,ar);
                test(N,N)=0;
                
                if max(test)<=1
                    nra{idx}=ar;
                    idx = idx+1;
                end
                
            end
            
            n_max = 0;
            i_max = 1;
            for i = 1:length(nra)
                if ~isempty(nra{i})
                    n_max = max(n_max,sum(nra{i}(:)));
                    i_max = i;
                end
            end
            nra_n = nra{i_max};
            
        end
        
        function nra6(N)
        %NRA6 6x6 Non-Redundant array
        %   MIP.NRA6(N) with 1<=N<=3 returns a 6x6 NRA
        %
        %
        % See also MIP.NRA
        
        % See González & Mejía,
        % "Nonredundant array of apertures to measure the
        %  spatial coherence in two dimensions with only one interferogram"
        %  http://www.docentes.unal.edu.co/ymejiab/docs/josaa01.pdf
        %  dx.doi.org/10.1364/JOSAA.28.001107
            nra1 = ...
                [1,0,0,0,1,1;
                0,0,1,0,0,0;
                0,0,1,0,0,1;
                0,0,0,0,0,0;
                0,1,0,0,0,0;
                1,0,1,0,0,0];
            
            nra2 = ...
                [0,0,1,1,0,0;
                1,0,1,0,0,1;
                0,0,0,0,0,0;
                1,0,0,0,1,0;
                0,0,0,0,0,1;
                0,1,0,0,0,0];
            
            nra3 = ...
                [1,0,0,0,1,0;
                1,0,1,0,0,1;
                0,0,0,0,0,0;
                1,1,0,0,0,0;
                0,0,0,0,1,0;
                0,0,0,1,0,0];
            
            if     N==1
                nra = nra1
            elseif N==2
                nra = nra2
            elseif N==3
                nra = nra3
            else
                error('Please use N<4. NRA is bad for the US of A')
            end
        end

        
        function [img_out] = rescale_nn( img_in, int_factor)
        %RESCALE_NN nearest neighbour rescale with integer factor
        %   [img_out] = rescale_nn( img_in, int_factor)
        %
        % See also MIP.RESCALE, MIP.BIN
        
        if round(int_factor)~=int_factor || int_factor<1
            error('MIP.rescale_nn: please use an integer factor')
        end
        
            N_x = size(img_in,1);
            N_y = size(img_in,2);
            img_out = zeros(N_x*int_factor,N_y*int_factor);
            for ix=1:N_x
                for iy=1:N_y
                    img_out(((ix-1)*int_factor+1):(ix*int_factor),...
                            ((iy-1)*int_factor+1):(iy*int_factor))...
                            = img_in(ix,iy);
                end
            end
        end
        
        function c = xcorr2(a,b)
        %XCORR2 Two-dimensional cross-correlation.
        %   XCORR2(A,B) computes the crosscorrelation of matrices A and B.
        %   XCORR2(A) is the autocorrelation function.
            
        %   Author(s): M. Ullman, 2-6-86
        %   	   J.N. Little, 6-13-88, revised
        %   Copyright 1988-2002 The MathWorks, Inc.
        %   $Revision: 1.9 $  $Date: 2002/11/21 15:47:07 $
            
            if nargin == 1
                b = a;
            end
            
            c = conv2(a, rot90(conj(b),2));
        end


        
        %FIXME make sure legit....
        function [ img_out ] = invert( img_in )
        %INVERT Invert the tone of an image
        %   img_inverse = MIP.invert(img)
            
            M = max(img_in(:));
            m = min(img_in(:));
            img_out = ((img_in-m)/(M-m)-0.5)*(-1)*(M-m)+m;
        end
        
        
        function [x_d,y_d] = register(img1,img2)
        %REGISTER Sub-pixel registration fo two images
        %   [xd, yd] = MIP.register(img1,img2) give the distance between
        %   two images

            if size(img1,1) ~= size(img2,1) && size(img1,2) ~= size(img2,2)
                error('MIP.register : the two images must be the same size')
            end
            
            [a,~] = MIP.dftregistration(fft2(img1),fft2(img2),15);
            x_d = a(4);
            y_d = a(3);
        end
        
        function g = gaussian(x_px,mean_px,fwhm_px)
        %GAUSSIAN Generates a gaussian function
        %   g = gaussian(x_px) generates a gaussian function with zero mean
        %       and 1 unit full-width at half-maximum
        %   g = gaussian(x_px,mean_px,fwhm_px)
            
            if nargin == 1
                mean_px = 0;
                fwhm_px = 1;
            end
            
            sigma_x = fwhm_px/2*sqrt(2*log(2));
            g = exp(-((x_px-mean_px)/(sqrt(2)*sigma_x)).^2);
        end
        
        
        
        
        %%% Display

        function h = imagespec(varargin)
        %IMAGESPEC Displays the direct spectrum intensity of an image
        %   handle = MIP.imagespec(img);
        %   handle = MIP.imagespec(img, gamma);
        %   handle = MIP.imagespec(x, y, img);
        %   handle = MIP.imagespec(x, y, img, gamma);
        %
        %
        %   See also MIP.imagec
        
            gamma = 1;
            if nargin ==1
                img = varargin{1};
                x = MIP.fs(1:size(img,1));
                y = MIP.fs(1:size(img,2));
            elseif nargin >= 3
                x   = varargin{1};
                y   = varargin{2};
                img = varargin{3};
            end
            if nargin==2
                gamma = varargin{2};
            elseif nargin==4
                gamma = varargin{4};
            end
            h = imagesc(x, y, abs(MIP.ft(img).^gamma));
            axis image off
        end
        
        %TODO add examples of use
        function img_rgb = hsv(img)
        %HSV hsv scaling of the data, returned as a rgb matrix of the image
        % (This is useful when one wants to write the png of complex image
        %   img_rgb = MIP.hsv(img)
        %   returns the MxNx3 matrix to a HSV scaling of the data.
        %
        % See also MIP.IMAGEC, MIP.WRITE_PNG
            
            [N_row,N_col] = size(img);
            
            map = @(x,channel) min(6*mod(x-(channel+1)/3,1),1).*...
                max(min(6*(2/3-mod(x-(channel+1)/3,1)),1),0);
            
            img_abs = abs(img(:)-min(img(:)))/((max(img(:))-min(img(:))));
            
            img_rgb(:,:,1)  = reshape(map(img_abs,1),N_row,N_col);
            img_rgb(:,:,2)  = reshape(map(img_abs,2),N_row,N_col);
            img_rgb(:,:,3)  = reshape(map(img_abs,3),N_row,N_col);
        end
        
        
        function img_rgb = hsv2(img)
        %HSV2 Alternative hsv scaling of the image
        %   img_rgb = MIP.hsv2(img)
        %   uses a triple sine-coded hsv mapping for the data
            
            [N_row,N_col] = size(img);
            
            %x=0:0.01:1;
            map = @(x,delta) cos(2*pi*(x-(delta-1)/3)).*(1+sign(cos(2*pi*(x-(delta-1)/3))))/2;
            %plot(x,map(x,1),x,map(x,2),x,map(x,3))
            
            
            img_abs = abs(img(:)-min(img(:)))/(max(img(:))-min(img(:)));
            
            %third dimension
            img_rgb(:,:,1)  = reshape(map(img_abs,3),N_row,N_col);
            img_rgb(:,:,2)  = reshape(map(img_abs,2),N_row,N_col);
            img_rgb(:,:,3)  = reshape(map(img_abs,1),N_row,N_col);
        end
        
        function  jetplot(N)
        %JETPLOT Plot several graphs with a nice colormap
        %   MIP.jetplot( N ) should be placed before the plot
            set(gca,'ColorOrder',jet(N),'NextPlot','ReplaceChildren')
        end
        
        
        function img_rgb = jet(img)
        %JET Returns an image as RGB using Matlab's jet colormap
        % (This is useful when one wants to write the png of on imagesc
        %   img_rgb = jet(img)
        %
        % See also MIP.write_png
            
            [N_row,N_col] = size(img);
            
            map = @(x,channel) max(min(4*(x+1/8),1),0).*...
                max(min(4*(5/8-x),1),0);
            
            map1 = @(x) map(1-x);
            map2 = @(x) map(1-x-1/4);
            map3 = @(x) map(x);
            
            img_abs = abs(img(:)-min(img(:)))/(max(img(:))-min(img(:)));
            
            %third dimension
            img_rgb(:,:,1)  = reshape(map1(img_abs),N_row,N_col);
            img_rgb(:,:,2)  = reshape(map2(img_abs),N_row,N_col);
            img_rgb(:,:,3)  = reshape(map3(img_abs),N_row,N_col);
        end
        
        
        function printmat(matrix,precision)
        %PRINTMAT Print matrices to the console in a copyable form
        %   MIP.printmat(matrix, precision)
        %
        %   This is useful when one wants to share hard coded matrices
        %   for copy paste, e.g MIP.printmat(humps(rand(3)),1) prints:
        %
        %   [+1.9e+01, +1.9e+01, +1.9e+01;...
        %    +2.1e+01, +2.1e+01, +2.1e+01;...
        %    +1.8e+01, +1.8e+01, +1.8e+01]
        %
            
            if nargin<2
                precision = 3;
            end
            format = sprintf('%s%de','%+1.',precision);
            fprintf('\n[')
            
            if size(matrix,1)==1 || size(matrix,2)==1
                for i=1:length(matrix)
                    if i>1
                        fprintf(', ');
                    end
                    fprintf(format,matrix(i));
                end
            else
                for i=1:size(matrix,1)
                    for j=1:size(matrix,2)
                        if i>1 && j==1
                            fprintf(' ');
                        end
                        if j>1
                            fprintf(', ');
                        end
                        fprintf(format,matrix(i));
                    end
                    if i<size(matrix,1)
                        fprintf(';...\n')
                    end
                end
            end
            fprintf(']\n')
        end

        function img_rgb = image_bw(img)
        %IMAGE_BW Returns a zero-centered RBG image
        %
        % See also IMAGE_BWR
            img = img./max(img(:));
            img_rgb(:,:,1)  = img;
            img_rgb(:,:,2)  = img;
            img_rgb(:,:,3)  = img;
        end

        
        function img_rgb = image_bwr(img)
        %IMAGE_BWR   Blue-White-Red color mapping of the data
        %   img_rgb = MIP.image_bwr(img)
        %   returs the MxNx3 matrix corresponding to the colormap
        %
        % See also IMAGE_BW

            img_norm = img./(max(abs(img(:))));
            
            img_n =-img_norm.*double((img_norm<=0));
            img_p = img_norm.*double((img_norm>=0));
            
            %max(img_arg)
            %third dimension
            img_rgb = zeros(size(img,1),size(img,2),3);
            img_rgb(:,:,1)  = (1-img_n);
            img_rgb(:,:,2)  = (1+(-img_p-img_n));
            img_rgb(:,:,3)  = (1-img_p);
        end
        
        
        function img_rgb = imagec( img_complex )
        %IMAGEC Image of complex image with hsv color mapping
        %   imagec(img_complex) displays an image where  amplitude is mapped 
        %   on the value and the phase is mapped on the hue.
        %   img_rgb =imagec(img_complex) returns the corresponding RGB data,
        %   
        %
        %   See also WRITE_PNG, MIP.HSV, MIP.COLORSCALE
            
            [N_row,N_col] = size(img_complex);
            
            map = @(x,channel) min(6*mod(x-(channel+1)/3,1),1).*...
                max(min(6*(2/3-mod(x-(channel+1)/3,1)),1),0);
            
            img_abs = abs(img_complex(:))/(max(abs(img_complex(:))));
            img_arg = angle(img_complex(:))/(2*pi)+0.5;
            
            %max(img_arg)
            %third dimension
            img_rgb(:,:,1)  = reshape(map(img_arg,1).*abs(img_abs),N_row,N_col);
            img_rgb(:,:,2)  = reshape(map(img_arg,2).*abs(img_abs),N_row,N_col);
            img_rgb(:,:,3)  = reshape(map(img_arg,3).*abs(img_abs),N_row,N_col);
            
            if nargout == 0
                image(img_rgb);
                axis image
            end
        end
        
        function img_rgb = imagecc( img_complex, ratio )
        %IMAGECC Image of complex image with hsv color mapping
        %   imagecc( img_complex, ratio ) displays an image where
        %   amplitude is mapped on the value and the
        %   phase is mapped on the hue.
        %
        %   See also MIP.HSV, MIP.COLORSCALE
            
            if nargin == 1
                ratio = 4; %default = quarterwave
            end
            
            [N_row,N_col] = size(img_complex);
            
            map = @(x,channel) min(6*mod(x-(channel+1)/3,1),1).*...
                max(min(6*(2/3-mod(x-(channel+1)/3,1)),1),0);
            
            img_abs = abs(img_complex(:))/(max(abs(img_complex(:))));
            img_arg = angle(img_complex(:))/(2*pi)*ratio+0.5;
            
            %third dimension
            img_rgb(:,:,1)  = reshape(map(img_arg,1).*abs(img_abs),N_row,N_col);
            img_rgb(:,:,2)  = reshape(map(img_arg,2).*abs(img_abs),N_row,N_col);
            img_rgb(:,:,3)  = reshape(map(img_arg,3).*abs(img_abs),N_row,N_col);
            
            if nargout == 0
                image(img_rgb);
            end
            
            axis image
        end
        
        
        %FIXME : axis argument
        function [ h ] = imagezc( img )
        %IMAGEZC Zero-Centered image display
        %   h= imagezc(img)
        %
        %   if the image is unipolar (min(img)>=0), the image is displayed
        %   in black and white, where black=0 and white=max(img)
        %
        %   if the image is bipolar, the image is displayed in the
        %   darkbluered color scheme, with white being the zero level
        %
        % See also MIP.IMAGE_BW, MIP.IMAGE_BWR, MIP.IMAGEC
        
            if min(img(:))>=0
                h = image(MIP.image_bw(img));
            else
                h = image(MIP.image_bwr(img));
            end
        end
        
        function cscale = colorscale(size,initial_phase,flip)
        %COLORSCALE Scaling of the amplitude+phase visaulization
        %   cscale = colorscale(size) outputs an image with the
        %   requested size
        %
        % See also MIP.IMAGEC
        
            if nargin < 1;
                size = 512;
                initial_phase = 0;
                flip = false;
            elseif nargin<2
                initial_phase = 0;
                flip = false;
            elseif nargin<3
                flip = false;
            end
            
            [X,Y]=meshgrid(linspace(-1,1,size));
            
            A = atan2(Y,X);
            R = sqrt(X.^2+Y.^2);
            if flip
                A = -(A);
            end
            
            B = R.*exp(1i*(A+initial_phase));
            mask = R<1;
            
            cscale = MIP.imagec(B.*mask);
        end

        
        function filename = print(idx)
        %PRINT Prints the current figure into the current folder
        %   MIP.print(idx);
        %   prints the current current figure using the unix timestamp
        %   and its index
        %
        % Example :
        % x=-2:0.1:2
        % for i=10
        %   plot(1:10, x.^i)
        %   MIP.print(i)
        % end
            
            if nargin==0
                idx = 0;
            end
            filename = sprintf('%s_%04.0f', datestr(now,'yyyymmddhhMMSS'), idx);
            
            chandles =  findobj(gcf,'type','axes');
            for i=1:length(chandles)
                if length(chandles)==1
                    set(chandles(i),'FontSize',12)
                else
                    set(chandles(i),'FontSize',12)
                end
            end
            
            
            print(filename,'-dpng');
            print(strcat(MIP.save_folder2,filename),'-dpng');
            h=gcf;
            saveas(h,strcat(MIP.save_folder2,'/',filename),'fig')
            saveas(h,strcat(MIP.save_folder2,'/',filename),'pdf')
            
        end
        
        function filename = print2(font_size, linewidth)
        %PRINT2 Prints the current figure into the current folder
        %   MIP.print(idx);
        %   prints the current current figure using the unix timestamp
        %   and its index
            
            filename = sprintf('%s', datestr(now,'yyyymmddhhMMSS'));
            
            chandles =  findobj(gcf,'type','axes');
            for i=1:length(chandles)
                if length(chandles)==1
                    set(chandles(i),'FontSize',font_size)
                else
                    set(chandles(i),'FontSize',font_size)
                end
            end
            
            chandles =  findobj(gcf,'type','line');
            for i=1:length(chandles)
                set(chandles,'LineWidth',linewidth)
            end
            
            print(filename,'-dpng');
            print(strcat(MIP.save_folder2,filename),'-dpng');
            h=gcf;
            saveas(h,strcat(MIP.save_folder2,'/',filename),'fig')
            saveas(h,strcat(MIP.save_folder2,'/',filename),'pdf')
            
        end
        
        %% I/O
        
        %FIXME : overwrite check
        function write_png(img, filename, folder)
        %WRITE_PNG
        %   MIP.write_png(img, filename, folder)
        %   saves the B&W or RGB image to a filename inide a specific folder
        %
        % See also MIP.WRITE_HDF5, MIP.WRITE_KIF
        
            current_folder = pwd;
            
            if nargin == 3
                cd(folder)
            end
            try
                if length(filename)<4 || ~strcmp(filename((end-3):end),'.png')
                    filename = strcat(filename,'.png');
                end
                imwrite(uint16(img/max(img(:))*(2^16-1)),filename,'png')
            catch err
                if nargin ==3
                    cd(current_folder)
                end
                rethrow(err)
            end
            if nargin == 3
                cd(current_folder)
            end
        end
        
        %FIXME Be careful with the size
        function thumbnail_jpg(img,filename)
            if size(img,1)==2048
                img = MIP.resize(img,512);
            end
            imwrite(uint8(img/max(img(:))*(2^8-1)),strcat(filename,'.jpg'),'jpg')
        end
        
        function write_hdf5(img, filename, folder)
        %WRITE_HDF5    Write complex data using HDF5 file format
        %   MIP.write_hdf5(img, filename)
        %   writes the complex image into the current folder
        %
        %   MIP.write_hdf5(img, filename, folder)
        %   writes the complex image into a specific folder
        %
        % See also MIP.READ_HDF5, MIP.WRITE_PNG, MIP.WRITE_KIF

        %going to the destination folder, if needed
            if nargin == 3
                current_folder = pwd;
                cd(folder)
            end
            
            try
                if ~strcmp(filename((end-2):end),'.h5')
                    filename = strcat(filename,'.h5');
                end
                hdf5write(filename,'/real part',real(img),'/imaginary part',imag(img))
                
            catch err
                if nargin ==3
                    cd(current_folder)
                end
                rethrow(err)
            end
            if nargin == 3
                cd(current_folder)
            end
        end
        
        
        function img = read_hdf5(filename)
        %READ_HDF5 Reads complex data encoded in HDF5 file format
        %   img = MIP.read_hdf5(filename)
        %   read the complex image using its filename
        %
        % See also MIP.WRITE_HDF5
            
            %info  = hdf5info(filename);
            try
                img_real = hdf5read(filename,'/real part');
                img_imag = hdf5read(filename,'/imaginary part');
                img = img_real + 1i*img_imag;
            catch err
                rethrow(err)
            end
        end
        
        %overwrite protection
        function k = write_kif( filename, x, y)
        %WRITE_KIF Write Ken-Iacopo Interchange Format files
        %   fid = MIP.write_kif(filename);
        %   x is the real part of the data and y is the imaginary part (if any)
        %
        % See also MIP.READ_KIF, MIP.WRITE_KIFC

        % I Mochi
            
        %converting to a row major environment
        % A wojdyla Sept'15
        
            x=x.';
            y=y.';
            
            if nargin == 2
                re = real(x);
                im = imag(x);
                x = re;
                y = im;
            end
            
            if length(filename)<4 || ~strcmp(filename(end-3:end),'.kif')
                filename = strcat(filename,'.kif');
            end
            
            fid = fopen(filename,'wb');
            fwrite(fid,size(x,1),'uint32');
            fwrite(fid,size(x,2),'uint32');
            
            fwrite(fid,x,'float32');
            fwrite(fid,y,'float32');
            
            k=fclose(fid);
        end
        
        
        function k = write_kifc( filename, aeiphi )
        %  MIP.WRITE_KIFC Write Ken-Iacopo interchange Format files
        %   fid = MIP.write_kifc(filename,img_complex);
        %   x is the real part of the data and y is the imaginary part (if any)
            
            k = MIP.write_kif( filename, real(aeiphi), imag(aeiphi));
        end
        
        function [x,y] = read_kif( filename )
        %READ_KIF Read Ken-Iacopo Interchange format
        %   [x,y] = MIP.read_kif(filename);
        %   x is the real part of the data and y is the imaginary part (if any)
        %
        % See also MIP.READ_KIFC, MIP.WRITE_KIF

        % I Mochi
            
            if length(filename)<4 || ~strcmp(filename(end-3:end),'.kif')
                filename = strcat(filename,'.kif');
            end
            
            fid = fopen(filename,'rb');
            r=fread(fid,1,'uint32');
            c=fread(fid,1,'uint32');
            
            x=fread(fid,[r,c],'float32');
            y=fread(fid,[r,c],'float32');
            
            fclose(fid);
            
            %converting to a column major environment
            % A wojdyla Sept'15
            x=x';
            y=y';
        end
        
        function aeiphi = read_kifc( filename )
        %READ_KIFC Read Ken-Iacopo Interchange format, complex output
        %   [img_complex] = MIP.read_kifc(filename);
        %
        %   See also MIP.WRITE_KIFC, MIP.READ_KIF
            
            [re,im] = MIP.read_kif(filename);
            aeiphi = re+1i*im;
        end

        %TODO document
        function img_resized = resize2(img,x_size_px,y_size_px)
        %RESIZE2  Resize an image using FFT
        %   img_resized = MIP.resize2(img, size_px)
        %   img_resized = MIP.resize2(img, x_size_px, y_size_px)
            
            if nargin==2
                if size(img,1)==size(img,2)
                    y_size_px = x_size_px;
                else
                    error('the image should be square if you only use two arguments.')
                end
            end
            
            if size(img,1)>x_size_px && size(img,2)>y_size_px
                IMG = MIP.ft(img);
                IMG_CROPPED = MIP.crop2(IMG,x_size_px,y_size_px);
                img_resized = real(MIP.ift(IMG_CROPPED));
            else
                IMG = MIP.ft(img);
                IMG_PADDED = MIP.pad2(IMG,x_size_px,y_size_px);
                img_resized = real(MIP.ift(IMG_PADDED));
            end
        end
        
        function [ centro ] = centroid( x,y )
        %CENTROID Centroid computation in one dimension (first moment)
        %   cent = MIP.centroid( x,y )
        %
        % See also MIP.CENTROID2
            
            centro = trapz(x.*y)./trapz(y);
        end
        
        function [ xc_px, yc_px] = centroid2( img )
        %CENTROID2 Centroid computation in two dimension
        %   cent = MIP.centroid2( x,y )
        %
        %   See also MIP.CENTROID

            if ~isa(img,'cell')
                im_tmp = img;
                img = cell(1);
                img{1} = im_tmp;
            end
            
            roi_size_x = size(img{1},1);
            roi_size_y = size(img{1},2);
            x_px = (-roi_size_x/2:roi_size_x/2-1);
            y_px = (-roi_size_y/2:roi_size_y/2-1);
            
            xc_px = zeros(1,length(img));
            yc_px = zeros(1,length(img));
            for i_m =1:length(img)
                img_x = sum(img{i_m},1)/mean(sum(img{i_m},1));
                img_y = sum(img{i_m},2)/mean(sum(img{i_m},2));
                xc_px(i_m) = MIP.centroid(x_px,img_x);
                yc_px(i_m) = MIP.centroid(y_px,img_y');
            end
        end
        
        % Fourier Optics
        
        function u_out=propTF(u_in,L,lambda,z)
        %PROPTF Fourier optics propagation using Transfer Function kernel method
        %   u_out=MIP.propTF(u_in,L,lambda,z)
        %
        %     propagation - transfer function approach
        %     assumes same x and y side lengths and
        %     uniform sampling
        %     u1 - source plane field
        %     L - source and observation plane side length
        %     lambda - wavelength
        %     z - propagation distance
        %     u2 - observation plane field
        %
        % See also MIP.PROPIR, MIP.PROPFF, MIP.TILT, MIP.LENS

        % December 2012
        % Antoine Wojdyla, CXRO/LBNL awojdyla@lbl.gov
        % adapted from 'Computational Fourier Optics' chap V.2
            
            [M,~]=size(u_in); %get input field array size
            dx=L/M; %sample interval
            k=2*pi/lambda; %wavenumber
            if dx<lambda.*z/L
                warning(sprintf('Fourier Transform propagation kernel is not appropriate (by a factor %1.3f).\nConsider using Impulse Response propagation kernel (propIR) to reduce aliasing\n',lambda.*z/(L*dx)))
            end
            
            fx=-1/(2*dx):1/L:1/(2*dx)-1/L; %freq coords
            [FX,FY]=meshgrid(fx,fx);
            
            H=exp(-1i*pi*lambda*z*(FX.^2+FY.^2)); %trans func
            H=fftshift(H); %shift trans func
            U1=fft2(fftshift(u_in)); %shift, fft src field
            U2=H.*U1; %multiply
            u_out=ifftshift(ifft2(U2)); %inv fft, center obs field
        end
        
        
        function u2 = propIR(u1,L,lambda,z)
        %PROPIR Fourier Optics propagation using Impulse Response kernel method
        %
        %     propagation - impulse response approach
        %     assumes same x and y side lengths and
        %     uniform sampling
        %     u1 - source plane field
        %     L - source and observation plane side length
        %     lambda - wavelength
        %     z - propagation distance
        %     u2 - observation plane field
        %
        % See also MIP.PROPTF, MIP.PROPFF, MIP.TILT, MIP.LENS

        %December 2012
        %Antoine Wojdyla, CXRO/LBNL awojdyla@lbl.gov
        %adapted from 'Computational Fourier Optics' chap V.2
            
            [M,~]=size(u1); %get input field array size
            dx=L/M; %sample interval
            if dx>lambda.*z/L
                fprintf('Impulse Response propagation kernel is not appropriate (by a factor %1.3f).\nConsider using Fourier Transform propagation kernel (propTF) to reduce aliasing\n',1/(lambda.*z/(L*dx)))
            end
            k=2*pi/lambda; %wavenumber
            
            x=-L/2:dx:L/2-dx; %spatial coords
            [X,Y]=meshgrid(x,x);
            
            h=1/(1i*lambda*z)*exp(1i*k/(2*z)*(X.^2+Y.^2)); %impulse
            H=fft2(fftshift(h))*dx^2; %create trans func
            U1=fft2(fftshift(u1)); %shift, fft src field
            U2=H.*U1; %multiply
            u2=ifftshift(ifftshift(ifft2(U2))); %inv fft, center obs field
        end
        
        
        function [uout] = tilt(uin,L,lambda,alpha,theta)
        %TILT Tilts the phase modulation equivalent to a tilt
        %   u_out = MIP.tile(u_in,L,lambda,alpha,theta)
        % where u_in is the inital EM field
        %          L is the screen diameter
        %          lambda is the wavelength (in screen units)
        %          alpha is the azimuthal angle(rad)
        %          theta is the other angle (rad)
        %
        % See also MIP.LENS, MIP.PROPTF, MIP.PROPIR, MIP.PROPFF

        % tilt phasefront
        % uniform sampling assumed
        % uin - input field
        % L - side length
        % lambda - wavelength
        % alpha - tilt angle
        % theta - rotation angle (x axis 0)
        % uout - output field
            
            %from 'Computational Fourier Optics' chap VI.1
            %should be alpha<lambda*.(0.5/dx-B), where B is the bandwidth of the source
            
            [M,N]=size(uin); %get input field array size
            dx=L/M; %sample interval
            k=2*pi/lambda; %wavenumber
            
            x=-L/2:dx:L/2-dx; %coords
            [X,Y]=meshgrid(x,x);
            
            uout=uin.*exp(1i*k*(X*cos(theta)+Y*sin(theta)) *tan(alpha)); %apply tilt
        end
        
        
        function u_out = lens(u_in,L,lambda,zf,diam)
        %LENS Creates the phase modulation of a lens
        %   u_out = MIP.lens(u_in,L,lambda,zf,diam)
        %   where u_in is the inital EM field
        %          L is the screen diameter
        %          lambda is the wavelength (in screen units)
        %          zf is the ocal distance  (in screen units)
        %          diam is the diameter of the lens (in screen units)
        %
        % See also MIP.TILT, MIP.PROPTF, MIP.PROPIR, MIP.PROPFF
            
            % should be zf/diam>dx/lambda
            
            % converging or diverging phase-front
            % uniform sampling assumed
            % uin - input field
            % L - side length
            % lambda - wavelength
            % zf - focal distance (+ converge, - diverge)
            % diam - lens diameter
            % uout - output field
            
            [M,N]=size(u_in); %get input field array size
            dx=L/M; %sample interval
            k=2*pi/lambda; %wavenumber
            %
            x=-L/2:dx:L/2-dx; %coords
            [X,Y]=meshgrid(x,x);
            
            u_out=u_in.*exp(-1i*k/(2*zf)*(X.^2+Y.^2)).*double(sqrt(X.^2+Y.^2)<diam/2); %apply focus
        end
        
        function [uq, xq_m] = propHF(up, xp_m, xq_m, lambda_m, z_m, option)
        %PROPFH Huygens-Fresnel propagation
        % u_out = propHF(u_in, x_in_m, x_out_m, lambda, z_m)
        % u_out = propHF(u_in, x_in_m, x_out_m, lambda, z_m, 'off')
        %   removes the normalization
        %
        %   This function can be quite slow, consider MIP.PROPTF for speed
        %
        % See also MIP.PROPTF, MIP.PROPIR
            Np = length(xp_m);
            Nq = length(xq_m);
            uq = zeros(1,Nq);
            
            % propagation
            for iq = 1:Nq
                r_m = sqrt((xq_m(iq)+xp_m).^2+z_m.^2);
                uq(iq) = sum(up.*exp(1i*2*pi*r_m/lambda_m)./r_m);
            end
            
            % normalization
            if nargin<6 || ~strcmp(option,'off')
                uq = uq./sqrt(sum(abs(uq).^2));
                uq = uq.*sqrt(sum(abs(up).^2));
            end
        end
        
        %%%% Aberrations

        function [a, order_x, order_y] = generate_polynomials(order, Np_x, Np_y)
        %GENERATE_POLYNOMIALS Generate 2D polynomials
        %   [a, order_x, order_y] = generate_polynomials(order, Np_x, Np_y)
        %
        %   See also MIP.GENERATE_BASIS, MIP.GENERATE_WEIGHED_POLYNOMIALS
        
            Nxy = ((order + 1) * (order + 2) / 2);
            
            % % Generate a polynomial basis
            index = 1;
            xn = linspace(-1, 1, Np_x);
            yn = linspace(-1, 1, Np_y);
            
            [Xn, Yn] = meshgrid(xn, yn);
            order_x = zeros(1,Nxy);
            order_y = zeros(1,Nxy);
            a = zeros(Np_y, Np_x, Nxy);
            for i=0:order
                for q=0:i
                    
                    ix = i-q;
                    order_x(index) = ix;
                    
                    iy = q;
                    order_y(index) = iy;
                    
                    a(:,:,index) = power(Xn, ix).*power(Yn, iy);
                    index = index + 1;
                end
            end
        end
        
        function weight_ngauss = generate_gaussian_weight(X_m, Y_m, sigma_x, sigma_y)
        %GENERATE_GAUSSIAN_WEIGHT Generate a HHLO beam amplitude profile
        %   weight_ngauss = generate_gaussian_weight(X_m, Y_m)
        %       generates a standard HHLO profile
        %   weight_ngauss = generate_gaussian_weight(X_m, Y_m, sigma_x, sigma_y)
        %       generates a gaussian profile with specified variance
        %
        %   See also MIP.GENERATE_BASIS
        
            if ~exist('sigma_x',1)
                sigma_x = 1;
            end
            if ~exist('sigma_y',1)
                sigma_y = 1;
            end
            mu_x = 0;
            mu_y = 0;
            
            % Define a gaussian
            gauss = exp(-((X_m-mu_x)./(sqrt(2)*sigma_x)).^2-((Y_m-mu_y)./(sqrt(2)*sigma_y)).^2);
            % take the square root (intensity->amplitude)
            gauss_weight = sqrt(gauss);
            % normalize to get a proper weight definition
            weight_ngauss = gauss_weight./sum(gauss_weight(:));
        end
        
        
        function w_basis = generate_basis(polynomials, weight)
        %GENERATE_BASIS Generates an orthonormal basis from basis function
        %   basis = generate_basis(polynomials)
        %       creates an orthogonal basis with unitary weighting
        %   w_basis = generate_basis(polynomials, weight)
        %       creates an orthogonal basis using defined weighting 
        %
        % See also MIP.PROJECT_ON_BASIS, MIP.CHECK_ORTHOGONALITY,
        %          MIP.GENERATE_POLYNOMIALS

            if nargin<2
                weight = 1;
            end
            
            %orthonormalization
            Nxy = size(polynomials,3);
            w_basis = polynomials;
            w_basis(:,:,1) = w_basis(:,:,1)./sqrt(sum(sum(weight.*w_basis(:,:,1).* w_basis(:,:,1))));
            for i = 2:Nxy
                for j=1:(i-1)
                    %orthonality
                    w_basis(:,:,i) = w_basis(:,:,i)-w_basis(:,:,j).*sum(sum(weight.*w_basis(:,:,i).*w_basis(:,:,j)))./sum(sum(weight.*w_basis(:,:,j).* w_basis(:,:,j)));
                    %unitarity
                    w_basis(:,:,i) = w_basis(:,:,i)./sqrt(sum(sum(weight.*w_basis(:, :, i).*w_basis(:,:,i))));
                end
            end
        end
        
        function projection = project_on_basis(vector, basis, weight)
        %PROJECT_ON_BASIS Projects a vector on an (orthonormal) basis
        %   projection = project_on_basis(vector, basis)
        %   projection = project_on_basis(vector, basis, weight)   
        %       does the same with defined weighing
        %
        % See also MIP.GENERATE_BASIS, MIP.CHECK_ORTHOGONALITY,
        %          MIP.GENERATE_POLYNOMIALS
        
            if nargin<3
                weight = 1;
            end
            
            Nxy = size(basis,ndims(basis));
            projection = zeros(1,Nxy);
            for i = 1:Nxy
                projection(i) = sum(sum(weight.*vector.*basis(:,:,i)));
            end
        end
        
        function auto_product = check_orthogonality(basis, weight)
        %CHECK_ORTHOGONALITY Checks the orthogonality of the basis
        %   auto_product = check_orthogonality(basis)
        %   auto_product = check_orthogonality(basis, weight)
        %
        % See also MIP.GENERATE_BASIS, MIP.PROJECT_ON_BASIS,
        %          MIP.GENERATE_POLYNOMIALS
        
            if nargin<2
                weight = 1;
            end
            
            Nxy = size(basis,ndims(basis));
            auto_product = zeros(Nxy);
            for i = 1:Nxy
                auto_product(:,i) = MIP.project_on_basis(weight.*basis(:,:,i),basis);
            end
        end
        
        function w_polynomials = generate_weighed_polynomials(polynomials, weight)
        %GENERATE_WEIGHED_POLYNOMIALS Generate polynomials with a weight
        %   w_polynomials = generate_weighed_polynomials(polynomials, weight)
        % 
        %   See also MIP.GENERATE_POLYNOMIALS, MIP.GENERATE_BASIS, MIP.PROJECT_ON_BASIS,     
            
            Np_y = size(polynomials,1);
            Np_x = size(polynomials,2);
            Nxy  = size(polynomials,3);
            w_polynomials = zeros(Np_y, Np_x, Nxy);
            for i=1:Nxy
                w_polynomials(:,:,i) = polynomials(:,:,i).*weight;
            end
            
        end
        
        function opde_wave = opde(aberration_rad, weight)
        %OPDE Computes the weighed RMS Optical path difference error
        %   opde_wave = opde(aberration_rad, weight)
        
            weight = weight(:)./sum(weight(:));
            opde_wave=sqrt(sum(weight(:).*(aberration_rad(:)-sum(weight(:).*aberration_rad(:))).^2));
        end
        
        function [coeff, polyquad] = best_fit(X_m, Y_m, aberration_rad, weight)
        %BEST_FIT Best 2D quadratic fit (X^2) of the wavefront
        %   [coeff, polyquad] = best_fit(X_m, Y_m, aberration_rad)
        %       performs the fit with uniform weighing
        %   [coeff, polyquad] = best_fit(X_m, Y_m, aberration_rad, weight)
        %       performs the fit with defined weighing
       
            x = X_m(:);
            y = Y_m(:);
            w = weight(:);
            
            polynomial = [w.*x.^0, w.*x.^2.*y.^2];
        
            coeffs =  polynomial \ (aberration_rad(:).*weight(:));
            coeff  = coeffs(2);
            polyquad = X_m.^2;
        end
        
        %%% image measurements
        
        function [fwhm_px, xl, xr] = fwhm( signal, thr )
        %FHWM determines the FWHM
        %   [fwhm_px, xl, xr] = MIP.fwhm( signal, thr )

            % error if 2D
            if (size(signal,1)==1 && size(signal,2)==1)
                error('data has to be 1D')
            end
            
            % reshaping if needed
            if size(signal,1)>size(signal,2)
                signal = signal';
            end
            
            % allow a different threshold
            if nargin==1
                thr = max(signal(:)/2);
            end
            
            bounds = find(signal>thr);
            % here are the coarse estimate of threshold crossing
            ixl_e = min(bounds);
            ixr_e = max(bounds);
            
            % refine the threasold crossing estimate using
            % explicit linear interpolation
            % left edge
            if ixl_e>1 %make sure there is a left edge
                xl = ixl_e-(signal(ixl_e)-thr)/(signal(ixl_e)-signal(ixl_e-1));
            else %otherwise, pupulate missing edge as NaNs
                xl = NaN;
            end
            % right edge
            if ixr_e<length(signal)
                xr = ixr_e-(signal(ixr_e)-thr)/(signal(ixr_e+1)-signal(ixr_e));
            else
                xr = NaN;
            end
            
            fwhm_px = abs(xr-xl);
        end

        function [enc, x_px, fraction] = enc_energy(Intensity, size_px)
        %ENC_ENERGY Encircled energy in one dimension
        %   enc = enc_energy(intensity_1D)
        %   [enc, x_px, fraction] = enc_energy(Intensity, size_px)
            
            if size(Intensity,1)>size(Intensity,2)
                I = Intensity';
            else
                I = Intensity;
            end
            
            Npx = length(I);
            
            % if the signal is center-symmetric
            if mod(Npx,2)==0
                enc = [0 cumsum(I((Npx/2+1):Npx)+I((Npx/2):-1:1))]./sum(I);
            else
                enc_left  = ( cumsum(I(((Npx+1)/2): 1:(Npx-1)))+ cumsum( I(((Npx+1)/2+1): 1: Npx)) )/2;
                enc_right = ( cumsum(I(((Npx+1)/2):-1:      2))+ cumsum( I(((Npx+1)/2-1):-1:   1)) )/2;
                enc = [0 enc_left+enc_right];
                enc = enc./sum(I);
            end
            
            %radius
            x_px = 0:(length(enc)-1);
            
            fraction = 1;
            if nargin==2 && size_px<x_px(end)
                fraction  = interp1(x_px,enc,size_px);
            end
        end

        
        function [output, Greg] = dftregistration(buf1ft,buf2ft,usfac)
        % function [output Greg] = dftregistration(buf1ft,buf2ft,usfac);
        % Efficient subpixel image registration by crosscorrelation. This code
        % gives the same precision as the FFT upsampled cross correlation in a
        % small fraction of the computation time and with reduced memory
        % requirements. It obtains an initial estimate of the crosscorrelation peak
        % by an FFT and then refines the shift estimation by upsampling the DFT
        % only in a small neighborhood of that estimate by means of a
        % matrix-multiply DFT. With this procedure all the image points are used to
        % compute the upsampled crosscorrelation.
        % Manuel Guizar - Dec 13, 2007
            
        % Portions of this code were taken from code written by Ann M. Kowalczyk
        % and James R. Fienup.
        % J.R. Fienup and A.M. Kowalczyk, "Phase retrieval for a complex-valued
        % object by using a low-resolution image," J. Opt. Soc. Am. A 7, 450-458
        % (1990).

        % Citation for this algorithm:
        % Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup,
        % "Efficient subpixel image registration algorithms," Opt. Lett. 33,
        % 156-158 (2008).

        % Inputs
        % buf1ft    Fourier transform of reference image,
        %           DC in (1,1)   [DO NOT FFTSHIFT]
        % buf2ft    Fourier transform of image to register,
        %           DC in (1,1) [DO NOT FFTSHIFT]
        % usfac     Upsampling factor (integer). Images will be registered to
        %           within 1/usfac of a pixel. For example usfac = 20 means the
        %           images will be registered within 1/20 of a pixel. (default = 1)

        % Outputs
        % output =  [error,diffphase,net_row_shift,net_col_shift]
        % error     Translation invariant normalized RMS error between f and g
        % diffphase     Global phase difference between the two images (should be
        %               zero if images are non-negative).
        % net_row_shift net_col_shift   Pixel shifts between images
        % Greg      (Optional) Fourier transform of registered version of buf2ft,
        %           the global phase difference is compensated for.

            % Default usfac to 1
            if exist('usfac')~=1, usfac=1; end
            
            % Compute error for no pixel shift
            if usfac == 0,
                CCmax = sum(sum(buf1ft.*conj(buf2ft)));
                rfzero = sum(abs(buf1ft(:)).^2);
                rgzero = sum(abs(buf2ft(:)).^2);
                error = 1.0 - CCmax.*conj(CCmax)/(rgzero*rfzero);
                error = sqrt(abs(error));
                diffphase=atan2(imag(CCmax),real(CCmax));
                output=[error,diffphase];
                
                % Whole-pixel shift - Compute crosscorrelation by an IFFT and locate the
                % peak
            elseif usfac == 1,
                [m,n]=size(buf1ft);
                CC = ifft2(buf1ft.*conj(buf2ft));
                [max1,loc1] = max(CC);
                [max2,loc2] = max(max1);
                rloc=loc1(loc2);
                cloc=loc2;
                CCmax=CC(rloc,cloc);
                rfzero = sum(abs(buf1ft(:)).^2)/(m*n);
                rgzero = sum(abs(buf2ft(:)).^2)/(m*n);
                error = 1.0 - CCmax.*conj(CCmax)/(rgzero(1,1)*rfzero(1,1));
                error = sqrt(abs(error));
                diffphase=atan2(imag(CCmax),real(CCmax));
                md2 = fix(m/2);
                nd2 = fix(n/2);
                if rloc > md2
                    row_shift = rloc - m - 1;
                else
                    row_shift = rloc - 1;
                end
                
                if cloc > nd2
                    col_shift = cloc - n - 1;
                else
                    col_shift = cloc - 1;
                end
                output=[error,diffphase,row_shift,col_shift];
                
                % Partial-pixel shift
            else
                
                % First upsample by a factor of 2 to obtain initial estimate
                % Embed Fourier data in a 2x larger array
                [m,n]=size(buf1ft);
                mlarge=m*2;
                nlarge=n*2;
                CC=zeros(mlarge,nlarge);
                CC(m+1-fix(m/2):m+1+fix((m-1)/2),n+1-fix(n/2):n+1+fix((n-1)/2)) = ...
                    fftshift(buf1ft).*conj(fftshift(buf2ft));
                
                % Compute crosscorrelation and locate the peak
                CC = ifft2(ifftshift(CC)); % Calculate cross-correlation
                [max1,loc1] = max(CC);
                [max2,loc2] = max(max1);
                rloc=loc1(loc2);cloc=loc2;
                CCmax=CC(rloc,cloc);
                
                % Obtain shift in original pixel grid from the position of the
                % crosscorrelation peak
                [m,n] = size(CC); md2 = fix(m/2); nd2 = fix(n/2);
                if rloc > md2
                    row_shift = rloc - m - 1;
                else
                    row_shift = rloc - 1;
                end
                if cloc > nd2
                    col_shift = cloc - n - 1;
                else
                    col_shift = cloc - 1;
                end
                row_shift=row_shift/2;
                col_shift=col_shift/2;
                
                % If upsampling > 2, then refine estimate with matrix multiply DFT
                if usfac > 2,
                    %%% DFT computation %%%
                    % Initial shift estimate in upsampled grid
                    row_shift = round(row_shift*usfac)/usfac;
                    col_shift = round(col_shift*usfac)/usfac;
                    dftshift = fix(ceil(usfac*1.5)/2); %% Center of output array at dftshift+1
                    % Matrix multiply DFT around the current shift estimate
                    CC = conj(MIP.dftups(buf2ft.*conj(buf1ft),ceil(usfac*1.5),ceil(usfac*1.5),usfac,...
                        dftshift-row_shift*usfac,dftshift-col_shift*usfac))/(md2*nd2*usfac^2);
                    % Locate maximum and map back to original pixel grid
                    [max1,loc1] = max(CC);
                    [max2,loc2] = max(max1);
                    rloc = loc1(loc2); cloc = loc2;
                    CCmax = CC(rloc,cloc);
                    rg00 = MIP.dftups(buf1ft.*conj(buf1ft),1,1,usfac)/(md2*nd2*usfac^2);
                    rf00 = MIP.dftups(buf2ft.*conj(buf2ft),1,1,usfac)/(md2*nd2*usfac^2);
                    rloc = rloc - dftshift - 1;
                    cloc = cloc - dftshift - 1;
                    row_shift = row_shift + rloc/usfac;
                    col_shift = col_shift + cloc/usfac;
                    
                    % If upsampling = 2, no additional pixel shift refinement
                else
                    rg00 = sum(sum( buf1ft.*conj(buf1ft) ))/m/n;
                    rf00 = sum(sum( buf2ft.*conj(buf2ft) ))/m/n;
                end
                error = 1.0 - CCmax.*conj(CCmax)/(rg00*rf00);
                error = sqrt(abs(error));
                diffphase=atan2(imag(CCmax),real(CCmax));
                % If its only one row or column the shift along that dimension has no
                % effect. We set to zero.
                if md2 == 1,
                    row_shift = 0;
                end
                if nd2 == 1,
                    col_shift = 0;
                end
                output=[error,diffphase,row_shift,col_shift];
            end
            
            % Compute registered version of buf2ft
            if (nargout > 1)&&(usfac > 0),
                [nr,nc]=size(buf2ft);
                Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
                Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
                [Nc,Nr] = meshgrid(Nc,Nr);
                Greg = buf2ft.*exp(1i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc));
                Greg = Greg*exp(1i*diffphase);
            elseif (nargout > 1)&&(usfac == 0)
                Greg = buf2ft*exp(1i*diffphase);
            end
            return
        end
        
        function out=dftups(in,nor,noc,usfac,roff,coff)
        % function out=dftups(in,nor,noc,usfac,roff,coff);
        % Upsampled DFT by matrix multiplies, can compute an upsampled DFT in just
        % a small region.
        % usfac         Upsampling factor (default usfac = 1)
        % [nor,noc]     Number of pixels in the output upsampled DFT, in
        %               units of upsampled pixels (default = size(in))
        % roff, coff    Row and column offsets, allow to shift the output array to
        %               a region of interest on the DFT (default = 0)
        % Recieves DC in upper left corner, image center must be in (1,1)
        % Manuel Guizar - Dec 13, 2007
        % Modified from dftus, by J.R. Fienup 7/31/06

        % This code is intended to provide the same result as if the following
        % operations were performed
        %   - Embed the array "in" in an array that is usfac times larger in each
        %     dimension. ifftshift to bring the center of the image to (1,1).
        %   - Take the FFT of the larger array
        %   - Extract an [nor, noc] region of the result. Starting with the
        %     [roff+1 coff+1] element.

        % It achieves this result by computing the DFT in the output array without
        % the need to zeropad. Much faster and memory efficient than the
        % zero-padded FFT approach if [nor noc] are much smaller than [nr*usfac nc*usfac]
            
            [nr,nc]=size(in);
            % Set defaults
            if exist('roff')~=1, roff=0; end
            if exist('coff')~=1, coff=0; end
            if exist('usfac')~=1, usfac=1; end
            if exist('noc')~=1, noc=nc; end
            if exist('nor')~=1, nor=nr; end
            % Compute kernels and obtain DFT by matrix products
            kernc=exp((-1i*2*pi/(nc*usfac))*( ifftshift([0:nc-1]).' - floor(nc/2) )*( [0:noc-1] - coff ));
            kernr=exp((-1i*2*pi/(nr*usfac))*( [0:nor-1].' - roff )*( ifftshift([0:nr-1]) - floor(nr/2)  ));
            out=kernr*in*kernc;
            return
        end
    end
    
    %%TODO
    %Bossung plot
    
    
    methods(Hidden)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%HANDLE CLASS METHODS THAT SHOULD BE HIDDEN TO MAKE
        %%AUTO-COMPLETION EASIER
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function lh = addlistener(varargin)
            lh = addlistener@handle(varargin{:});
        end
        function notify(varargin)
            notify@handle(varargin{:});
        end
        function delete(varargin)
            delete@handle(varargin{:});
        end
        function Hmatch = findobj(varargin)
            Hmatch = findobj@handle(varargin{:});
        end
        function p = findprop(varargin)
            p = findprop@handle(varargin{:});
        end
        function TF = eq(varargin)
            TF = eq@handle(varargin{:});
        end
        function TF = ne(varargin)
            TF = ne@handle(varargin{:});
        end
        function TF = lt(varargin)
            TF = lt@handle(varargin{:});
        end
        function TF = le(varargin)
            TF = le@handle(varargin{:});
        end
        function TF = gt(varargin)
            TF = gt@handle(varargin{:});
        end
        function TF = ge(varargin)
            TF = ge@handle(varargin{:});
        end
    end
end