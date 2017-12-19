# MIP - Matlab Image Processing

## About
MIP is a Matlab class provided for basic image processing and a few optical calculations based on Fourier optics.

It does not require any Matlab toolbox to run (most importantly, a few function equivalent to those in [Image Procesing Toolbox](https://www.mathworks.com/products/image.html) are implemented; it comes handy in cases where your are porting some of your code to a computer that does not have this toolbox and need some basic functions, but it is by no means intended to replace the toolbox itself.)

If you are interested in using open programming languages, check out [JLo](https://github.com/awojdyla/jlo), an equivalent for [Julia](http://julialang.org/), a [really powerful language](http://antoine.wojdyla.fr/blog/2017/12/17/julia-language/) with a syntax very similar to Matlab.

## Usage
To use it, download MIP.m
```git clone https://github.com/awojdyla/mip```, make sure `MIP.m` is in your current folder (or in your path, e.g. `addpath(/Users/awojdyla/Documents/MATLAB/mip)`). Then you can start using it by calling function like this:

```matlab
img_rot = MIP.rotate(img, angle_deg)
```

rotates an image `img` by an angle `angle_deg` in degree.

The code should be reasonably documented, and you can help for most functions by typing in the Command Window (e.g. if you need help on `MIP.extract_ler`):

```matlab
>>help MIP.extract_ler
```

## Things to know
MIP can be used to process images stacks. The image stacks are 1D-cells.

Some functions are not implemented for image stacks, but they can be batch processed using `MIP.batch`, e.g.:

```matlab
img_rot = MIP.batch(img,sprintf('MIP.rotate(x,%c)',angle_deg));
```

where `img` is a cell of images, and `angle_deg` is an angle in degree.

## Examples
(Please refere to `test_mip.m` for many other examples)

+ Rotate and rescale an image
+ Crop and pad an image
+ Cell manipulation for display
+ Background subtraction and detrending
+ Cross-section, centroids
+ Displaying and saving complex data (`.kif` format)
+ Optical propagation using Fourier optics
+ Zernike polynomials and projection
+ Measuring the line edge roughess
+ Fourier ptychography reconstruction
+ Fourier Ring Coefficients (FRC)

## Disclaimer
The software is provided "as is", without warranty of any kind, express or implied.

