%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Example Code for Reflection Intensity Phase Microscopy
%
%Purpose: This script provides an example of the image processing and
%object reconstruction pipeline from the 2020 work "Inverse scattering for 
%reflection intensity phase microscopy" by Matlock, A., Sentenac, A., 
%Chaumet, P., Yi, J., and Tian, L. in Biomedical Optics Express. Please cite
% this publication in future publications using this code. 
%
%Script Structure:
%
% 1. Variable and Handle Declarations and Grid Generation
% 2. Loading, Normalizing, and Background-subtracting Raw Microscope Images
% 3. Illumination Angle Calculation and Transfer Function Generation
% 4. Object Reconstruction
% 5. Data Saving and Plotting
%
% If you have any comments or concerns regarding this code, please email
% Alex Matlock (amatlock@bu.edu) or post an issue to our github repository.
%
%-------------------------------------------------------------------------%

clear all
close all

%--------1. Variable and Handle Declarations and Grid Generation----------%
%% Handle Declarations
FT = @(x) fftshift(fft2(ifftshift(x)));  % Defines FFT
iFT = @(x) fftshift(ifft2(ifftshift(x)));  % Defines iFFT

%% Set data file location, save date, and save folder name
file.dPath = 'Insert data location pathway here';
% file.dPath = 'D:\BU\Programs\Matlab Programs\rIDT\Example Data\Buccal_1';
file.svDate = '201022';
file.svLbl = 'Example';

%% Set microscope, object parameters

%Set distance-dependent parameters
scope.pL = 5.05;  % Camera Pixel Size (um)
scope.f = 20e3;  % Microscope focal length (um)
scope.lambda = 0.53;  % Imaging wavelength (um)
scope.sPin = 150;  % Light Source Radius (um)
scope.zVal = 0;  % Microscope focal plane (um)

% Set distance-independent parameters
scope.Np = [1080 1080];  % Lateral image size in pixels
scope.Mag = 10;  % Microscope Detection Path Magnification
scope.NA = 0.25;  % Objective NA
scope.sMag = 3.3;  % Microscope Source Path Magnification
scope.rngNA = [0 0.21];   % sets allowed NA range[minimum, maximum]
scope.pNA = [0.04 pi/8 0 scope.NA];  % illum. NA [radial NA step, angular NA step, min. NA, max. NA];
scope.DOF = scope.lambda./(scope.NA^2);  % Microscope Depth-of-Field

%Set object processing parameters
obj.RI = [1, 1.52];  % [Imaging medium RI, boundary RI]
obj.tau = [1e3 1e3];  % [Real reg. weight, Imag. Reg. weight]


%% Generate real-space, Fourier space grids for TF Generation
[grid] = makeGrids(scope);

%-------------------------------------------------------------------------%
%----2. Load, Normalize, and Background-subtract Raw Microscope Images----%

%Obtain data from original files
load([file.dPath '\Raw Images_flip.mat']);

%filter images, normalize, remove DC
data = LPFilter_FT(data);
data = data - mean(mean(data));

%-------------------------------------------------------------------------%
%---3. Illumination Angle Calculation and Transfer Function Generation----%

%Generate illumination angles
iNA = calc_illumNA(grid,scope);

%% Generate Transfer Functions
[Hi,Hr,idx] = generate_TF(grid, obj, scope, iNA);

%-------------------------------------------------------------------------%
%-----------------------4. Object Reconstruction--------------------------%

% Reconstruct object phase (permittivity * height)
obj.phase = apply_tikhonov(FT(data(:, :, idx)), Hr, Hi, obj.tau);

% Convert phase to Refractive index contrast
obj.n = perm2RI(real(obj.phase)/scope.DOF, imag(obj.phase)/scope.DOF, obj.RI(1));

%-------------------------------------------------------------------------%
%---------------------5. Data Saving and Plotting-------------------------%
%% Save Reconstruction Results, reconstruction parameters
mkdir([file.dPath '\Processed\' file.svDate '\' file.svLbl]);
save([file.dPath '\Processed\' file.svDate '\' file.svLbl ...
      '\Recon_Tau=' num2str(obj.tau(1)) '_' num2str(obj.tau(2)),'.mat'],...
      'obj','scope','grid','-v7.3');

%% Plot Results
rCent = [652 522];  
rCentSz = [600 600];
rReg = [(rCent(1)-rCentSz(1)/2),(rCent(1)+rCentSz(2)/2-1),...
        (rCent(2)-rCentSz(2)/2),(rCent(2)+rCentSz(2)/2-1)];

figure(1);
subplot(1,2,1);
imagesc(flipud(real(obj.n)));axis image off;
title('\Delta n_{real}');
axis([rReg(1) rReg(2) rReg(3) rReg(4)]);
colormap(gray);
colorbar;
caxis([-0.003 0.003]);
subplot(1,2,2);
imagesc(flipud(imag(obj.n)));axis image off;
title('\Delta n_{imag}');
axis([rReg(1) rReg(2) rReg(3) rReg(4)]);
colormap(gray);
colorbar;
caxis([-0.001 0.003]);

%-------------------------------------------------------------------------%