function [Hi,Hr,idx] = genTF_rIDT(grid,obj,scope,iNA)
%%
%genTF function
%
%Purpose: This function performs the transfer function generation for a set
%of input angles based on scalar reflection-mode imaging in a semi-infinite
%boundary condition.
%
%Inputs: grid - structure containing parameters:
%               U,V - MxM matrices holding horizontal, vertical spatial
%                     frequencies of image space, respectively
%        obj - structure containing parameters:
%              RI - 1x2 vector containing refractive index values (RI)
%                   Form: [Imaging medium RI, boundary material RI]
%        scope - structure containing parameters:
%                lambda - scalar of imaging wavelength (um)
%                rngNA  - 1x2 vector containing min. and max. allowed
%                         illumination NA.
%                         Form: [min. illumination NA, max. illum. NA]
%                zVal - scalar defining Microscope focal plane (um)
%                NA - scalar defining objective NA
%                sMag - scalar defining source illumination magnification
%                       (1 if no magnification is used on source size).
%                f - scalar defining objective focal length (um)
%                sPin - scalar defining light source's pre-magnified 
%                       radius (um)
%        iNA - Nx2 matrix containing illumination angle NA values
%              Form: [NAy, NAx], where NAy,x are Nx1 vectors
%Outputs: Hi, Hr = M x M x N 3D matrices containing rIDT absorption and
%                  real Transfer functions, respectively.
%         idx - Nx1 vector defining which images will be used for
%               reconstruction based on defined NA range.
%-------------------------------------------------------------------------%

%--------------------Argument Check and Handle Initialization-------------%

%% Handle Initialization
FT = @(x) fftshift(fft2(ifftshift(x)));
iFT = @(x) fftshift(ifft2(ifftshift(x)));

%% Variable Initialization

%Generate spatial frequency values from illumination angles
iSF = iNA ./ scope.lambda;
iSF(:,3) = sqrt((1/scope.lambda)^2 - (iSF(:,1).^2 + iSF(:,2).^2));

%Calculate radially symmetric illumination NA
NA_rad = sqrt(iNA(:,2).^2 + iNA(:,1).^2);

%Generate binary array for telling whether to include image or not
useAng = zeros(size(iNA,1),1);
useAng(NA_rad >= scope.rngNA(1) & NA_rad <= scope.rngNA(2)) = 1;

%Calculate Incident illumination wavevector
K = ((2 * pi) ./ scope.lambda) * obj.RI(1);

%Calculate obliquity factor
obl = real(sqrt(K.^2 - (2*pi)^2 .* (grid.U.^2 + grid.V.^2)));

%Calculate Radial NA for full U,V Grid
gridNA = scope.lambda .* sqrt(grid.U.^2 + grid.V.^2);

%Calculate Propagation Exponential Factor
prop = exp(1i .* obl .* scope.zVal);

%Calculate fresnel coefficients for entire SF grid
Rv = fresnelCoeff_PC(obj.RI(1),obj.RI(2),asind(gridNA));

% Generate pupil, source function
P = double(circ(sqrt((grid.U./(scope.NA/scope.lambda)).^2 + ...
                (grid.V./(scope.NA/scope.lambda)).^2))); 
rS = scope.sPin*scope.sMag / (scope.lambda * scope.f);
S = double(circ(sqrt((grid.U./rS).^2 + (grid.V./rS).^2)));

%Clip fresnel coefficients to pupil region
RP = P .* Rv;

%Calculate B term for generating Transfer Function
B = conj(P) .* prop .* (1./obl);

%Remove NaN from dividing by zero
B(isnan(B)) = 0;
%-------------------------------------------------------------------------%

%-----------------------------Variable Preallocation----------------------%

%Preallocate Transfer Function Matrices
Hi = zeros(size(grid.U,2),size(grid.V,1), sum(useAng'));
Hr = Hi;

%Define counter for indexing TF within acceptable illumination NA range
cnt = 0;

%Preallocate variable saving image index values used during TF generation
idx = [];
%-------------------------------------------------------------------------%

%--------------------------Generate Transfer Functions--------------------%

%Create wait bar while TF is generated
doneBar = waitbar(0,'Calculating Transfer function');

%Generate Transfer functions in for loop
for k = 1:size(iSF,1)
    %Update waitbar
    waitbar(k/size(iSF,1),doneBar);
    
    %Only make TF for angles within illumination NA range
    if(useAng(k))
        cnt = cnt + 1;
        
        %%%%TF Generation
        
        %Move source mask to current illumination angle
        shift_S = generateSource(S,iSF(k,:),grid.U,grid.V,0);
        
        %Obtain reflection coefficients for current illumination
        RS = shift_S .* Rv;
        
        %Generate A terms for TF convolution
        A = shift_S .* abs(P).^2 .* RS .* P .* conj(prop);
        
        %Generate DC Intensity Value
        DC = sum(sum(RS.^2 .* abs(P).^4 .* shift_S)); 
        
        %Generate Ha TF for current illumination angle
        Hi(:,:,cnt) = iFT(real(FT(RS .* A) .* FT(B))) + iFT(real(FT(A) .* FT(RP .* B)));
        Hi(:,:,cnt) = -(K^2 .* Hi(:,:,cnt))/DC;
        %Generate Hb TF for current illumination angle
        Hr(:,:,cnt) =iFT(1i.*imag(FT(RS.* A).*FT(B))) + iFT(1i.*imag(FT(A).*FT(RP.*B)));
        Hr(:,:,cnt) = -(1i * K^2 * Hr(:,:,cnt))/DC;  
        %Save image index values
        idx = [idx,k];
              
    end %End of use angle check
end %End of illumination spatial freq. loop

%Close Wait bar
close(doneBar);

%-------------------------------------------------------------------------%

end %End of genTF function
