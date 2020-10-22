function [imOut,LP] = LPFilter_FT(img)
%%
%LPFilter_FT
%
%Purpose: This function filters the input image by a low-pass filter
%matching the image size. Filtering is done in the Fourier space and can
%contain edge artifacts, but this is not always the case. The 
%
%Inputs: img - NxNxM image stack to be filtered
%       
%Outputs: - imOut - NxNxM normalized image stack without DC Removal! This can
%                   be done separately by subtracting the filtered image mean from the image.
%                   It should effectively be 1.
%           LP - NxNxM low-passed filtered image stack. This should look
%                like slowly varying blobs when plotted and is mostly used
%                for quality control 

%% Handle Declaration
FT = @(x) fftshift(fft2(ifftshift(x))); %Performs rFPM Fourier Transform
iFT = @(x) fftshift(ifft2(ifftshift(x))); %Performs rFPM Inverse Fourier Transform

%% Variable Initialization
%Convert img to double if it's uint16 or uint32, get image stack size
img = double(img); 
nImg = size(img,3);

%Generate filter
H = FT(fspecial('average',[size(img,1),size(img,2)])); %Create LP Filter
weighty = waitbar(0,'Normalizing Image'); %Create waitbar for image processing
%% Filter images
for k = 1:nImg
    waitbar(k/nImg,weighty);
    disp(['Filtering Image ' num2str(k) ' of ' num2str(nImg) '...']);
    
    %Filter, Normalize Image
    bkgnd = iFT(FT(img(:,:,k)).*H);
    imOut(:,:,k) = img(:,:,k)./bkgnd;
    LP(:,:,k) = bkgnd;
    
end %End of image filter loop
close(weighty);

end%End of Function
