function obj = apply_tikhonov(fIm,Hr,Hi,w)
%%
%apply_tikhonov Function
%
%Purpose: This function performs tikhonov regularization for reconstructing
%         an object from a set of spectral images under different oblique
%         illuminations and their simulated transfer functions
%
%Inputs: fIm - M x M x N matrix of Fourier-transformed normalized images.
%        Hr - M x M x N matrix or real transfer functions.
%
%        Hi - M x M x N matrix of imaginary transfer functions.
%
%        w - 2x1 vector of [real, imaginary] regularization weights for
%            tikhonov regularization.
%
%Outputs: im = MxMx2 reconstructed real and imaginary object component
%
%-------------------------------------------------------------------------%
iFT = @(x) fftshift(ifft2(ifftshift(x))); %Performs rFPM Inverse Fourier Transform

%Calculate square L-2 Norm for Ha, Hb
sumHi = sum(abs(Hi).^2,3);
sumHr = sum(abs(Hr).^2,3);

%Calculate unweighted denominator term
denom_n = sumHi .* sumHr - (sum(conj(Hr) .* Hi,3) .* sum(conj(Hi) .* Hr,3));

%Add tikhonov regularizers to sumHa. sumHb
sumHr = sumHr + w(1);
sumHi = sumHi + w(2);

%Perform convolution of exp. results with TF conjugates in fourier space
aIm = sum(conj(Hi) .* fIm,3);
bIm = sum(conj(Hr) .* fIm,3);

%Calculate common denominator for Image Reconstruction
denom = sumHi .* sumHr - (sum(conj(Hr) .* Hi,3) .* sum(conj(Hi) .* Hr,3));
denom(denom < eps) = 0;

%Perform imaginary object reconstruction
iIm = aIm .* sumHr;
iIm = iIm - (bIm .* sum(conj(Hi) .* Hr,3));
iIm = iIm ./ denom;

%Perform real object reconstruction
rIm = bIm .* sumHi;
rIm = rIm - (aIm .* sum(conj(Hr) .* Hi,3));
rIm = rIm ./ denom;

%Correct divby0 object reconstruction errors due to finite TF support
iIm(isinf(iIm)) = 0;
iIm(isnan(iIm)) = 0;
rIm(isinf(rIm)) = 0;
rIm(isnan(rIm)) = 0;

% Convert to real-space
obj = real(iFT(rIm)) + 1i .* real(iFT(iIm));

end %End of recon_TR function
