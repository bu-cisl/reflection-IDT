function dn = perm2RI(pR,pI,n0)

%%Purpose: This function converts complex permittivities into refractive
%%index values.
%
%Input: pR - MxN matrix of real permittivity values
%       pI = MxN matri of imag. permittivities
%       n0 = assumed imaging medium RI
%Output: dn = complex MxN matrix containing real and imaginary refractive
%             index contrast values.
%
%-------------------------------------------------------------------------%
    nr = sqrt(0.5 * ((n0^2 + pR) + sqrt((n0^2+pR).^2 + pI.^2)));
    ni = pI./(2.*nr);
    dn = nr + 1i .* ni - n0;
end