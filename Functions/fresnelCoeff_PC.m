function [Rav] = fresnelCoeff_PC(n1,n2,aMat)
%%
%fresnelCoeff Function
%
%Purpose: This function calculates the fresnel transmission and reflection
%coefficients for a given incident angle on a surface for a field of s and
%p polarization propagating from medium 1 to medium 2. 
%
%-------------------------------------------------------------------------%
%Calculate transmission angle using snell's law
thetaT = asin((n1./n2) .* sind(aMat));
%find the s-polarization reflection coefficient matrix
Rs = (n1*cosd(aMat) - n2*cos(thetaT))./(n1*cosd(aMat) + n2*cos(thetaT));

%Find p-polarization coefficient matrix
Rp = (n2 * cosd(aMat) - n1 * cos(thetaT))./(n2*cosd(aMat) + n1 * cos(thetaT));

%Find average reflection coefficient
Rav = (Rs - Rp)./2;

Rs(isnan(Rs)) = 0; %Removes NaN due to divby0 errors
end%End of fresnelCoeff function