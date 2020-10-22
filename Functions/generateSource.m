function [sOut] = generateSource(sMask,SFval,U,V,flip)
%%
%generateSource Function
%
%Purpose: This is an alternate circular source generation code that aims to
%         generate identical source masks for each illumination angle. This
%         function is basically just a shell for interpolating an input
%         mask to a new position. This function won't work well if you have
%         a mask that has nonzero values near the image edge or if the
%         mask is not 1s and 0s.
%
%Inputs: sMask: 2-D matrix describing the generic mask to be shifted
%        SFval: 1 x 3 vector describing illumination angles to shift mask
%               by
%        U: 2-D matrix of spatial frequencies along horizontal direction
%        V: 2-D matrix of spatial frequencies along vertical direction
%        flip: scalar toggling whether values need to flip along y axis
%Outputs: sOut: 2-D Matrix containing new source mask
%-------------------------------------------------------------------------%
if(nargin < 5)
    flip = 1;
end
iTog.floor = 0; %Floor all resulting pixel values
dfx = abs(U(1,2) - U(1,1));
dfy = abs(V(1,1) - V(2,1));
% dfx = abs(U(1,1) - U(2,1));
% dfy = abs(V(1,2) - V(1,1));
if(iTog.floor == 0)
    if(flip)
        sOut = imtranslate(sMask,[round(SFval(2)./dfx), -round(SFval(1)./dfy)]);
    else
        sOut = imtranslate(sMask,[round(SFval(2)./dfx), round(SFval(1)./dfy)]);
    end
else
    sgns = sign(SFval); %Obtain signs for each term
    if(flip)
        sOut = imtranslate(sMask,[sgns(2) * floor(abs(SFval(2)./dfx)), -sgns(1) * floor(abs(SFval(1)./dfy))]);
    else
        sOut = imtranslate(sMask,[sgns(2)*floor(abs(SFval(2)./dfx)), sgns(1)*floor(abs(SFval(1)/dfy))]);
    end
end
% sOut = interp2(U,V,sMask,U-SFval(2),V-SFval(1));
% sOut(isnan(sOut)) = 0;
% sOut(sOut > 0) = 1; 
end %End of generateSource function