function [NAOut] = calc_illumNA(grid, scope)
%%
% calc_illumNA function
%
%Purpose: This function generates the expected illumination angles from the
%         rIDT experimental setup under an assumption of a circular
%         illumination grid. Using a provided radial NA stepsize and
%         angular NA stepsize, this code determines the predicted
%         illumination angles and outputs them as a vector of [NAx, NAy]
%         values.
%
%Input:  grid - structure containing parameters:
%           U - MxM matrix of horizontal Spatial Frequencies
%           V - MxM matrix of vertical Spatial Frequencies
%           u,v - 1xM vectors containing horizontal, vertical spatial
%                 frequencies of image, respectively
%        scope - structure containing parameters:
%           pNA - 1x4 vector containing illumination imaging parameters
%                     [radial NA step, angular NA step, min. NA, max. NA]
%
%Output: NAOut - Nx2 matrix containing illumination angle NA values 
%                   Form: [NAy, NAx]
%-------------------------------------------------------------------------%

%% Extract relevant structure variables
U = grid.U;
V = grid.V;
rNA = scope.pNA(3:4);
pNA = scope.pNA(1:2);
lambda = scope.lambda;
center = [(size(V,1)/2+1), (size(U,1)/2+1)];
imsz = size(U);

dSF = [abs(grid.u(2) - grid.u(1)), abs(grid.v(2) - grid.v(1))];

rSF = rNA ./ lambda; %Convert NA range to spatial freq. min and max

%% Calculate radial, angular NA pitch
theta = 0:pNA(2):(2*pi);  % Calculate angular pitch with circular illum. grid
stepR = 1/lambda * (rNA(1):pNA(1):rNA(2));  % Calculate radial NA pitch with circular illum. grid

% Convert illum. coordinates to pixel coordinates
nux = round((stepR' * cos(theta)) ./ dSF(1));  
nuy = round((stepR' * sin(theta)) ./ dSF(2));

% Obtain closest radial spatial frequencies in Fourier space grid to illuminations
stepR_c = [U(1,round(stepR./dSF(1))+center(2))', ...
           V(round(stepR./dSF(2))+center(1),1)]; %Grabs spatial frequencies within U,V grid for illumination

%% Generate mask containing all circular illumination positions

ptGrid = zeros(imsz);
for k = 1:size(stepR_c,1)
  ring_mask = zeros(imsz);
  pos_mask = zeros(imsz);
  if(stepR_c(k) == 0)
      ring_mask(center(1),center(2), k) = 1;
  else
      ring_mask = circ(sqrt((U./(stepR_c(k,1)+dSF(1))).^2 + ...
                            (V./(stepR_c(k,2)+dSF(2))).^2)) - ...
                  circ(sqrt((U./(stepR_c(k,1)-dSF(1))).^2 + ...
                            (V./(stepR_c(k,2)-dSF(2))).^2));    
  end
  pos_mask(nuy(k,:)+center(1), nux(k,:)+center(2)) = 1;
  ptGrid = ptGrid + ring_mask .* pos_mask;
end

% Obtain illumination NA values
allNA = lambda .* (sqrt((U .* ptGrid).^2 + (V .* ptGrid).^2));    %Finds NA at each illumination position to filter out larger than acceptable NA

% Remove illuminations below minimum set illum. NA
if(rNA ~= 0)
    ptGrid(allNA < rNA(1)) = 0;
end
% Remove illuminations above maximum set illum. NA
ptGrid(allNA > rNA(2)) = 0; 

%Flips coordinates to align with rIDT inverted camera plane
U = U'; V = -V'; 

% Convert illuminations into NA index matching rIDT system raster scan
NAOut(:,1) = lambda .* V(find(ptGrid == 1));
NAOut(:,2) = lambda .* U(find(ptGrid == 1)); 


end% End of function generateGrid
