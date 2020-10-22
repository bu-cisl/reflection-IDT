function [NAOut,polOut] = genIllum_rIDT(grid,scope)
%%

%% Extract relevant structure variables
U = grid.U;
V = grid.V;
rNA = scope.pNA(3:4);
pNA = scope.pNA(1:2);
lambda = scope.lambda;
gType = 'circle';
%-------------------------Convert NA to spatial frequency-----------------%
dSF = [abs(U(1,2) - U(1,1)), abs(V(1,1) - V(2,1))]; %Finds U,V grid spatial frequency step size
rSF = rNA ./ lambda; %Convert NA range to spatial freq. min and max
%Checks if square grid is requested 
if(strcmp(lower(gType),'square') || strcmp(lower(gType),'sqr') || strcmp(lower(gType),'squarena'))
    dSF_NA = pNA ./ lambda; %Convert NA pitch into spatial frequency values

    %Round NA SF to be equal integer number of pixels
    dSF_NA = round(dSF_NA ./ dSF); %rounds pixel count to next integer, converts back into spatial frequencies
%-------------------------------------------------------------------------%
%-----------------------Generate Source Centerpoint Grid------------------%
    cntr = [size(U,2)/2+1 size(V,1)/2+1]; %Note: This only works with 
    uPts = [fliplr(U(1,(cntr(1)-dSF_NA(1)):-dSF_NA(1):1)) U(1,cntr(1):dSF_NA(1):end)]; %Identifies horizontal spatial freq. at all possible pitch positions
    vPts = [fliplr((V((cntr(2)-dSF_NA(2)):-dSF_NA(2):1,1))') (V(cntr(2):dSF_NA(2):end,1))']; %Finds vertial spatial freq. at all possible pitch positions

    ptGrid = checkVals(U,uPts) .* checkVals(V,vPts);
    uOut = U .* ptGrid;
    vOut = V .* ptGrid;
    allNA = lambda .* (sqrt(uOut.^2 + vOut.^2));

    if(rNA ~= 0)
        ptGrid(allNA < rNA(1)) = 0;
    %     allNA(allNA < rNA(1)) = 0;
    end
    ptGrid(allNA > rNA(2)) = 0;
    % allNA(allNA > rNA(2)) =0;

%     axis([850 1200 850 1200]);
    NAOut(:,2) = lambda .* U(find(ptGrid == 1)); %Obtain horizontal source spatial freq. positions
    NAOut(:,1) = lambda .* V(find(ptGrid == 1)); %Obtain vertical source spatial freq. positions

%-------------------------------------------------------------------------%
elseif (strcmp(lower(gType),'circle') || strcmp(lower(gType),'circ') || strcmp(lower(gType),'radial'))
    theta = 0:pNA(2):(2*pi);    %Generates all angles based on input angle pitch
    rPitch = pNA(1) ./ lambda;   %Convert radial pitch into spatial frequency
    
    
    %G
    stepR = rSF(1):rPitch:rSF(2); %Generates radial stepsize over NA range of objective

    nux = stepR' * cos(theta);    %Converts value to x and y coordinates based on angular pitch
    nuy = stepR' * sin(theta);
      stepR_c = [U(1,round(stepR./dSF(1))+(size(U,1)/2+1))' V(round(stepR./dSF(2))+(size(V,1)/2+1),1)]; %Grabs spatial frequencies within U,V grid for illumination
      ptGrid = zeros(size(U));  %Preallocates point grid for describing positions of illumination
%       if(showGrid)
%           fig = figure(120);
%           hold on;
      for k = 1:size(stepR_c)
          if(stepR_c(k,1) ~= 0 && stepR_c(k,2) ~= 0)
            tmp = circ(sqrt((U./(stepR_c(k,1)+dSF(1))).^2 + (V./(stepR_c(k,2)+dSF(2))).^2)); %Generates circle with radius equal to current radial stepsize plus an additional pixel
            tmp = tmp - circ(sqrt((U./(stepR_c(k,1)-dSF(1))).^2 + (V./(stepR_c(k,2)-dSF(2))).^2)); %Generates single pixel ring of available values for a given radius
          else
              tmp = ones(size(U));
          end
          nux(k,:) = round(nux(k,:)./dSF(1));nuy(k,:) = round(nuy(k,:)./dSF(2));  %Finds pixels values for each illumination being evaluated
          uPts = U(1,nux(k,:)+(size(U,1)/2+1));vPts = V(nuy(k,:) + (size(V,1)/2+1),1); %Snags U,V coordinates for each illumination
          ptGrid = ptGrid + tmp .* checkVals(U,uPts) .* checkVals(V,vPts);  %Adds illumination points to overall point grid mask
      end
      ptGrid(ptGrid >= 1) = 1; %Removes any values summed together indicating duplicate illuminations 

    uOut = U .* ptGrid;        %Mask U,V grid to obtain cartesian SF corresponding to each illumination
    vOut = V .* ptGrid;
    allNA = lambda .* (sqrt(uOut.^2 + vOut.^2));    %Finds NA at each illumination position to filter out larger than acceptable NA
    
    if(rNA ~= 0)
        ptGrid(allNA < rNA(1)) = 0;
    end
    ptGrid(allNA > rNA(2)) = 0; %Low-pass filter anything above maximum allowed NA in filter
    U = U'; V = -V';            %Flips coordinates to provide proper illumination angles
    NAOut(:,1) = lambda .* V(find(ptGrid == 1)); %Obtain vertical source spatial freq. positions
    NAOut(:,2) = lambda .* U(find(ptGrid == 1)); %Obtain horizontal source spatial freq. positions
    if(nargout > 1)
        thGrid = atan2d(V',U');
        rGrid = sqrt(U.^2 + V.^2);
        polOut(:,1) = rGrid(find(ptGrid == 1));
        polOut(:,2) = thGrid(find(ptGrid == 1));
    end
end %End of gType grid statement
if(showGrid)
    figure(121)
    imagesc(U(1,:),V(:,1),ptGrid);
    axis image;axis([-rSF(2) rSF(2) -rSF(2) rSF(2)]);
    disp(['Number of Points: ' num2str(sum(sum(ptGrid)))]);
end
end% End of function generateGrid
