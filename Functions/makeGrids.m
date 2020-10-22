function [grid] = makeGrids(par)
%%
%makeGrids Function
%
%Purpose: This is a small function for generating spatial frequency grids
%using the meshGrid built-in matlab function.
%       
%
%Inputs: par - Structure with following parameters
%              nPix = 2x1 vector containing number of pixels for simulating
%                     grid in y and x
%              pixL = length of pixel for generating SF sampling rate
%              Mag = magnification if pixel length needs to be
%                    demagnified/magnified first
%Outputs: grid - structure containing following parameters
%              x = Nx1 vector containing x values of real image space 
%                   (assumes square image)
%                  
%              y = Nx1 vector containing y values of real image space
%                  (assumes square image)
%              u = Nx1 vector containing horizontal values of Fourier image
%                  space. (assumes square image)
%              v = Nx1 vector containing vertical values of Fourier image
%                  space. (assumes square image)
%              U = NxN matrix of horizontal Fourier image space values.
%              V = NxN matrix of vertical Fourier image space values.
%              X = NxN matrix of horizontal real image space values.
%              Y = NxN matrix of vertical real image space values.
%-------------------------------------------------------------------------%

idx = [-par.Np(1)/2:par.Np(1)/2-1];
grid.x = par.pL .* idx;
grid.y = par.pL .* idx;
grid.u = par.Mag/(par.Np(1).*par.pL) .* idx;
grid.v = par.Mag/(par.Np(2).*par.pL) .* idx;

[grid.U,grid.V] = meshgrid(grid.u,grid.v);
[grid.X,grid.Y] = meshgrid(grid.x,grid.y);

end %End of Function