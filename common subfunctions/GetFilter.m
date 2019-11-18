function Filt = GetFilter(nx0,ny0,Nx,Ny, D_in)

% function created by maimouna bocoum 18-11-2019

% Ix = 1:size(D_in,2);
% Iy = 1:size(D_in,1);

Filt = 0*D_in;

Filt( ny0:(ny0+Ny) , nx0:(nx0+Nx) ) = 1 ;





end
