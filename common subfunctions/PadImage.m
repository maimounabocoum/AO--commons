function Iout = PadImage(Iin,Nx,Ny)

N_xin = size(Iin,2);
N_yin = size(Iin,1);
Nc_xin = floor(N_xin/2) + 1 ; % image center X
Nc_yin = floor(N_yin/2) + 1 ; % image center Y


Iout = zeros(Ny,Nx);


[IX, IY] = meshgrid(( -(floor(N_xin/2)):(floor(N_xin/2)-1) ) + Nc_xin,...
                   ( -(floor(N_yin/2)):(floor(N_yin/2)-1) ) + Nc_yin ) ;
 
 
               IX = repmat(IX,ceil(Ny/N_yin),ceil(Nx/N_xin));
               IX = IX(1:Ny,1:Nx);
               IY = repmat(IY,ceil(Ny/N_yin),ceil(Nx/N_xin));
               IY = IY(1:Ny,1:Nx);
linearInd = sub2ind( [N_yin,N_xin] , IY(:) , IX(:) ) ;
% Ifill_x = ( -(floor(N_xin/2)):(floor(N_xin/2)-1) ) + Nc_xin ;


Iout(:) = Iin(linearInd );






end

