function Iout = PoissonNoise(Iin)



% Iin = (F.dx*F.dy)*Iin/Ephoton; % convert J/m^2 to photon count 

% gaussian distribution for N >100

% Iout = Iin + randn(size(Iin)).*sqrt(Iin) ;

Iout = floor( Iin + randn(size(Iin)).*sqrt(Iin) );



% Poisson noise for  N <= 100

I_poisson = find( Iin <= 0 );


Iout( I_poisson ) = poissrnd( Iin( I_poisson ) );

% Iout = Iout*Ephoton/(F.dx*F.dy) ;  % convert photon count to J


end

