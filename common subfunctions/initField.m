function [E0_tag,E0_untag,Eref] = initField(P0,Pref,ModeWidth,eta,F)


[X,Z]   = meshgrid(F.x,F.z);
[FX,FZ] = meshgrid(F.fx,F.fz);


%% ====================== Phase randomization induced by scattering phantom


% IS0         = exp(-(FX.^2+FZ.^2)/(ModeWidth)^2); % initial Main laser beam 
% % generation of tagged photons
% PHASE       = 2*pi*rand(size(X)) ;
% IS0_tag     = IS0.*exp(1i*PHASE);
% % generation of untagged photons
% PHASE       = 2*pi*rand(size(X)) ;
% IS0_untag   = IS0.*exp(1i*PHASE);
% 
% % ============================= field out of Phantom
% E0_tag      = F.ifourier(IS0_tag);
% E0_untag    = F.ifourier(IS0_untag);
% 
% %% old code
% IS0ref      = exp(-(FX.^2+FZ.^2)/(0.00001*ModeWidth)^2);                 % initial Main laser beam 
% % speckle generation in reference beam
% PHASE       = 2*pi*rand(size(X)) ;
% Iref        = IS0ref.*exp(1i*PHASE);
% Eref        = F.ifourier(Iref);


  E0_tag        = exp( -(X.^2+Z.^2)/(1e-3)^2 ) ;  	% initial Ref laser beam 
  E0_untag      = E0_tag;  	% initial Ref laser beam 
  Eref          = exp(-(X.^2+Z.^2)/(1.5e-3)^2).*exp(1i*(15e4)*Z);  	% initial Ref laser beam 
 % Eref       = exp(1i*15e4*X);  	% initial Ref laser beam 

 
 
 
% % convert corresponding intensitZ to m^{-2}

E0_tag     = E0_tag/sqrt( trapz(F.x, trapz(F.z,abs(E0_tag).^2 )) );
E0_untag   = E0_untag/sqrt( trapz(F.x,trapz(F.z,abs(E0_untag ).^2)));
Eref       = Eref/sqrt( trapz(F.x,trapz(F.z,abs(Eref).^2)) );


% % convert intensitZ to W / m^{-2}
 E0_tag    = sqrt(eta*P0).*E0_tag  ;
 E0_untag  = sqrt((1-eta)*P0).*E0_untag ;
 Eref      = sqrt(Pref).*Eref ;




















end

