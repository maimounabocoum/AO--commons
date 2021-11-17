%% plot curves as function of wavelength
clearvars;

lambda = 800;
concAgar = 2.79/182 ;
concInkSolution = 0; % en pourcentage massique
concIntra = 1.0149/182 ;

for nplot = 1:length(concIntra)
% inpout masse volumique

[Youngmodulus,mu_s,mu_sprime,mu_a]=gelind(lambda,concAgar,concIntra(nplot),concInkSolution)

end



%% concentration at 800nm for increasing mu_s and mu_s_prim
Youngmodulus        = 0;        % pure dilution in water
lambda              = 800 ;     % center pulse wavelength
mu_a                = 0;        % ink absorption contribution
attenuationBase10   = 0:10;     % ballistic absoption
L                   = 1 ;       % cell length in cm
mu_sprime           = attenuationBase10*log(10)/L ; %cm-1
mu_sprime           = mu_sprime ; %m-1 to cm-1
massGel             = 3;

[massAgar,massIntralipid,massInkSolution,massWater] = geldir(lambda,Youngmodulus,mu_sprime,mu_a,massGel);

Values = [massIntralipid'/0.9729,massWater'/1]*1e3; % in mul

figure;
plot(attenuationBase10, massIntralipid/massGel)
title('volume of intralipid in 3ml of water')
xlabel('attenuation')
ylabel('mass dillution coefficient')
%%