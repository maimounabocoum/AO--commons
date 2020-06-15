
% calculation of Fresnel number:
lambda = 800e-9;


%% variation of d distance to object

a = 1e-3;%m
d = 0.5:0.1:10;%m

Np = a^2./(lambda*d);

figure(1)
plot(d,Np)
xlabel('d(m)')
ylabel('N_p')
title('Fresnel Number')

%% variation of radius a at focus of a lens

a = 1e-6*(0.1:0.1:500);%m
f = 10e-2;%m

Np = a.^2/(lambda*f);

figure(1)
plot(a*1e6,Np)
xlabel('a(\mu m)')
ylabel('N_p')
title('Fresnel Number')

%% caracteristic lenth in filter
lambda=800e-9;
f = 10e-2;
D = 1e-6*(10:0.02:100);
L = (1.22*lambda*f)./D;

figure(1)
plot(D*1e6,L*1e3)
xlabel('D(\mu m)')
ylabel('Length (mm)')
title('Caracteristic length for filter convolution')

%%


