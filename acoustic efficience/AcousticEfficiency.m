%% acoustic efficiency

x = 0:0.1:10;



P = (0.1:0.1:10)*1e6;
lambda = 780e-9;
dndp = 1.32e-10;
L = 1e-3;

Phi = 2*pi*dndp*L*P/lambda;

J1 = besselj(1,Phi);

nopt = 1.28*lambda/(2*pi*L);

L = (1:0.1:10)*1e-3;
Gamma = 1540/(6e6);

G = 2*pi*lambda*L/(1.33*Gamma^2);
figure(1)
plot(L,G)

figure(2);plot(Phi,J1)
xlabel('Pression(MPa)')
ylabel('J_1(\phi)')