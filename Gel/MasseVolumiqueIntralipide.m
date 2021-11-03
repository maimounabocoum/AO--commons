%% masse volumique intralipide

% mesurement 01-11-2021
M = [0.0188,0.0383,0.0584,0.0778,0.0970,0.9865,0.9841,1.9714,2.9603]; % in g
V = [20,40,60,80,100,1000,1000,2000,3000]*1e-3; % in mL

figure;plot(V,M,'-o')

mean(M./V) % masse volumique = 9.7293e-04 g/ml