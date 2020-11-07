%%% choose lens focal

%% relation of conjugated plane :

% w0' = (lambda*f)/(pi*w0)
clearvars


lambda = 780e-9;
f = (950:1:1050)*1e-3 ;
d0 = 1.5e-3;

d0_p = 2*(lambda*f)/(pi*(d0/2));
Zr = pi*(d0_p/2).^2/lambda ;

figure(3);
subplot(121)
plot(f*1e3,Zr*1e3)
xlabel('focal(mm)')
ylabel('raigley length(mm)')
subplot(122)
plot(f*1e3,d0_p*1e3)
xlabel('focal(mm)')
ylabel('diameter(mm)')

% Zr = pi*(d0/2)^2/lambda
%% diverge due to iris
lambda = 780e-9;
w0 =  (1:0.1:100)*1e-6; % iris size
theta = atan(lambda/pi./w0); % half angle


figure(1);
plot(w0*1e6,(180/pi)*theta)
xlabel('waist(\mu m)')
ylabel('\theta (°)')