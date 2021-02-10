%%   this is a loop analysis routine %%%

%% add relative path to current path

addpath('..\shared functions folder');
addpath('..\read and write files');

%% load data and initialize classes
clearvars

%% import LogFile
[NB,nuX0,nuZ0] = ReadLogFile('E:\2021-01-04') ;

%% load sequence
[Filename,Foldername] = uigetfile('E:\2021-01-04\*.tiff','MultiSelect','on');
Nfiles = length(Filename);
if Nfiles==1
   Filename = {Filename};
end

%% load images of interest

%REF = double( importdata('Q:\datas\2019-12-10\EXP100\Ref_OnlyP40uW_nofiler.tif') );
IM = double( importdata([Foldername,Filename{1}]) );
BG = double( importdata('E:\2021-01-04\refs\BG.tiff') );
MAIN = double( importdata('E:\2021-01-04\refs\MAIN.tiff') );
REF = double( importdata('E:\2021-01-04\refs\REF.tiff') );

%% define a camera
MyXimea = camera('xiB-64')    ;
MyXimea.format = 8;
MyXimea.wavelength = 780e-9;
MyXimea.IntegrationTime = 100e-6;
MyXimea = MyXimea.ResizePixels(1024,1024);

%% set references need for deconvolution
[Frame,P_tot_main] = GetIntensity( IM , BG , MyXimea );

%% fourier class
N   = 2^10;
F = TF2D( N , N , 1/(MyXimea.dpixel) , 1/(MyXimea.dpixel) );

%% define filter for FFT in pixel
myFilter = ImageFilter( [613 513 20 20] );
myFilter.DrawROI ;
BW = myFilter.getROI(N,N);



%% -------------- begin loop 
IndexRecord0 = ExtractIndex(Filename{1},'labview');
MU  = zeros(N,N);
MU2 = zeros(N,N);
P_tot = zeros(1,Nfiles);
mu  =  0;
mu2 = 0;
ImageCorr = zeros(1,Nfiles);
[X,Y] = meshgrid(MyXimea.x_cam*MyXimea.dpixel-mean(MyXimea.x_cam*MyXimea.dpixel),MyXimea.y_cam*MyXimea.dpixel-mean(MyXimea.y_cam*MyXimea.dpixel));
figure(3)
hold off
%PhaseNoise = rand(1,Nfiles);
%AmpNoise = rand(1,Nfiles);
for loop = 1:Nfiles
    
%     if loop == 1
%     IM = double( importdata([Foldername,Filename{loop}]) );
%     [Frame0,P_frame0] = GetIntensity( IM , BG , MyXimea ) ;     
%     end
    
    IM = double( importdata([Foldername,Filename{loop}]) );
    IndexRecord = ExtractIndex(Filename{loop},'labview') - IndexRecord0 ;
    [Frame,P_frame] = GetIntensity( IM , BG , MyXimea );
    
 %   Frame = abs(1 + exp(-(X.^2+Y.^2)/(5e-3)^2+ 0.3*AmpNoise(loop)).*exp(1i*(1.52e5*X + PhaseNoise(loop) )) ).^2 + rand(1024,1024)*0.1;
    
    % filtering of tagged photons
    
    FrameFFT      = F.fourier( Frame )           ; % image FTT

    %figure ; imagesc(log(abs(FrameFFT)))        ;
    FilteredFrame = F.ifourier( FrameFFT.*BW )   ; % fourier filtering
    %FilteredFrame = imag(FilteredFrame)     ; % intensity of fourier window
 
    % evaluation of total power
    temp = Frame ;
    P_tot(loop) = sum( abs(temp(:))*MyXimea.dpixel*MyXimea.dpixel );


    mu  = mu  + P_tot(loop)/Nfiles;
    mu2 = mu2 + P_tot(loop)^2/Nfiles;
    MU  = MU  + temp/Nfiles;    % iterative average
    MU2 = MU2 + temp.^2/Nfiles; % iterative variance
    
    intgram = Frame(500:512,200:800);
    r = real(FilteredFrame(500:512,200:800));
    Im = imag(FilteredFrame(500:512,200:800));
    a = angle(FilteredFrame(500:512,200:800));
    muav = MU_av(500:512,200:800);
%     m = abs(FilteredFrame(400:512,500:512));
%imagesc( log(abs(FrameFFT)) ); myFilter.DrawROI ;

% %scatter(unwrap(a(:))/(2*pi),intgram(:).^2)

scatter3(r(:), Im(:) ,intgram(:)- muav(:),1,intgram(:)- muav(:))
view([0 90])
 hold on
axis equal
% title('IM - <IM>')
%   subplot(2,2,3)
%   hold on
%   scatter( a(:) , r(:))
%  xlabel('phase(Rad)')
%  title('Real part inteferogram')
%  subplot(224)
%  hold on
%  scatter( a(:) , Im(:) )
%  xlabel('phase(Rad)')
%  title('Imag part inteferogram')


% subplot(122)
% scatter(a(:),m(:))
% hold on
% imagesc( F.x*1e3 , F.z*1e3 , temp )
% cb = colorbar ;
% xlabel('mm')
% ylabel('mm')
% ylabel(cb,'Intensity in \mu W /cm^2')
% title(['P_{tot}',num2str(1e6*P_tot(loop))])
% ImageCorr(loop) = corr2(Frame,Frame0)    ;

drawnow
%saveas(gcf,['Q:\datas\simulated datas\2021-01-24\IM_',num2str(loop)],'png')
end


%%





%% obtain STD in 2-steps
%1-step: get average
% MU_  = zeros(N,N);
% MU2_ = zeros(N,N);
% for loop = 1:Nfiles
% 
%     IM = double( importdata([Foldername,Filename{loop}]) );
%     IndexRecord = ExtractIndex(Filename{loop},'labview') - IndexRecord0 ;
%     [Frame,P_frame] = GetIntensity( IM , BG , MyXimea );
%     MU_  = MU_  + Frame;        % iterative average    
% end
% MU_ = MU_/Nfiles;
% % 2-step: get std
% for loop = 1:Nfiles
% 
%     IM = double( importdata([Foldername,Filename{loop}]) );
%     IndexRecord = ExtractIndex(Filename{loop},'labview') - IndexRecord0 ;
%     [Frame,P_frame] = GetIntensity( IM , BG , MyXimea );
%     MU2_  = MU2_ + ( Frame - MU_ ).^2;    % iterative average
%     
% end
% MU2_ = MU2_/Nfiles;
% 
% figure(3)
% imagesc(F.x*1e3 , F.z*1e3 , sqrt(MU2_))
% cb = colorbar ;
% xlabel('mm')
% ylabel('mm')

%% plot result
sigma = sqrt( mu2 - mu^2 );
sigma_sn = MyXimea.GetShotNoise(mu);
SIGMA =  sqrt(MU2 - MU.^2) ;
MU_av = MU;
MU_fourier = F.fourier( MU );
FilteredMU = F.ifourier( MU_fourier.*BW )   ; % fourier filtering
    figure(1)
    imagesc( 1e3*MyXimea.x_cam*MyXimea.dpixel , 1e3*MyXimea.x_cam*MyXimea.dpixel , 100*MU)
    cb = colorbar ;
    xlabel('mm')
    ylabel('mm')
    ylabel(cb,'Intensity in \mu W /cm^2')
    title(['\mu =',num2str(1e6*mu),'\mu W ; \sigma_{sn} = ',num2str(num2str(1e6*sigma_sn)),'\mu W'])
    drawnow
    figure(2)
    imagesc( 1e3*MyXimea.x_cam*MyXimea.dpixel , 1e3*MyXimea.x_cam*MyXimea.dpixel , 100*SIGMA)
    cb = colorbar ;
    xlabel('mm')
    ylabel('mm')
    ylabel(cb,'Intensity in \mu W /cm^2')
    title(['\sigma = ',num2str(num2str(1e6*sigma)),'\mu W'])
    drawnow
    
  

%     %% plot results
%     Ptot = [0.0056421, 0.025105, 0.15004,0.71961,3.5872,13.5005];
%     Sigma = [9.3106e-5,0.00070026,0.0023577,0.017529,0.091051,0.23117];
%     Sigma_sn = [3.7919e-6,7.9987e-6,1.9554e-5,4.2824e-5,9.5613e-5,0.00018549];
%     
%%  

    figure(3);
    hold on
    plot(1e3*MyXimea.x_cam*MyXimea.dpixel , 100*MU(512,:));
    hold on
    plot(1e3*MyXimea.x_cam*MyXimea.dpixel , 100*SIGMA(512,:) );
    xlim([1.5 2.5])
    xlabel('mm')
    ylabel('\mu W /cm^2')
    title('fringes comparison')
    legend('\mu','\sigma')

    
    %%
    
    
    
    
    
    
    
    
    
    