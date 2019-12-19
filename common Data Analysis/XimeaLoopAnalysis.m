%%   this is a loop analysis routine %%%

%% add relative path to current path

addpath('..\shared functions folder');

%% load data and initialize classes
clearvars

%% load sequence
[Filename,Foldername] = uigetfile('*.tiff','MultiSelect','on');
Nfiles = length(Filename);
if Nfiles==1
   Filename = {Filename};
end

%% load images of interest

%REF = double( importdata('Q:\datas\2019-12-10\EXP100\Ref_OnlyP40uW_nofiler.tif') );
BG = double( importdata('Q:\datas\2019-12-10\EXP100\BG_EXP100us.tif') );
MAIN = double( importdata('Q:\datas\2019-12-10\EXP100\Main_OnlyP9uW_EXP100us.tif') );
REF = double( importdata('Q:\datas\2019-12-10\EXP100\Ref_OnlyP40uW_EXP100us.tif') );

%% define a camera
MyXimea = camera('xiB-64')    ;
MyXimea.format = 8;
MyXimea.wavelength = 780e-9;
MyXimea.IntegrationTime = 100e-6;

%% set references need for deconvolution
[Frame_ref,P_tot_ref] = GetIntensity( REF , BG , MyXimea );
[Frame_main,P_tot_main] = GetIntensity( MAIN , BG , MyXimea );


%% fourier class
N   = 2^10;
F = TF2D( N , N , 1/(MyXimea.dpixel) , 1/(MyXimea.dpixel) );
%% define filter for FFT in pixel
myFilter = ImageFilter( [653 531 95 95]  );
%myFilter = ImageFilter( [470 560 470 560] ); % center
BW = myFilter.getROI(N,N);

%% -------------- begin loop 
MU  = zeros(N,N);
MU2 = zeros(N,N);
P_tot = zeros(1,Nfiles);
mu  =  0;
mu2 = 0;
ImageCorr = zeros(1,Nfiles);

for loop = 1:Nfiles
    
%     if loop == 1
%     IM = double( importdata([Foldername,Filename{loop}]) );
%     [Frame0,P_frame0] = GetIntensity( IM , BG , MyXimea ) ;     
%     end
    
    IM = double( importdata([Foldername,Filename{loop}]) );
    [Frame,P_frame] = GetIntensity( IM , BG , MyXimea );
    
    
    % filter out tagged photons
    FrameFFT      = F.fourier( Frame );
    FilteredFrame = F.ifourier( FrameFFT.*BW) ;
    FilteredFrame = 2*abs(FilteredFrame).^2./(Frame_ref);
    
    % evaluation of total power
    temp = Frame;
    P_tot(loop) = sum( abs(temp(:))*MyXimea.dpixel*MyXimea.dpixel );

    mu  = mu  + P_tot(loop)/Nfiles;
    mu2 = mu2 + P_tot(loop)^2/Nfiles;
    MU  = MU  + temp/Nfiles; % iterative average
    MU2 = MU2 + temp.^2/Nfiles; % iterative variance
    
    
    % ImageCorr(loop) = corr2(Frame,Frame0)    ;
end

%% plot result
sigma = sqrt( mu2 - mu^2 );
sigma_sn = MyXimea.GetShotNoise(mu);
SIGMA = sqrt( MU2 - MU.^2 );
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
    
%     figure;
%     hold on
%     plot(Ptot,100*Sigma./Ptot);
  %  hold on
  %  semilogx(Ptot- Ptot(1),(Ptot- Ptot(1))./Sigma);
    
    
    