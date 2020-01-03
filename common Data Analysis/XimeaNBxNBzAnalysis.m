%%   this is a loop analysis routine %%%

%% add relative path to current path

addpath('..\shared functions folder');
addpath('..\read and write files');

%% load data and initialize classes
clearvars

%% import LogFile
[NB,nuX0,nuZ0] = ReadLogFile('');
G = TF2D( 2^5 , 2^5 , (2^5-1)*nuX0 , (2^5-1)*nuZ0 );
ObjectFFT = zeros(2^5 , 2^5);

%% load sequence
[Filename,Foldername] = uigetfile('*.tiff','MultiSelect','on');
Nfiles = length(Filename);
if Nfiles==1
   Filename = {Filename};
end

%% load images of interest

%REF = double( importdata('Q:\datas\2019-12-10\EXP100\Ref_OnlyP40uW_nofiler.tif') );
BG      = double( importdata('D:\Data\Mai\2020-01-03\buffer\BG\BG_51211.tiff') );
MAIN    = double( importdata('D:\Data\Mai\2020-01-03\buffer\MAIN\MAIN_51719.tiff') );
REF     = double( importdata('D:\Data\Mai\2020-01-03\buffer\REF\ref_51473.tiff') );

%% define a camera
MyXimea = camera('xiB-64')    ;
MyXimea.format = 8;
MyXimea.wavelength = 780e-9;
MyXimea.IntegrationTime = 100e-6;
MyXimea = MyXimea.ResizePixels(1024,1024);

%% set references need for deconvolution
[Frame_ref,P_tot_ref]   = GetIntensity( REF , BG , MyXimea );
[Frame_main,P_tot_main] = GetIntensity( MAIN , BG , MyXimea );


%% fourier class
N   = 2^10;
F = TF2D( N , N , 1/(MyXimea.dpixel) , 1/(MyXimea.dpixel) );
%% define filter for FFT in pixel
myFilter = ImageFilter( [410 720 120 120] );
%myFilter = ImageFilter( [470 560 470 560] ); % center
BW = myFilter.getROI(N,N);

%% -------------- begin loop 
IndexRecord0 = ExtractIndex(Filename{1});
IM_rec  = zeros(N,N);

ImageCorr = zeros(1,Nfiles);

figure(3)
for loop = 1:Nfiles
    
%     if loop == 1
%     IM = double( importdata([Foldername,Filename{loop}]) );
%     [Frame0,P_frame0] = GetIntensity( IM , BG , MyXimea ) ;     
%     end
    
    IM = double( importdata([Foldername,Filename{loop}]) );
    IndexRecord = ExtractIndex(Filename{loop}) - IndexRecord0 ;
    [Frame,P_frame] = GetIntensity( IM , BG , MyXimea );
    
    
    % filter out tagged photons
    FrameFFT      = F.fourier( Frame );
    FilteredFrame = F.ifourier( FrameFFT.*BW) ;
    FilteredFrame = 2*abs(FilteredFrame).^2;%./(Frame_ref);
    FilteredFrame(abs(FilteredFrame) > 10 ) = 0 ;
    FilteredFrame(800:end,:)=0;
    % evaluation of total power
    temp = FilteredFrame;
    P_tot(loop) = sum( abs(temp(:))*MyXimea.dpixel*MyXimea.dpixel );
    

    
    % indirect reconstruction
    Nbx = NB(mod(IndexRecord,280)+1,2) ;
    Nbz = NB(mod(IndexRecord,280)+1,3) ;
    PHASE= NB(mod(IndexRecord,280)+1,4);
    
    DecalZ  =   0.3; % ??
    s = exp(2i*pi*DecalZ*Nbz);
    ObjectFFT((2^4+1)+Nbz,(2^4+1)+Nbx) = s*exp(1i*2*pi*PHASE)*P_tot(loop);
    ObjectFFT((2^4+1)-Nbz,(2^4+1)-Nbx) = conj( ObjectFFT((2^4+1)+Nbz,(2^4+1)+Nbx));

    % direct reconstruction
    %IM_rec = IM_rec + exp(1i*2*pi*PHASE)*FilteredFrame;
    
    
%imagesc( F.x*1e3 , F.z*1e3 , 100*FilteredFrame )
% imagesc( F.x*1e3 , F.z*1e3 , real(IM_rec) )
% cb = colorbar ;
% xlabel('mm')
% ylabel('mm')
% ylabel(cb,'Intensity in \mu W /cm^2')
% title(['P_{tot}',num2str(1e6*P_tot(loop)),'\mu W , (NBx,NBz)=(',num2str(Nbx),',',num2str(Nbz),')'])
% % ImageCorr(loop) = corr2(Frame,Frame0)    ;
% drawnow
end

%%

%Reconstruct = G.ifourier( ObjectFFT );
I = ifft2(fftshift(ObjectFFT));
I = I - ones(2^5,1)*I(1,:);
I = ifftshift(I,2);
figure;imagesc(I)




    