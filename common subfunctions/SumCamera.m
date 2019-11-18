function IcamOut = SumCamera(x_cam,y_cam,x,y,I,np)

dx = x(2) - x(1) ;
dy = y(2) - y(1) ;

% average power over np consecutive point in (F.x,F.y) grid
M = movmean(I,[np,np]);

Ix = interp1(x,1:length(x),x_cam,'nearest',0);
Iy = interp1(y,1:length(y),y_cam,'nearest',0);
[IX,IY] = meshgrid(Ix,Iy);

% row    : Iy
% column : Ix
linearInd = sub2ind(size(I), IY(:), IX(:)) ;
 
IcamOut = M;
IcamOut = reshape(IcamOut(linearInd),[length(y_cam),length(x_cam)]);
IcamOut = IcamOut.*(dx*dy*np^2);


% figure(1);
% subplot(121) ;imagesc(x*1e3,y*1e3,I) ; colorbar
% hold on
% subplot(121) ;plot(X(linearInd)*1e3,Y(linearInd)*1e3,'x','color','red')
% % axis([-0.05 0.05 -0.05 0.05])
% hold off
% figure(1);
% subplot(122) ; imagesc(IcamOut); colorbar






















end

