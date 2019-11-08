function Ireconstruct = RetroprojectionOS_shared(I,X_m,ActiveLIST,z_out,theta, M0 , X0 , Z0 , Kx, H )
% function created by maimouna bocoum 13/09/2017

L = max(X_m) - min(X_m) ;
z_out = z_out(:)';

% check consistancy of data dimensions : 
if size(I,1)~=length(z_out)
    error('inconsistent data size')
end

% generation of mask function M + t ut
% followed by interpolation on grid 


% retroprojection : 


 [X,Z]= meshgrid(X_m,z_out);
 Ireconstruct = zeros(size(X,1),size(X,2),'like',X);

%  H = figure;
%  A = axes ;
 set(H,'WindowStyle','docked');
 
  for i= 1:2:length(theta)
       
         
        T =   (X - M0(i,1)).*sin( theta(i) ) ...
            + (Z - M0(i,2)).*cos( theta(i) ) ;
        % distance to origine
        S =   (X - M0(i,1)).*cos( theta(i) ) ...
            - (Z - M0(i,2)).*sin( theta(i) ) ;
      % common interpolation:  
        
        %Mask = double( interp1(X_m,ActiveLIST(:,i),X,'linear',0) );
        
        projContrib = interp1(z_out,I(:,i),T(:),'linear',0);
        projProfil1 = interp1(X_m,double(ActiveLIST(:,i)),S(:),'linear',0);
        projProfil2 = interp1(X_m,double(ActiveLIST(:,i+1)),S(:),'linear',0);
        
        projProfil = (projProfil1 + 1i*projProfil2) ;
        %projKx = sum(double(ActiveLIST(:,i))'.*exp(1i*(Kx(i)/L)*X_m));
        %projProfil = interp1(X_m,exp(1i*(Kx(i)/L)*X_m),S(:),'linear',0);
        
        projContrib = projContrib.*projProfil ;
       % retroprojection:  
        Ireconstruct = Ireconstruct + reshape(projContrib,length(z_out),length(X_m)); 
        %%% real time monitoring %%%   
       imagesc( X_m*1e3,z_out*1e3,real(Ireconstruct))
       colormap(parula)
       cb = colorbar ;
       title(['angle(�): ',num2str(theta(i)*180/pi),'Kx',num2str(Kx(i))])
       xlabel('x (mm)')
       ylabel('z (mm)')
       caxis( [ min(real(Ireconstruct(:))) , real(max(Ireconstruct(:))) ] )
       drawnow 

  end
  

    
    %title('Reconstruction')
    ylabel(cb,'AC tension (mV)')
    colormap(parula)
    set(findall(H,'-property','FontSize'),'FontSize',15) 

     
%     figure
%     imagesc(theta*180/pi,t*1e3,profil_1D)
%     figure;
% imagesc(theta*180/pi,z_out*1e3,I)
% hold on
% plot(theta*180/pi,...
%     1e3*(19.2e-3- M0(:,1)').*sin(theta)+1e3*(20e-3- M0(:,2)').*cos(theta),...
%     '--','linewidth',3,'color','red')

end









