function [angle, M0 , X0 , Z0] = EvalDelayLaw_shared( X_m , DelayLAWS  , ActiveLIST , c)
%creation of function 11/09/2017 for delay law retreival using sequence
%parameters. 

screenPlot = 0 ;

% DelayLAWS : each column represents the delay law in s for all probe
% elements. Note: non active elements are set to zero
 for i = 1:size(DelayLAWS,2)
      ct(i,:) =    DelayLAWS(:,i)*c;
 end
 
 % ct = M0 M(t)
Nangle = size(ct , 1) ;

% initialization:
angle = zeros(1,Nangle);
M0 = zeros(Nangle,2);

% 0 is defined by the (0,0) on probe linear plane

if screenPlot == 1
  Hf = figure;
  set(Hf,'WindowStyle','docked');
  subplot(211)
    cc = jet(Nangle);
    colormap(jet);
end
    
    for i = 1:Nangle     
       % find index of minimum active element :
       Nmin = find(ActiveLIST(:,i) == 1, 1 );
       Nmax = find(ActiveLIST(:,i) == 1, 1, 'last' );
        

       angle(i) = asin( (ct(i,Nmax)-ct(i,Nmin))/(X_m(Nmax) - X_m(Nmin)) );
       u = [sin(angle(i)),cos(angle(i))];
       % definition ut vector :
              % affichage des coordonn�e initiales :
       % (X_m,0) - ct*ut
       X0 = X_m - u(1)*DelayLAWS(:,i)'*c ;
       Z0 = 0   - u(2)*DelayLAWS(:,i)'*c;
       M0(i,1) = 0 - u(1)*DelayLAWS(1,i)'*c; 
       M0(i,2) = 0   - u(2)*DelayLAWS(1,i)'*c;
       if screenPlot == 1
       plot( X0*1e3 , Z0*1e3 ,'linewidth',3,'color',cc(i,:)) 
       hold on
       end
       
    end

if screenPlot == 1    
    cb = colorbar ;
    ylabel(cb,'angular index (a.u.)')
    xlabel('X (mm)')
    ylabel('Z(mm)')
    set(gca,'YDir','reverse')
    
subplot(212)
plot(180/pi*angle,'-o','linewidth',3)
xlabel('shoot index')
xlabel('angle (�)')
    
end







end

