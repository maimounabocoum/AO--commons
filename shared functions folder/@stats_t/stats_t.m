classdef stats_t
    
    % this class applies to matrix input
    
    properties
        Fs
    end
    
    methods
        function obj = stats_t( Fs )
            % definition of stats_t class
            obj.Fs = Fs ;
        end
              
        function [t,mu] = average(obj,X)
            
            if size(X,1)== 1
                X = X' ;
            end
            
            mu = mean(X , 1) ;
            t = (1:size(X,1));
            
        end

        function mu = standard_dev(obj,X)
            mu = sqrt( var(X, 0 , 1) ) ;
        end
        
        function [f,PSD] = PowSpecDens(obj, x )
        
            % tranpose if vector
            if size( x , 1 ) == 1
                x = x' ;
            end
        N               = size( x ,1) ;
        f               = ( (obj.Fs)/(N) )*( 0:(N/2) );
        XDFT            = fft( x );
        XDFT            = XDFT(1:((N/2)+1),:);
        
        % homogeneous to [x^2/Hz]
        % PSD             = (2/(N*obj.Fs))*abs(XDFT).^2;  
        
        % homogeneous to [x^2/Hz]
        PSD             = (2/(obj.Fs)^2)*abs(XDFT).^2;  
        
        end
        
        function e_t = Energy_t(obj, x )
            
            if size(x,1)==1
                x = x' ;
            end
           
            e_t = (1/(obj.Fs))*sum( abs(x).^2 , 1 ) ;
            
        end
        
        function e_psd = Energy_psd(obj, PSD_in )
            
            % no frequency
            % e_psd = sum( PSD_in , 1 ) - PSD_in(1,:)/2- PSD_in(end,:)/2 ;
            % accounting for frequency integration
            
            N       = (2*size(PSD_in,1)-2);
            df      = (obj.Fs)/N ;
            %e_psd   = df*( sum( PSD_in,1) ) ;
            e_psd   = df*( sum( PSD_in , 1 ) - PSD_in(1,:)/2 - PSD_in(end,:)/2 ) ;
            
        end        
        
    end
end

