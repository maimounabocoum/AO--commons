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
        
        
        function mu = average(X)
            mu = mean(X) ;
        end
        
        function PSD = PowSpecDens(obj, x )
        
            % tranpose if vector
            if size( x , 1 ) ==1
                x = x' ;
            end
        N               = size( x ,1) ;
        XDFT            = fft( x );
        XDFT            = XDFT(1:((N/2)+1),:);
        PSD           = (2/(N*obj.Fs))*abs(XDFT).^2;   
        
        end
        
        function e_t = Energy_t(obj, x )
           
            e_t = (1/(obj.Fs))*sum( abs(x).^2 ) ;
            
        end
        
        function e_psd = Energy_psd(obj, PSD_in )
           
            e_psd = sum( PSD_in ) - PSD_in(1)/2- PSD_in(end)/2 ;
            
        end        
        
    end
end

