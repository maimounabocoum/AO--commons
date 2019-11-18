classdef TF2D
    %TF2D Summarz of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        Nx       % number of point for fourier transform
        Nz
        xRange  
        x       % position in m
        fx       %frequencz varaiable
        kx       %frequencz variable in rad
        dx
        dfx
        zRange  
        z       % position in m
        fz       %frequencz varaiable
        kz       %frequencz variable in rad
        dz
        dfz
    end
    
    methods
        
        function obj = TF2D(Nx,Nz,Fmax_x,Fmax_z)
            % Fmax corresponds to the mas sampling of image
            
            obj.Nx = Nx;
            obj.Nz = Nz;
            obj.dx = 1/Fmax_x; % in m
            obj.dz = 1/Fmax_z; % in m
            obj.x = (-Nx/2:1:Nx/2-1)*obj.dx;
            obj.z = (-Nz/2:1:Nz/2-1)*obj.dz;
            obj.xRange = (Nx-1)*obj.dx;
            obj.zRange = (Nz-1)*obj.dz;
            obj.dfx = 1/obj.xRange;
            obj.dfz = 1/obj.zRange;
            obj.fx = (-Nx/2:1:Nx/2-1)*obj.dfx;
            obj.fz = (-Nz/2:1:Nz/2-1)*obj.dfz;
            obj.kx = 2*pi*obj.fx;
            obj.kz = 2*pi*obj.fz;
            
        end
        
        function Ekxkz = fourier(obj, Exz)
            %fftshift(Et);    %real(F) sera toujours positif pour phi=0
            Ekxkz=fft2(ifftshift(Exz))*(obj.xRange/obj.Nx)*(obj.zRange/obj.Nz) ;
            Ekxkz=fftshift(Ekxkz);
        end
        
        function Exkz = fourierz(obj, Exz)
            %fftshift(Et);    %real(F) sera toujours positif pour phi=0
            Exkz = fft(ifftshift(Exz,1),obj.Nz,1)*(obj.zRange/obj.Nz) ;
            Exkz = fftshift(Exkz,1);
        end
         
        function Exz = ifourierz(obj, Exkz)
            %fftshift(Et);    %real(F) sera toujours positif pour phi=0
            Exz = ifft(fftshift(Exkz,1),obj.Nz,1)*(obj.Nz/obj.zRange) ;
            Exz = ifftshift(Exz,1);
        end
        
        function Ekxz = fourierx(obj, Exz)
            %fftshift(Et);    %real(F) sera toujours positif pour phi=0
            Ekxz = fft(ifftshift(Exz,2),obj.Nx,2)*(obj.xRange/obj.Nx) ;
            Ekxz = fftshift(Ekxz,2);
        end
        
        function Ekxz = ifourierx(obj, Exz)
            %fftshift(Et);    %real(F) sera toujours positif pour phi=0
            Ekxz = ifft(fftshift(Exz,2),obj.Nx,2)*(obj.Nx/obj.xRange)  ;
            Ekxz = ifftshift(Ekxz,2);
        end
        
        function Exz = ifourier(obj, Ekxkz)
            %fftshift(Et);    %real(F) sera toujours positif pour phi=0
            Exz=ifft2(ifftshift(Ekxkz))*(obj.Nz/obj.zRange)*(obj.Nx/obj.xRange)  ;
            Exz=fftshift(Exz);
        end
        
    end
    
end

