classdef camera
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        FPS
        AD
        Idark
        ReadNoise
        bit
        Nx_cam 
        Ny_cam 
        dpixel
        x_cam
        y_cam
    end
    
    methods
        function obj = camera(model)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            switch model
                case 'PCO.edge'
            obj.FPS = 50;
            obj.AD = 4.3 ; 
            obj.Idark = 1 ; %electron/px/s
            obj.ReadNoise = 1 ; % electrons
            obj.bit = 16 ;
            obj.Nx_cam = 2048;
            obj.Ny_cam = 2048;
            obj.dpixel = 5e-6;
            
                case 'xiB-64'
            obj.FPS = 300;
            obj.AD = 9.77 ; 
            obj.Idark = 1e4 ; %electron/px/s
            obj.ReadNoise = 10 ; % electrons
            obj.bit = 10 ;
            obj.Nx_cam = 4096;
            obj.Ny_cam = 4096;    
            obj.dpixel = 5e-6;
            end
        end
        

    end
end

