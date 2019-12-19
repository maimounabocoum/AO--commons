classdef ImageFilter
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Coord
    end
    
    methods
        function obj = ImageFilter(points)
            % [centerX,centerZ,WidthX,LengthZ]
            obj.Coord = points ;
        end
        
        function H = getROI(obj,Nx,Ny)

            H = ones(Ny,Nx);
            
                  
            H( : ,1:Nx > obj.Coord(1) + obj.Coord(3)/2 )  = 0;
            H( :,1:Nx < obj.Coord(1) - obj.Coord(3)/2 )  = 0;
            H( 1:Ny > obj.Coord(2) + obj.Coord(4)/2 , : ) = 0;
            H( 1:Ny < obj.Coord(2) - obj.Coord(4)/2, :  ) = 0;

        end
        
        function [] = DrawROI(obj)
            BOX = [obj.Coord(1) - obj.Coord(3)/2,...
                   obj.Coord(1) + obj.Coord(3)/2,...
                   obj.Coord(2) - obj.Coord(4)/2,...
                   obj.Coord(2) + obj.Coord(4)/2 ];
            rectangle('Position',[BOX(1)-0.5,BOX(3)-0.5,BOX(2)-BOX(1)+0.5,BOX(4)- BOX(3)+0.5])
        end
    end
end

