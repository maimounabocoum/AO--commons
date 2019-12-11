classdef ImageFilter
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Coord
    end
    
    methods
        function obj = ImageFilter(points)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Coord = points ;
        end
        
        function H = getROI(obj,Nx,Ny)

            H = ones(2^10,2^10);
            BOX = obj.Coord ;
            H( : ,1:Nx > BOX(2) )  = 0;
            H( :,1:Nx < BOX(1)  )  = 0;
            H( 1:Ny > BOX(4) , : ) = 0;
            H( 1:Ny < BOX(3), :  ) = 0;

        end
        
        function [] = DrawROI(obj)
            rectangle('Position',[obj.Coord(1),obj.Coord(3),obj.Coord(2)-obj.Coord(1)+0.5,obj.Coord(4)-obj.Coord(3)+0.5])
        end
    end
end

