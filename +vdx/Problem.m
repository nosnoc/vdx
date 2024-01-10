classdef Problem < handle 
    properties (Access=public)
        % Primal variabiles
        w
        % Constraints
        g
        % Parameters
        p
        % Objective
        f
    end
    properties (Access=private)
        % Variable fun
        variables struct
        constraints struct
    end

    methods (Access=public)
        function obj = Problem()
            obj.w = vdx.Vector(obj);
            obj.p = vdx.Vector(obj);
            obj.g = vdx.Vector(obj);
            obj.f = 0;
        end
    end
    
    methods (Access=protected)
    end
end
