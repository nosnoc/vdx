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
            obj.w = vdx.Vector(obj, -inf, inf, 0);
            obj.p = vdx.Vector(obj, -inf, inf, 0);
            obj.g = vdx.Vector(obj, 0, 0, 0);
            obj.f = 0;
        end
    end
    
    methods (Access=protected)
    end
end
