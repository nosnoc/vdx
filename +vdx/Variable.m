classdef Variable < matlab.mixin.indexing.RedefinesParen &...
        matlab.mixin.indexing.RedefinesBrace &...
        handle
    properties
        %
        indices

        %
        vector

        %
        solver
    end

    methods(Access=public)
        function obj = Variable(vector, solver)
            obj.indices = {};
            obj.vector = vector;
            obj.solver = obj.solver;
        end
        function out = cat(dim,varargin)
            error('Concatenation not supported')
        end

        function varargout = size(obj,varargin)
            varargout{1} = size(obj.indices);
            %TODO(anton) needs to return correct values for varargin
        end
    end

    methods (Access=protected)
        function varargout = parenReference(obj, index_op)
            if isscalar(index_op)
                % TODO(@anton) Decide whether we squeeze, or concatenate with sorted indices.
                %              This is in my opinion purely a decision that should be made and stuck to.
                symbolics = cellfun(@(x) obj.vector.w(x), obj.indices, 'uni', false);
                out = squeeze(symbolics.(index_op));
                if isscalar(out)
                    varargout{1} = out{1};
                else
                    varargout{1} = out;
                end
            else
                if index_op(2).Type == 'Dot'
                    
                else
                    error('unsupported indexing');
                    % TODO(@anton) better error here.
                end
            end
        end

        function obj = parenAssign(obj,index_op,varargin)
            if isscalar(index_op)
                % get the cell array of args (x,x0,lbx,ubx)
                arg = varargin{1};
                indices = obj.vector.add_variable(arg{:});
                obj.indices{index_op.Indices{:}} = indices;
            else
                % TODO(anton) this can acually be used to allow mpcc.x(1,1).lb = lbx
                %             which is awesome syntactic sugar!
                error("Chained indexing not supported for NosnocVariable")
            end
        end

        function n = parenListLength(obj,index_op,ctx)
           n = 1;
        end

        function varargout = braceReference(obj, index_op)
            if isscalar(index_op)
                values = cellfun(@(x) obj.solver.results.w(x), obj.indices, 'uni', false);
                varargout{1} = squeeze(obj.symbolics.(index_op));
                return;
            else
                error('Brace indexing only accesses current values of the solution and no chained indexing')
            end
        end

        function obj = braceAssign(obj,index_op,varargin)
            if isscalar(index_op)
            end
        end

        function n = braceListLength(obj,index_op,ctx)
           n = 1;
        end

        function obj = parenDelete(obj,index_op)
            error('Deletion of symbolics is not supported through the variable view')
        end
    end

    methods (Static, Access=public)
        function obj = empty()
            obj = NosnocVariable([],[]);
        end
    end
end
