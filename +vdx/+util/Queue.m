classdef Queue < handle
    properties
        data
    end

    methods
        function obj = Queue()
            obj.data = cell(0,1);
        end

        function put(obj, elem)
            obj.data{end + 1} = elem;
        end

        function elem = get(obj)
            if ~obj.isempty()
                elem = obj.data{end};
                obj.data(end) = [];
            else
                errID = 'Queue:empty';
                msgtext = 'Cant get from empty queue.';

                err = MException(errID,msgtext);
                throw(err);
            end
        end

        function elems = get_n(obj, n)
            if n < obj.length
                elems = obj.data((end-n):end);
                obj.data((end-n):end) = [];
            else
                errID = 'Queue:empty';
                msgtext = ['Cant get ' num2char(n) ' elements from queue.'];

                err = MException(errID,msgtext);
                throw(err);
            end
        end

        function empty = isempty(obj)
            empty = obj.length == 0;
        end

        function n = length(obj)
            n = length(obj.data);
        end
    end
end
