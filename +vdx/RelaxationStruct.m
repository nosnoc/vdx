classdef RelaxationStruct
    properties (SetAccess=private)
        mode vdx.RelaxationMode
        slack_name
        weight_name
    end

    methods
        function obj = RelaxationStruct(mode, slack_name, weight_name);
            obj.mode = mode;
            obj.slack_name = slack_name;
            obj.weight_name = weight_name;
        end

        function is_relaxed = is_relaxed(obj)
            is_relaxed = ~(obj.mode == vdx.RelaxationMode.NONE);
        end
    end
end
