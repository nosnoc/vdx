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
    end
end
