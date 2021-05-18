% Copyright (C) 2020 Frank E. Curtis
%
% All Rights Reserved.
%
% Authors: Frank E. Curtis

% MeritParameterComputationModelReduction: addOptions
function addOptions(options,reporter)

% Add options
options.addBoolOption(reporter,'MCMR_linear_model',true);
options.addBoolOption(reporter,'MCMR_quadratic_model',false);
options.addDoubleOption(reporter,'MCMR_curvature_threshold',5e-09,0,inf);
options.addDoubleOption(reporter,'MCMR_model_reduction_factor',1e-01,0,1);
options.addDoubleOption(reporter,'MCMR_parameter_reduction_factor',1e-02,0,1);
options.addDoubleOption(reporter,'MCMR_normal_progress_factor',1-1e-04,0,1);

end % addOptions