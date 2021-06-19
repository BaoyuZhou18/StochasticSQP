% Copyright (C) 2020 Frank E. Curtis
%
% All Rights Reserved.
%
% Authors: Frank E. Curtis

% StepsizeComputationAdaptive: computeStepsize
function computeStepsize(S,options,quantities,reporter,strategies)

% Set scaling
stepsize_scaling = quantities.stepsizeScaling;
if quantities.stepsizeDiminishing == true
    stepsize_scaling = stepsize_scaling / quantities.iterationCounter;
end

% Check whether primal space update being zero or not
% make 1e-16 as an option
if norm(quantities.directionPrimal,inf) < S.direction_norm_tolerance_
    
    alpha = 1;
    
else
    
    % Compute preliminary values
    denominator = (quantities.meritParameter * quantities.objectiveLipschitzConstants + quantities.constraintLipschitzConstants) * norm(quantities.directionPrimal)^2;
    alpha_1 = min(1 , 2 * (1 - S.sufficient_decrease_) * stepsize_scaling * quantities.modelReduction / denominator);

    % Set projection bounds
    lower_bound = 2 * (1 - S.sufficient_decrease_) * stepsize_scaling * quantities.ratioParameter * quantities.meritParameter / (quantities.meritParameter * quantities.objectiveLipschitzConstants + quantities.constraintLipschitzConstants);
    upper_bound = lower_bound + S.projection_width_ * stepsize_scaling^2;
    
    if S.forward_lengthening_
        if alpha_1 < 1
            while alpha_1 < 1
                alpha_ext = S.lengthening_ratio_ * alpha_1;
                
                Ufunc = alpha_ext * (S.sufficient_decrease_ - 1) * stepsize_scaling * quantities.modelReduction ...
                    + norm(quantities.currentIterate.constraintFunctionEqualities(quantities) + alpha_ext * quantities.currentIterate.constraintJacobianEqualities(quantities) * quantities.directionPrimal) ...
                    - quantities.currentIterate.constraintNorm2(quantities) ...
                    + alpha_ext * (quantities.currentIterate.constraintNorm2(quantities) - quantities.residualDualNorm2) ...
                    + 0.5 * alpha_ext^2 * denominator;
                
                if Ufunc > 0
                    break;
                else
                    alpha_1 = alpha_ext;
                end
            end
            % alpha_1 = min(alpha_1,1);
        end
    end
    % Project values
    alpha = max(lower_bound,min(alpha_1,upper_bound));
end

% Set stepsize
quantities.setStepsize(alpha);

% Create trial iterate
trial_iterate = Point(quantities.currentIterate,quantities.currentIterate.primalPoint + quantities.stepsize * quantities.directionPrimal);

% Set trial iterate
quantities.setTrialIterate(trial_iterate);

end % computeStepsize