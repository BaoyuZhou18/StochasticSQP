% Copyright (C) 2020 Frank E. Curtis
%
% All Rights Reserved.
%
% Authors: Frank E. Curtis

% DirectionComputationIEQP: printIterationValues
function printIterationValues(S,quantities,reporter)

% Get multipliers
[yE,yI] = quantities.currentIterate.multipliers;

% Print information
reporter.printf(Enumerations.R_SOLVER,Enumerations.R_PER_ITERATION,...
  ' %+e %+e %+e    %d    %+e %+e ',...
  norm(quantities.directionPrimal,inf),...
  norm([yE;yI],inf),...
  quantities.currentIterate.stationarityMeasure(quantities),...
  quantities.terminationTestNumber,...
  norm(quantities.residualPrimal,inf),...
  norm(quantities.residualDual,inf));

end % printIterationValues