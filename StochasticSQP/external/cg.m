function x = cg(A,b,x,e)

% function x = cg(A,b,x,e)
%
% Author      : Frank E. Curtis
% Description : Conjugate Gradient method (preliminary version).
% Input       : A ~ symmetric matrix
%               b ~ right-hand side vector
%               x ~ initial iterate
%               e ~ solution tolerance
% Output      : x ~ final iterate
% Revised by  : (your name here)

% % Store output strings
% out_line = '==============================+=====================================';
% out_data = '  k       ||x||       ||r||   |    ||p||        alpha        beta';
% out_null =                                 '----------  -----------  -----------';
% 
% % Print output header
% fprintf('%s\n%s\n%s\n',out_line,out_data,out_line);

% Evaluate initial iterate norm
norms.x = norm(x);

% Evaluate initial residual
r = A*x-b;

% Evaluate residual norm
norms.r = norm(r);

% Store initial residual norm
norms.r0 = norms.r;

% Evaluate initial direction
p = -r;

% Evaluate direction norm
norms.p = norm(p);

% Initialize iteration counter
k = 0;

% Main CG loop
while 1
  
%   % Print iterate information
%   fprintf('%5d  %.4e  %.4e | ',k,norms.x,norms.r);
  
  % Check for termination
  if norms.r <= e*max(norms.r0,1), break; end;
  
  % Evaluate matrix-vector product
  Ap = A*p;
  
  % Evaluate vector-vector product
  pAp = p'*Ap;
  
  % Evaluate steplength
  alpha = -(r'*p)/pAp;
  
  % Update iterate
  x = x + alpha*p;
  
  % Evaluate iterate norm
  norms.x = norm(x);
  
  % Evaluate residual
  r = A*x - b;
  
  % Evaluate residual norm
  norms.r = norm(r);
  
  % Evaluate CG multiplier
  beta = (r'*Ap)/pAp;
  
%   % Print step information
%   fprintf('%.4e  %+.4e  %+.4e\n',norms.p,alpha,beta);
  
  % Update direction
  p = -r + beta*p;
  
  % Evaluate direction norm
  norms.p = norm(p);
  
  % Increment iteration counter
  k = k + 1;
  
end

% % Print output footer
% fprintf('%s\n%s\n',out_null,out_line);