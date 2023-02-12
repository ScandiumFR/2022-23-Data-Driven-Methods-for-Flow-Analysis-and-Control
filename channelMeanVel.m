%%  Code to calculate turbulent mean velocity profile
% Written on 06/07/2013 by Rashad Moarref
% Modified on 06/09/2014 by Mitul Luhar
% Uses Reynolds-Tiederman eddy viscosity profile

% Inputs: 
% Re: friction Reynolds number 
% N: number of wall-normal collocation points

% Outputs
% U: mean velocity profile in plus-units (i.e. scaled by u_tau)
% y: wall-normal Chebyshev collocation points
% dy: integration weights for wall normal coordinate
% D1: First differential matrix
% D2: Second differential matrix

%%
function [y,D1,D2,dy,U0] = channelMeanVel(Re,N)

kappa = 0.426; % parameter in R-T eddy-viscosity
Amp   = 25.4;  % parameter in R-T eddy-viscosity

% Coordinate system and differentiation matrices
% Based on the chebyshev methods toolbox developed by Weidemann and Reddy
[y,DM] = chebdif(N,2);  
D1 = DM(:,:,1);
D2 = DM(:,:,2);
% Integration weights
[~,dy] = clencurt(N-1); 

% Reynolds-Tiederman eddy viscosity %
nuT = 0.5*( (1 + ( (1/3)*kappa*Re*(1 - 1*y.^2).*(1 + 2*y.^2).*(1 - exp(-(1 - abs(y))*Re/Amp)) ).^2 ).^(1/2) - 1  );
% Differential
DnuT = D1*nuT;
DnuT2 = D1*DnuT;

% Mean velocity profile moment balance
LHS = diag(1+nuT)*D2 + diag(DnuT)*D1;
RHS = -Re*ones(N,1);
% Impose boundary conditions on mean velocity profile
LHS(1,:)=0;  LHS(N,:)=0;
LHS(1,1)=1; LHS(N,N)=1;
RHS(1)=0; RHS(N)=0;
% Estimate mean velocity profile
%U0 = y; 
U0 = LHS\RHS;

uv_stress = -(nuT.*(D1*U0))/Re;
duv_dy = D1*(uv_stress);
end