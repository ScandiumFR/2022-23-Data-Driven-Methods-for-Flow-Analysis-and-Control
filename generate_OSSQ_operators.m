function [L,M,Q,U0] = generate_OSSQ_operators(kx,kz,N,Re,eddy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs:
% L: Linear operator (OS and SQ)
% M: Mass matrix
% Q: Weight matrix for enforcing energy norm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


I = eye(N-2);
Z = zeros(N-2,N-2);

[y,DM] = chebdif(N,2);

D1 = DM(2:N-1,2:N-1,1);
D2 = DM(2:N-1,2:N-1,2);
[~,D4]=cheb4c(N); 

% Turbulent mean from eddy-viscosity model
[~,~,~,~,U0] = channelMeanVel(Re,N); 
%U0 = 1-y.^2; % Uncomment to use laminar base flow
U0y = DM(:,:,1)*U0;
U0yy = DM(:,:,2)*U0;

U0 = U0(2:N-1);
U0y = U0y(2:N-1);
U0yy = U0yy(2:N-1);

% Add eddy viscosity related details

kappa = 0.426; % parameter in R-T eddy-viscosity
Amp   = 25.4;  % parameter in R-T eddy-viscosity

% Reynolds-Tiederman eddy viscosity 
nuT = 0.5/Re*( (1 + 1*( (1/3)*kappa*Re*(1 - y.^2).*(1 + 2*y.^2).*(1 - exp(-(1 - abs(y))*Re/Amp)) ).^2 ).^(1/2))+0.5/Re;
DnuT = DM(:,:,1)*nuT;
DnuT2 = DM(:,:,1)*DnuT;

% Linear operators and mass matrix

k2 = kx^2 + kz^2;
k4 = k2^2;
if eddy == 1
    L_OS = 1i*kx*diag(U0)*(k2*I-D2) + 1i*kx*diag(U0yy) + ...
        diag(nuT(2:N-1))*(D4 - 2*k2*D2 + k4*I) + ...
        1*2*diag(DnuT(2:N-1))*D1*(D2-k2*I) + ...
        1*diag(DnuT2(2:N-1))*(D2 + k2*I);
    L_SQ =  1i*kx*diag(U0) + diag(nuT(2:N-1))*(k2*I-D2) - 1*diag(DnuT(2:N-1))*D1;
else
    L_OS = 1i*kx*diag(U0)*(k2*I-D2) + 1i*kx*diag(U0yy) + 1/Re*(D4 - 2*k2*D2 + k4*I);
    L_SQ =  1i*kx*diag(U0) + 1/Re*(k2*I-D2);
end
L_C = 1*1i*kz*diag(U0y);

L = [L_OS, Z; L_C, L_SQ];
M = [(k2*I-D2), Z; Z, I];

% Generate weight matrix for energy norm

[~,IWT] = clencurt(N-1);
IWT = diag(IWT);
IW = IWT(2:end-1,2:end-1);
D = DM(:,:,1);
QvT = (D'*IWT*D/k2 + IWT);
Qv = QvT(2:end-1,2:end-1);
Qeta = IW/k2;
Q = [Qv, Z; Z, Qeta];


