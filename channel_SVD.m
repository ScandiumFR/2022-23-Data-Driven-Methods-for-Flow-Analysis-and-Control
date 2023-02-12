%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resolvent analysis using the OSSQ formulation for a turbulent channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
clc

% Input parameters

Re = 180; % Friction Reynolds number
kx = 0; % Streamwise wavenumber
kz = 4; % Spanwise wavenumber
cP = 0; % The wave speed
omega = cP*kx; % The angular frequency
N = 201; % Number of wall-normal discretization points. Original 201

[y,~] = clencurt(N-1); % The grid in the wall-normal direction

% Get operators
[L,M,Q,U0] = generate_OSSQ_operators(kx,kz,N,Re,0);

k = kx^2+kz^2;
I = eye(N-2);
Z = zeros(N-2,N-2);
[y,DM] = chebdif(N,2);
[~,w] = clencurt(N-1);
Mass = diag([w(2:end-1) w(2:end-1)]);
Mass2 = diag([w(2:end-1) w(2:end-1) w(2:end-1)]);
DY = DM(2:N-1,2:N-1,1);
sqrtQ = sqrtm(Q); 

C = (1/k)*([1i*kx*DY -1i*kz*I;...
            k*I       Z; ...
            1i*kz*DY  1i*kx*I]);
        
I2 = eye(2*(N-2));

sigm1 = zeros(1,21);
sigm2 = zeros(1,21);
sigm3 = zeros(1,21);

for omega = 0:20
invH = (-1i*omega)*I2  + M\L; % resolvent part. change omega
[U, S, V] = svds((sqrtQ*invH)/(sqrtQ),20,'smallest');
sigm = flip(1./diag(S));
sigm1(omega+1) = sigm(1);
sigm2(omega+1) = sigm(2);
sigm3(omega+1) = sigm(3);

end
figure(2)
set(gca, 'YScale', 'log')
xlabel('\omega') 
ylabel('Singular Values')
hold on
plot((0:20),sigm1);
hold on
plot((0:20),sigm2);
hold on
plot((0:20),sigm3);
legend({'1st singular value','2nd singular value','3rd singular value'},'Location','northeast')
