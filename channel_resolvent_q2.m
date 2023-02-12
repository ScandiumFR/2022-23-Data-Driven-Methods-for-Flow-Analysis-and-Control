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

sigm1 = zeros(1,11);
sigm2 = zeros(1,11);
sigm3 = zeros(1,11);

% for omega = 0:10
invH = (-1i*omega)*I2  + M\L; % resolvent part. change omega
[U, S, V] = svds((sqrtQ*invH)/(sqrtQ),10,'smallest');
sigm = flip(1./diag(S));
% sigm1(omega+1) = sigm(1);
% sigm2(omega+1) = sigm(2);
% sigm3(omega+1) = sigm(3);

for i = 9:-2:8
mode = i; %note

phi1 = C*1/sqrtQ*U(:,i);
psi1 = C*1/sqrtQ*V(:,i);
length = N-2;


% Plot the mode shapes 

figure(1)
subplot(1,3,1)
title('Psi U')
plot(y(2:end-1),psi1(1:length).*conj(psi1(1:length)),'LineWidth',2) %removed 'b' for line colour
hold on
subplot(1,3,2)
title('Psi V')
plot(y(2:end-1),psi1(length+1:2*length).*conj(psi1(length+1:2*length)),'r--','LineWidth',2)
hold on
subplot(1,3,3)
title('Psi W')
plot(y(2:end-1),psi1(2*length+1:3*length).*conj(psi1(2*length+1:3*length)),'k-.','LineWidth',2)
hold on

end
% end
% figure(2)
% set(gca, 'YScale', 'log')
% xlabel('\omega') 
% ylabel('Singular Values')
% hold on
% plot((0:10),sigm1);
% hold on
% plot((0:10),sigm2);
% hold on
% plot((0:10),sigm3);
% legend({'1st singular value','2nd singular value','3rd singular value'},'Location','northeast')
