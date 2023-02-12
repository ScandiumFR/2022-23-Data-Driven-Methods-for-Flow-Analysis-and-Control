% Simulate a turbulent channel flow using stochastic forcing

clear
clc

Re = 180; % Friction Reynolds number
kx = 0;   % Streamwise wavenumber
ky = 4;   % Spanwise wavenumber
N = 121;  % Number of discretisation points

% This gives us the LNS operator, mass matrix, weighting matrix and the 
% turbulent mean velocity profile

[L,M,W,U0] = get_operators(kx,ky,N,Re,0);

% This is the weighting matrix for the non-uniform grid that we are using

sqrtW = sqrtm(W); 

% These are operators we need in order to form the state-space model

k = kx^2+ky^2;
I = eye(N-2);
Z = zeros(N-2,N-2);
[y,DM] = chebdif(N,2);
DY = DM(2:N-1,2:N-1,1);

% Here we form the state-space operators: the state matrix, the output
% matrix and the control matrix
A = - M \ L;
C = (1/k)*([1i*kx*DY -1i*ky*I;...
            k*I       Z; ...
            1i*ky*DY  1i*kx*I]);
B =  inv(M)*[-1i*kx*DY, -k*I, -1i*ky*DY ; ...
               1i*ky*I,     Z,  -1i*kx*I ];

D = zeros(size(C,1),size(B,2));
sys = ss(A,B,C,0);

% Let's compute the matrix exponential as a function of time

t = [0:0.25:50];
for i=1:length(t)
    X = expm(A*t(i)); % This is the matrix exponential! Note it's expm() 
    [~, Sigma, ~] = svds((sqrtW*X/sqrtW),10); % Compute SVD of X
    Sig1(i) = Sigma(1,1); % We are interested in sigma1 only
end

figure(1)
hold on
xlabel('Time (t)','interpreter', 'latex')
ylabel('$\textup{First Singular Value of } e^{\textbf{A}t}$','interpreter', 'latex')
plot(t,Sig1,'LineWidth',2)
hold off

% Let's investigate the structures (and the corresponding disturbance) 
% that lead to the peak in transient growth

tpeak = 6.25;
X = expm(A*tpeak);
[U, Sigma, V] = svds((sqrtW*X/sqrtW),10);
psi1 = C*1/sqrtW*U(:,1);
phi1 = C*1/sqrtW*V(:,1);

proj = U'*A*U;

% Here we plot the velocity structures

figure(3)
subplot(1,3,1)
title('U')
plot(y(2:end-1),psi1(1:N-2,1).*conj(psi1(1:N-2,1)),'LineWidth',2)
hold on
subplot(1,3,2)
title('V')
plot(y(2:end-1),psi1(N-1:2*N-4,1).*conj(psi1(N-1:2*N-4,1)),'LineWidth',2)
hold on
subplot(1,3,3)
title('W')
plot(y(2:end-1),psi1(2*N-3:3*N-6,1).*conj(psi1(2*N-3:3*N-6,1)),'LineWidth',2)
hold on