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

tpeak = 6.25;
X = expm(A*tpeak);
[U, Sigma, V] = svds((sqrtW*X/sqrtW),5);

Abar = U'*A*V;
Bbar = U'*B;
Cbar = C*V;

D = zeros(size(C,1),size(B,2));
sys = ss(A,B,C,0);

x0 = randn(398,1);
t = [0:0.2:500];
u = randn(length(t),size(B,2));
G = lsim(sys,u,t);
for i=1:length(t)
    uu(i,:) = G(i,1:N-2).*conj(G(i,1:N-2));
    vv(i,:) = G(i,N-1:2*N-4).*conj(G(i,N-1:2*N-4));
    ww(i,:) = G(i,2*N-3:3*N-6).*conj(G(i,2*N-3:3*N-6));
end

Q = 0.1*eye(5); % first number indicates how strong the controller is
R = 0.1*eye(3*N-6);
[X,K,L] = icare(Abar,Bbar,Q,R,[],[],[]);
sys2 = ss(Abar-Bbar*K,Bbar,Cbar,[]);
G2 = lsim(sys2,u,t);
% Compute the variance of the velocity
for i=1:length(t)
    uu2(i,:) = G2(i,1:N-2).*conj(G2(i,1:N-2));
    vv2(i,:) = G2(i,N-1:2*N-4).*conj(G2(i,N-1:2*N-4));
    ww2(i,:) = G2(i,2*N-3:3*N-6).*conj(G2(i,2*N-3:3*N-6));
end

Q = 0*eye(5); % first number indicates how strong the controller is
[X,K,L] = icare(Abar,Bbar,Q,R,[],[],[]);
sys2 = ss(Abar-Bbar*K,Bbar,Cbar,[]);
G2 = lsim(sys2,u,t);
% Compute the variance of the velocity
for i=1:length(t)
    uu3(i,:) = G2(i,1:N-2).*conj(G2(i,1:N-2));
    vv3(i,:) = G2(i,N-1:2*N-4).*conj(G2(i,N-1:2*N-4));
    ww3(i,:) = G2(i,2*N-3:3*N-6).*conj(G2(i,2*N-3:3*N-6));
end

figure(4)
subplot(1,3,1)
title('U')
plot(y(2:end-1),mean(uu,1),'Linewidth',2)
hold on
plot(y(2:end-1),mean(uu2,1),'Linewidth',2)
hold on
plot(y(2:end-1),mean(uu3,1),'Linewidth',2)
legend({'No Controller','Controller','Penalised Controller'},'Location','northeast')
subplot(1,3,2)
title('V')
plot(y(2:end-1),mean(vv,1),'Linewidth',2)
hold on
plot(y(2:end-1),mean(vv2,1),'Linewidth',2)
hold on
plot(y(2:end-1),mean(vv3,1),'Linewidth',2)
subplot(1,3,3)
title('W')
plot(y(2:end-1),mean(ww,1),'Linewidth',2)
hold on
plot(y(2:end-1),mean(ww2,1),'Linewidth',2)
hold on
plot(y(2:end-1),mean(ww3,1),'Linewidth',2)


