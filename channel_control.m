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

B_wall = [B(1:end, 1) zeros(238,117) B(1:end, N-2) B(1:end, N-1) ...
    zeros(238,117) B(1:end, 2*N-4) B(1:end, 2*N-3) zeros(238,117) ...
    B(1:end, 3*N-6)];
% crop original B matrix right before icare using index notation for the
% desired rows

D = zeros(size(C,1),size(B,2));
sys = ss(A,B,C,0);

% We will simulate the linearised Navier-Stokes equations that are forced
% by random noise. 

x0 = randn(398,1);
t = [0:0.2:500];
u = randn(length(t),size(B,2));
G = lsim(sys,u,t);

% G contains the velocity field of the system forced by noise. Now, we
% plot the variance of the velocity field. These would also be the Reynolds
% stresses of the turbulence for this wavenumber pair. 

for i=1:length(t)
    uu(i,:) = G(i,1:N-2).*conj(G(i,1:N-2));
    vv(i,:) = G(i,N-1:2*N-4).*conj(G(i,N-1:2*N-4));
    ww(i,:) = G(i,2*N-3:3*N-6).*conj(G(i,2*N-3:3*N-6));
end

% Let's now design the controller! We need the Q matrix which penalises
% the energy of the states and the R matrix which penalises the control

Q = 0.1*eye(2*N-4); % first number indicates how strong the controller is
R = 0.1*eye(3*N-6);

% Now to design the controller that solves the algebraic Riccati equation,
% we call the command icare

[X,K,L] = icare(A,B,Q,R,[],[],[]);

% The closed-loop system is A-B*K which results in a new state-space model
% and we force it with the same noise that we used for the original system

sys2 = ss(A-B*K,B,C,[]);
G2 = lsim(sys2,u,t);

% Compute the variance of the velocity
for i=1:length(t)
    uu2(i,:) = G2(i,1:N-2).*conj(G2(i,1:N-2));
    vv2(i,:) = G2(i,N-1:2*N-4).*conj(G2(i,N-1:2*N-4));
    ww2(i,:) = G2(i,2*N-3:3*N-6).*conj(G2(i,2*N-3:3*N-6));
end

Qp = 0.00001*eye(2*N-4); % first number indicates how strong the controller is
[X,K,L] = icare(A,B,Qp,R,[],[],[]);
sys3 = ss(A-B*K,B,C,[]);
G3 = lsim(sys3,u,t);
for i=1:length(t)
    uu3(i,:) = G3(i,1:N-2).*conj(G3(i,1:N-2));
    vv3(i,:) = G3(i,N-1:2*N-4).*conj(G3(i,N-1:2*N-4));
    ww3(i,:) = G3(i,2*N-3:3*N-6).*conj(G3(i,2*N-3:3*N-6));
end

% QUESTION 5

[X,K,L] = icare(A,B_wall,Q,R,[],[],[]);
sysb = ss(A-B*K,B,C,[]);
Gb = lsim(sysb,u,t);
for i=1:length(t)
    uub(i,:) = Gb(i,1:N-2).*conj(Gb(i,1:N-2));
    vvb(i,:) = Gb(i,N-1:2*N-4).*conj(Gb(i,N-1:2*N-4));
    wwb(i,:) = Gb(i,2*N-3:3*N-6).*conj(Gb(i,2*N-3:3*N-6));
end

Qw = 1000*eye(2*N-4); % first number indicates how strong the controller is
[X,K,L] = icare(A,B_wall,Qw,R,[],[],[]);
sysw = ss(A-B*K,B,C,[]);
Gw = lsim(sysw,u,t);
for i=1:length(t)
    uuw(i,:) = Gw(i,1:N-2).*conj(Gw(i,1:N-2));
    vvw(i,:) = Gw(i,N-1:2*N-4).*conj(Gw(i,N-1:2*N-4));
    www(i,:) = Gw(i,2*N-3:3*N-6).*conj(Gw(i,2*N-3:3*N-6));
end

% Let's plot a random wall-normal location first to see the streamwise
% velocity at this location
figure(1)
plot(real(G(:,40)))
hold on
plot(real(G2(:,40)))
hold on
plot(real(G3(:,40)))

% Now let's plot the variance profiles where we see that the variance of
% the streamwise velocity is much lower in the controlled flow which means
% that the cost function has successfully penalised peturbations.

figure(2)
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

figure(3)
subplot(1,3,1)
title('U')
plot(y(2:end-1),mean(uu,1),'Linewidth',2)
hold on
plot(y(2:end-1),mean(uub,1),'--','Linewidth',2)
hold on
plot(y(2:end-1),mean(uuw,1),'Linewidth',2)
legend({'No Controller','Wall Controller','Wall Controller Amplified'},'Location','northeast')
subplot(1,3,2)
title('V')
plot(y(2:end-1),mean(vv,1),'Linewidth',2)
hold on
plot(y(2:end-1),mean(vvb,1),'--','Linewidth',2)
hold on
plot(y(2:end-1),mean(vvw,1),'Linewidth',2)
subplot(1,3,3)
title('W')
plot(y(2:end-1),mean(ww,1),'Linewidth',2)
hold on
plot(y(2:end-1),mean(wwb,1),'--','Linewidth',2)
hold on
plot(y(2:end-1),mean(www,1),'Linewidth',2)

