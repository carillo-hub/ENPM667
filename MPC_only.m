% ζ =[Δx,Δy, ̇Δψ, e_̇ψ , ey , s, θ, ω].transpose

% initial conditions
global Tfin N x0 xfin tmpc time 
Tfin = 10;
N = 12;
Z0 = [11,0, 0, 0, -1.5, 0, 0, 0];
phi_init = [0, 0, 0, 0.2];
tmpc = 0.1;
time = N*tmpc;

% Optimal Controller I.C 
seq_u_e_init = ones(1,N+1);
seq_u_init = seq_u_e_init(1,1:N);
winit = 0; thetainit = 0;
[w0, theta0] = Steering(seq_u_init, winit, thetainit, phi_init);
Steer_seq0 = [w0; theta0];
Z_seq_init = Z_pred(Steer_seq0, Z0')

%constraints 
lb = zeros(1,N+1); ub = zeros(1,N+1);
lb(1,1:N) = -5;
ub(1,1:N) = 5; 
lb(1,N+1) = 0;
ub(1,N+1) = inf ;
global dTmin dTmax
dTmax = 10;
dTmax = 10;


%% Optimal Control problem --> Minimization
%-------------------------------------------------
nonlcon = @dU_bounds ;
u_e_optimal = fmincon(@(seq_u_e)Minim(Z_seq_init, seq_u_e(:,1:N), seq_u_e(:,N+1)), seq_u_e_init, [], [], [], [], lb, ub, nonlcon );
u_opt = u_e_optimal(:,1:N)
e_opt = u_e_optimal(:,N+1)
[w_opt, theta_opt] = Steering(u_optimal, winit, thetainit, phi_init);
Steer_seq = [w_opt; theta_opt]
Z_opt = Z_pred(Steer_seq, Z0')
Cost = Minim(Z_opt, u_optimal, e)

% Plot results 
figure(1)
t=[0:tmpc:tmpc*(N-1)];
plot(t, Z_opt(5,:),'x','MarkerEdgeColor','black')
grid on 
hold on
plot(t,u_opt)
hold on 
scatter(t,e_opt,'filled')
title('Lateral deviation from lane center line')
xlabel('time')
ylabel('ey')
ylim([-2,2])
legend('ey','u_seq','e')
hold off 

% non linear constraints 
function [c,ceq] = dU_bounds(seq_u_e)
global N dTmin dTmax 
% Generate dU from seq_u
for ctr = 2:1:N
    dU(:,ctr) = seq_u_e(:,ctr) - seq_u_e(:,ctr-1); 
end
c1 = dU - dTmax;
c2 = dTmin - dU;
c = c1 + c2;
ceq = [];

end 

%% SS to describe the system --> outputs the Steering theta and w sequence for u_seq
%-------------------------------------------------------------------------
function [w, theta] = Steering(seq_u, wcurr, thetacurr, phi_list)
global tmpc time N
phi1 = phi_list(1,2);
phi2 = phi_list(1,3);
phi3 = phi_list(1,4);
theta = ones(1,N); w = ones(1,N);                  %initialize the output arrays

for ctr1 = 1:1:N
    u = seq_u(:,ctr1);                             %get first elem of u seq
    wnew = phi1*thetacurr + phi2*wcurr +phi3*u ;   %use it to find wnew
    thetanew = thetacurr + wnew*tmpc;              %use wnew to find thetanew
    w(:,ctr1) = wnew;                              %get the w_seq for u_seq
    theta(:,ctr1) = thetanew;                      %get the theta_seq for u_seq
    thetacurr = thetanew;                          %reset thetacurr for next iter
    wcurr = wnew;                                  %reset wcurr for next iter
end

end 

%% Generate entire Z sequence for N prediction horizon length 
% Requires an input sequence u and the current state Z_curr 
%  ζ =[Δx,Δy, ̇Δψ, e_̇ψ , ey , s, θ, ω].transpose
%-----------------------------------------------------------
function Z_seq = Z_pred(Steer_seq, Z_curr)

% Define some variabes 
%ζ =[Δx,Δy, ̇Δψ, e_̇ψ , ey , s, θ, ω].transpose
global tmpc N time 
Ns = 14.5;
Z_seq = zeros(8,N);
xdot = ones(1,N) * 11;
m = 1830; lf = 1.152; Caf = 40703; Car = 64495 ; lr = 1.693;
xdot = Z_curr(1,:); ydot = Z_curr(2,:);
cactusdot = Z_curr(3,:); ecactus = Z_curr(4,:);
ey = Z_curr(5,:); s = Z_curr(6,:); head_f = Z_curr(7,:) ;
w = Steer_seq(1,:);
theta = Steer_seq(2,:);

% Generate the rest of the Z array from Zcurr and theta array 
for ctr3 = 1:1:N
    global tmpc
    head_f = head_f/Ns ;
    ydotdot = (1/m)*(-1*m*-1*ydot/lf)*(xdot) + 2*(1/m)*( (-Caf*cos(head_f)*( (ydot/xdot) + (lf/xdot*(-1*ydot/lf)) - head_f) ) - (Car*((ydot/xdot) - (lr/xdot*(-1*ydot/lf)))));
    ydot_next = ydotdot*tmpc + ydot ;
    cactusdot_next = -1*ydot_next/lf;
    ecactus_next = cactusdot*tmpc + ecactus;
    eydot = (xdot*sin(ecactus)) + (ydot*cos(ecactus));
    ey_next = eydot*tmpc + ey;
    sdot = (xdot*cos(ecactus)) - (ydot*sin(ecactus));
    s_next = sdot*tmpc + s;
    head_f_next = theta(:,ctr3) ;
    w_next = w(:,ctr3);

    Z_seq(:,ctr3) = [xdot; ydot_next; cactusdot_next; ecactus_next; ey_next; s_next; head_f_next; w_next;];

    ydot = ydot_next;
    cactusdot = cactusdot_next;
    ecactus = ecactus_next;
    ey = ey_next;
    s = s_next;
end

end


%%    
% Re-evaluate State prediction model (discrete-time dynamics model) 
%if Switch < 0;
%    Z_curr = Zold;
%    Z_curr = discretize(Zold);   
%end


%%  Generate Cost 
% Requires an input sequence u and the current state Z_curr 
%  ζ =[Δx,Δy, ̇Δψ, e_̇ψ , ey , s, θ, ω].transpose
%-----------------------------------------------------------

function Cost = Minim(Z_seq, seq_u, e)

% Given Variables 
N = 12; 
Q = 15; R = 10; M = 1000;
W = zeros(8,8); W(2,2) = 5; W(3,3) = 5; W(5,5) = 10; % Weight matrix for penalizing items in ζ 
Tc_term = 0; Tstep = 0; Z_term = 0; Zstep = 0;
dU = zeros(1,N); Zref = [0; 0; 0; 0; 0; 0; 0; 0;];

% Generate dU from seq_u
for ctr = 2:1:N
    dU(:,ctr) = seq_u(:,ctr) - seq_u(:,ctr-1); 
end
% generate Tc_term in Cost eq
for i = 1:1:N-1
    Tstep = seq_u(:,i)'*Q*seq_u(:,i) + dU(:,i)'*R*dU(:,i); 
    Tc_term = Tc_term + Tstep;
end
% generate Z_term in Cost eq
for j = 1:1:N
    Zstep = (Z_seq(:,j) - Zref)'*W*(Z_seq(:,j) - Zref)  ;
    Z_term = Z_term + Zstep;
end 
% Generate Final Cost 
Cost = Tc_term + Z_term + (M*e);
end