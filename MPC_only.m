% ζ =[Δx,Δy, ̇Δψ, e_̇ψ , ey , s, θ, ω].transpose

% initial conditions
global Tfin N x0 xfin tmpc time 
Tfin = 10;
N = 12;
x0 = [11,0, 0, 0, -1.5, 0, 0, 0];
xfin = [11, 0, 0, 0, 0, 0, 0, 0];
e = 0;
phi_init = [0, 0.32, 1.63, 4.98];
tmpc = 0.1;
time = N*tmpc;
seq_u = ones(1,N) * 4;
winit = 0;
thetainit = 0;
lb = seq_u*-5;
ub = seq_u*5; 
%% 
%-------------------------------------------------
% Optimal Control problem --> Minimization
%-------------------------------------------------
u_optimal = fmincon(@(seq_u)Minim(Z_seq, seq_u, e), seq_u, [], [], [], [], lb, ub )
[w_opt, theta_opt] = Statespacerep(u_optimal, winit, thetainit, phi_init);
Steer_seq = [w_opt; theta_opt]
Z_opt = Z_pred(Steer_seq, x0')
Cost = Minim(Z_opt, u_optimal, e)


figure(1)
t=[0:tmpc:tmpc*(N-1)];
plot(t, Z_opt(5,:),'x','MarkerEdgeColor','black')
grid on 
title('Lateral deviation from lane center line')
xlabel('time')
ylabel('ey')

%----------------------------------------------
% Optima Control problem --> One iteration 
%----------------------------------------------
%Switch = 1;
%global Tfin x0 ;
%[w, theta] = Statespacerep(seq_u, winit, thetainit, phi_init);
%Steer_seq = [w; theta];
%Z_seq = Z_pred(Steer_seq, x0');
%Cost = Minim(Z_seq, seq_u, e);
%Switch = Switch.*-1;

%%
%-------------------------------------------------------------------------
% SS to describe the system --> outputs the theta and w sequence for u_seq
%-------------------------------------------------------------------------
function [w, theta] = Statespacerep(seq_u, wcurr, thetacurr, phi_list)
global tmpc time N
phi1 = phi_list(1,2);
phi2 = phi_list(1,3);
phi3 = phi_list(1,4);
theta = ones(1,N); w = ones(1,N);
%u = get_seq(seq_u);                       %generate first elem of u_sequence 

for ctr1 = 1:1:N
    u = seq_u(:,ctr1)
    wnew = phi1*thetacurr + phi2*wcurr +phi3*u ;   %use it to find wnew
    thetanew = thetacurr + wnew*tmpc;              %use wnew to find thetanew
    w(:,ctr1) = wnew;                        %get the w_seq for u_seq
    theta(:,ctr1) = thetanew;                %get the theta_seq for u_seq
end

end 


%----------------------
% Construct u sequence
%----------------------
function u = get_seq(seq_u)
global Tfin N time
a = N*time/Tfin ;
if a==0
    u = seq_u(:,1);
else
    u = seq_u(:,abs(ceil(a)));
end
end

%%
%-----------------------------------------------------------
% Generate entire Z sequence for N prediction horizon length
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
    s_next = sdot + s;
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


%%
%-----------------------------------------------------------
% Generate Cost 
% Requires an input sequence u and the current state Z_curr 
%  ζ =[Δx,Δy, ̇Δψ, e_̇ψ , ey , s, θ, ω].transpose
%-----------------------------------------------------------

function Cost = Minim(Z_seq, seq_u, e)

% Given Variables 
N = 12; 
Q = 15; R = 10; M = 1000;
W = zeros(8,8); W(2,2) = 5; W(3,3) = 5; W(5,5) = 10; % Weight matrix for penalizing items in ζ 
Tc_term = 0; Tstep = 0; Z_term = 0; Zstep = 0;
dU = zeros(1,N); Zref = [11; 0; 0; 0; 0; 0; 0; 0;];

% Generate dU from seq_u
for ctr = 2:1:N
    dU(:,ctr) = seq_u(:,ctr) - seq_u(:,ctr-1); 
end
dU ;

% generate Tc_term in Cost eq
for i = 1:1:N-1
    Tstep = seq_u(:,i)'*Q*seq_u(:,i) + dU(:,i)'*R*dU(:,i); 
    Tc_term = Tc_term + Tstep;
end
Tc_term ;

% generate Z_term in Cost eq
for j = 1:1:N
    Zstep = (Z_seq(:,j) - Zref)'*W*(Z_seq(:,j) - Zref)  ;
    Z_term = Z_term + Zstep;
end 
Z_term;

% Generate Final Cost 
Cost = Tc_term + Z_term + (M*e);

end
%%  