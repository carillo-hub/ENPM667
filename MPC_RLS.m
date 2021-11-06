
%% Figure 4 Recreation - Create Experimental Data for RLS Param ID 
% Create Tc (pre-defined), Case1-3 curve manually (through x-y coords) 
% --------------------------------------------------------------------

% input overlay torque Tc and time info 
fc=0.3;
T_tot = 10;
delta_t = 0.01;
t=[0:delta_t:T_tot];
Tc = 30.*sin(2*pi.*fc*t);  %this is a 1x1000 array, ex. to index the last column val, do Tc(1,1001)

% output Case 1 - Hands Off - Td = 0 
x_1 = [0  1 1.25 1.5 2  2.75 2.85   3    3.1    3.2      3.25   3.75  4.25 4.35 4.5   4.75  4.9    5   6 6.1    6.25    6.5  6.75    6.8  7.1 7.5 7.7  7.9    8  8.1 8.25  9   9.25 9.5 10];
y_1 = [0 42 41   39  0 -35  -34.75 -34.5 -34.25  -34.17   -34.1   10   32   31.9 31.75  31.5 31.25 31  -31 -30.9 -30.75 -30.5 -30.25 -30   0   28.5  29  28.75 28.5 28.4 27.5  -20  -30 -30 -30];
spl_1 = spline(x_1,y_1);
theta_1 = ppval(spl_1,t);  %this is a 1x1000 array, ex. to index the last column val, do theta_1(1,1001)

% output Case 2 - Compliant Driver - Tds = 0
%x_2 = [];
%y_2 = [];
%spl_2 = spline(x_2,y_2);
%theta_2 = ppval(spl_2,t);

% output Case 1 - Hands Off - Td = 0 
%x_3 = [];
%y_3 = [];
%spl_3 = spline(x_3,y_3);
%theta_3 = ppval(spl_3,t);

% -----------------------
% Figure 4 Creation
% Plot the input and cases
% -----------------------
figure(1)
plot(x_1,y_1,'x','MarkerEdgeColor','black')
hold on
plot(t,theta_1,'r-','LineWidth',2)
plot(t,Tc,'m','LineWidth',2)
hold off
grid on;
title('Hand wheel angle θ response to the overlay torque input Tc for the test cases')
xlabel('time')
ylabel('Angular Pos. [deg] / Torque [Nm]')

%% Extract X_k from the input data
% X_k = [1, θ, w, Tc] @ k-1
% There will be a X_k vect for each Case (each θ) 1-3
% θ is fed from the arrays theta_<1,2,3> from the prev section
% w = deltaθ/delta_t 
% Tc is fed from the Tc array from the prev section
%---------------------------------------------------------------

% output Case 1 - X_k_1 creation 
w_1 = zeros(1001,1);
for ctr = 2:1:1000
    prev = ctr -1 ;
    w_1(ctr,1) = (( theta_1(1,ctr) - theta_1(1,prev) ) / delta_t );
end

X_k_1 = ones( (T_tot/delta_t)+1 , 4 );
X_k_1(:,2) = theta_1';
X_k_1(:,3) = w_1;
X_k_1(:,4) = Tc';

% output Case 2 - X_k_2 creation 


% output Case 3 - X_k_3 creation 



%% RLS Parameter Identification w/ EFRA Algorithm
% Use eqs 18-21 from Section V-B. 
% -------------------------------------------------------------------------
% curr_phi_k -- (Φk) 1x4 param vect w/ bias term & 3 terms from eq16
% phi_k_init -- set to 0 initially
% new_phi_k -- new meas of phi_k
% curr_Y_k -- lin regress model output at time step k
% new_y_k -- the new measured (output) y from the current model
% X_k -- regress vector in lin regress model -- X_k = [1, θ, w, Tc] @ k-1
% K_k -- filter gain used to tune the calc'd parameter vector (Φk)
% P_k_curr -- cov matrix for calc'ing K, has defined params α,β,γ,λ
% P_k_init -- set to arb. lrg. # init. (P_k_init = Id matrix * sigma) 
% sigma -- large # .......but idk the matrix dims yet!!
% P_k_new -- updated for every iteration of EFRA 
% -------------------------------------------------------------------------

% initial conditions
Phi_k_init = [0 0 0 0];    %given
model_Y_k_init = -1;        % I chose this randomly, but knew it couldnt = 0
sigma = 50;              % I chose but it says should be >>0
P_k_init = sigma*eye(4);   % I think cov mat should be a 4x4  

% --------------------------
% EFRA algorithm description: 
% --------------------------
% 1. Start w/ Pk_init (lrg) and Φk_init (zero)
% 2. At each time step k, calc model output Yk & filt gain K_k 
% 3. This new meas of the model output is calc'd (yk), update Φk & Pk
    
% Step 1 --initialize the variables for EFRA 
Phi_k_curr = Phi_k_init;
y_k_curr = model_Y_k_init;
P_k_curr = P_k_init;

% Step 2&3 --at each time step k calc new y_k & K_k, then new Phi_k & P_k
k = t/delta_t;
Phi_k_list_case1 = zeros(1001,4);
Phi_k_list_case2 = zeros(1001,4);
Phi_k_list_case3 = zeros(1001,4);
[row_X_k,col] = size(X_k_1);   %X_k is same size for all cases


% output Case 1 - Parameter ID'ing via EFRA RLS method 
for each = 1:1:row_X_k
    X_k_in = X_k_1(each,:);
    [y_k_new, K_k, Phi_k_new, P_k_new] = EFRA(X_k_in, P_k_curr, Phi_k_curr,y_k_curr)  ; 
   
    y_k_curr = y_k_new; 
    Phi_k_curr = Phi_k_new;
    P_k_curr = P_k_new;

    Phi_k_list_case1(each,:) = Phi_k_curr;

end


% output Case 2 - Parameter ID'ing via EFRA RLS method 


% output Case 3 - Parameter ID'ing via EFRA RLS method 




% ----------------------------------
% Plot the output parameters (Fig 5)
% ----------------------------------
figure(2)
subplot(2,2,1)
x = t;
phi0_case1 = Phi_k_list_case1(:,1);
plot(x,phi0_case1)
title('Subplot 1: Phi-0)')

subplot(2,2,2)
phi1_case1 = Phi_k_list_case1(:,2);
plot(x,phi1_case1)
title('Subplot 2: Phi-1)')

subplot(2,2,3)
phi2_case1 = Phi_k_list_case1(:,3);
plot(x,phi2_case1)
title('Subplot 3: Phi-2)')

subplot(2,2,4)
phi3_case1 = Phi_k_list_case1(:,4);
plot(x,phi3_case1)
title('Subplot 4: Phi-3)')



%% Mechanical Impedance Property Determination
% Determine J_eq, b_eq, and k_eq from eq.15/16
% From the previous section, we calc'd phi1-3
% phi-1 = -(delta_t*k_eq)/J_eq --> k_eq = -(phi1*J_eq)/delta_T
% phi-2 = 1- (delta_t*b_eq)/J_eq --> b_eq = [J_eq - (phi2*J_eq)] / delta_t
% phi-3 = delta_t/J_eq --> J_eq = delta_t/phi3
% ---------------------------------------------

% Case 1 property determination
J_eq_case1 = delta_t./phi3_case1 ;
b_eq_case1 = (J_eq_case1 - (phi2_case1.*J_eq_case1)) ./delta_t ;
k_eq_case1 = -1.*(phi1_case1.*J_eq_case1) ./delta_t ;

% Case 2 property determination


% Case 3 property determination



% ----------------------------------
% Plot the output properties (Fig 6)
% ----------------------------------
figure(3)
subplot(3,1,1)
x = t;
plot(x,J_eq_case1)
title('Subplot 1: J_eq)')
ylim([-10,10])  %added limits for viewing purposes only

subplot(3,1,2)
plot(x,b_eq_case1)
title('Subplot 2: b_eq)')
ylim([-500,0])
 
subplot(3,1,3)
plot(x,k_eq_case1)
title('Subplot 3: k_eq)')
ylim([-10,10])



%% MPC Controller Design 
% Use the cost equation (24a) and constraints in (24c,d,e) to minimize Tc 
% & slack e while tracking the desired trajectory state Z (go straight) 
% ------------------------------------------------------------------------
% 1. Using initial state @ current time step k --> Z_k and u_k_prev
%    Note: Z_k = model curr state, u_k_prev = the last Tc applied to Z  
% 2. Recalculate sequence u_k = {u_0, u_1, u_2....u_N} for the pred. horiz.
%    Note: prediction horizon is N = 12 steps, u ~ Tc applied
% 3. Input this u_k seq to Z_ik = f(Z_i,u_i,c_i) 
%    Note: u_k ~ u_i in Z_ik --> is the ideal Tc seq for car = straight
%    Note: c_i = road curvature; for our experiment, road is straight,
%          so we take c_i ~ 0 always. 
% 4. Re-eval the cost eq --> minimum(u_k, e){fcn(u_ik,Z_ik,e)}
% 5. Repeat this process every ΔTc = 0.1s 
% ------------------------------------------------------------------------
% Some given variables...
% State reference variable Z_r = 0 b/c straight road w/ curv. c = 0
% M = 1000 to penalize slack e --> greater safety margin
% (||z||_G)^2 = [z.T][G][z] (matrix multiplication)
% Q & R = weights on Tc and ΔTc for smoothing (take Q , R = 15,10) 
% Weight W=diag(0,5,5,0,10,0,0,0) to track desired vehicle state 
% ΔTmpc = 0.2s for recalc'ing discrete model from contin model Z 
% --AKA, recalc'ing mech. impedance params phi-1,2,3
% ------------------------------------------------------------------------
% Safety constraints/variables...
% ey = lateral vehicle pos on road, ey = fcn(wlane,wc,dsafe)
% ey = [-0.67m, 4.07m] --> "american" road -- favors the right side
% width of road wlane = 3.4m, width of car wc = 1.86m 
% safety margin dsafe = 0.1m 
% ------------------------------------------------------------------------
% Input constraints/variables...
% Tc = [-5 Nm,5 Nm] ,  ΔTc = [-10 Nm/s , 10 Nm/s]
% ------------------------------------------------------------------------
% discrete-time aug. pred. model --> MPC controller 
% ζk+1 = f( ζk , uk , κk )
% ζ =[Δx,Δy, ̇Δψ, e_̇ψ , ey , s, θ, ω].transpose
% u_k = input Tc for that time step (from seq u_k)
% κk = c_k = road curvature = 0
% Δx = longitudinal velociy = ~11m/s 
% Δy = lateral velocity
% ̇Δψ = yaw rate about z axis 
% e_ψ = angle of vehicle wrt centerline tangent (heading angle)
% ey = vehicle lateral position wrt centerline = -1.5m initially
% s = longitudinal position of vehicle along lane centerline
% θ = steering wheel angle
% ω = rate of steering wheel angle change 
% ------------------------------------------------------------------------
% NOTE -- the lane boundaries, curvature of road and the heading of
% the vehicle WRT the lane are normally given by the vehicle vision sys.
% THEREFORE -- we need to FAKE the ey graph and e_ψ
% With this data we can get Δy (delta ey) and Δψ (delta e_ψ)
% ------------------------------------------------------------------------
% vehicle model variables
t_mpc = 0.2 %how often we estim the vehicle state 
m = 1830    %mass of car 
lf = 1.152  %length from front axle to CoG
lr = 1.693  %length from rear axle to CoG
Iz = 3477   %Car inertia (rotational ab Z axis)
Caf = 40703 %front tire cornering stiffness
Car = 64495 %rear tire cornering stiffness
wc = 1.86   %width of car
Ns = 14.5   %gear ratio
% ------------------------------------------------------------------------
% MPC controller variables
t_c = 0.1   %how often we solve the cost eq
t_tot = 10  %how long the test is for 
N = 12      %pred horizon
Q = 15      % penalty for u in cost eq
R = 10      % penalty for delta_u in cost eq
M = 1000    % slack penalty 
Tmin = -5
Tmax = 5
ey_min = 0.67
ey_max = 4.07
wlane = 3.4 
dsafe = 0.1 
dTmax = 10
dTmin = -10
W = diag(0,5,5,0,10,0,0,0)
% ------------------------------------------------------------------------
% Initial Conditions
U_init = ones(1,N)
Z_init = [11, 0, 0, 0, -1.5, 0, 0, 0]
Phi_init = [0, 0, 0, 0.2]

% Discretization of the continuous-time dynamics
% ζ =[Δx,Δy, ̇Δψ, e_̇Δψ , ey , s, θ, ω].transpose
[theta_pred, w_pred] = discretize(k_eq,b_eq,J_eq,theta_curr, w_curr, Tc_curr)

% Cost equation
Switch = 1
for time = 0:t_c:t_tot
    Tc_term = 0
    Zterm = 0
    Cost = 0

    for i = 1:1:N-1
        Tstep = U(i)'*Q*U(i) + dU(i)'*R*dU(i) 
        Tc_term = Tc_term + Tstep
    end 
    for j = 1:1:N
        Zstep = Z(j)'*W*Z(j)  
        Z_term = Z_term + Zstep
    end 

    New_Cost = Tterm + Zterm + (M*e)

    if Switch < 0
        Z_curr = Zold
        Z_curr = discretize(Zold)  %every 0.2s we have to re-estim the discrete-time dynamics model

    Switch = Switch.*-1

end 
