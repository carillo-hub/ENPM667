
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
x_1 = [0  1 1.25 1.5 2  2.75 3  3.25 4  4.25 4.5 5   6  6.25  6.5 7 7.5 7.75 8.25 9   9.25 9.5 10];
y_1 = [0 42 41   39  0 -35  -34 -33  10 32   31  31 -31 -30  -30  0 29  28   28  -20  -30 -30 -30];
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

% Plot the input and cases
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

w_1 = zeros(1001,1);
for ctr = 2:1:1000
    prev = ctr -1 ;
    w_1(ctr,1) = (( theta_1(1,ctr) - theta_1(1,prev) ) / delta_t );
end

X_k_1 = ones( (T_tot/delta_t)+1 , 4 );
X_k_1(:,2) = theta_1';
X_k_1(:,3) = w_1;
X_k_1(:,4) = Tc';

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
X_k_init = [1, 0, 0, 0];   % I chose the "0"'s, but the "1" is given
model_Y_k_init = 0;        % I chose

sigma = 1000;               % I chose but it says should be >>0
%P_k_init = sigma*eye(4)   % I think cov mat should be a 4x4  
P_k_init = sigma;
% --------------------------
% EFRA algorithm description: 
% --------------------------
% 1. Start w/ Pk_init (lrg) and Φk_init (zero)
% 2. At each time step k, calc model output Yk & filt gain K_k 
% 3. This new meas of the model output is calc'd (yk), update Φk & Pk
    
% Step 1 --initialize the variables for EFRA 
Phi_k_curr = Phi_k_init;
X_k = X_k_init;
y_k_curr = model_Y_k_init;
P_k_curr = P_k_init;

% Step 2&3 --at each time step k calc new y_k & K_k, then new Phi_k & P_k
k = t/delta_t;
Phi_k_list = zeros(1001,4)
len_X_k = (T_tot/delta_t)+1;
for each = 1:1:len_X_k
    X_k_in = X_k_1(each,:)
    [y_k_new, K_k, Phi_k_new, P_k_new] = EFRA(X_k_in, P_k_curr, Phi_k_curr,y_k_curr)   
   
    y_k_curr = y_k_new; 
    Phi_k_curr = Phi_k_new;
    P_k_curr = P_k_new;

    Phi_k_list(each,:) = Phi_k_curr;

end




















