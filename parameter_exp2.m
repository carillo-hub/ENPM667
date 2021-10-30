
%% Figure 4 Recreation - Create Experimental Data for RLS Param ID 
% Create Tc (pre-defined), Case1-3 curve manually (through x-y coords) 
% --------------------------------------------------------------------

% input overlay torque Tc 
fc=0.3;
t=[0:0.01:10];
Tc = 30.*sin(2*pi.*fc*t);


% output Case 1 - Hands Off - Td = 0 
x_1 = [0  1 1.25 1.5 2  2.75 3  3.25 4  4.25 4.5 5   6  6.25  6.5 7 7.5 7.75 8.25 9   9.25 9.5 10];
y_1 = [0 42 41   39  0 -35  -34 -33  10 32   31  31 -31 -30  -30  0 29  28   28  -20  -30 -30 -30];
spl_1 = spline(x_1,y_1);
case1_curve = ppval(spl_1,t);


% output Case 2 - Compliant Driver - Tds = 0
x_2 = [];
y_2 = [];
spl_2 = spline(x_2,y_2);
case2_curve = ppval(spl_2,xint);


% output Case 1 - Hands Off - Td = 0 
x_3 = [];
y_3 = [];
spl_3 = spline(x_3,y_3);
case3_curve = ppval(spl_3,xint);


% Plot the input and cases
figure(1)
plot(x_1,y_1,'x','MarkerEdgeColor','black')
hold on
plot(t,case1_curve,'r-','LineWidth',2)
plot(t,Tc,'m','LineWidth',2)
hold off
grid on;
title('Hand wheel angle θ response to the overlay torque input Tc for the test cases')
xlabel('time')
ylabel('Angular Pos. [deg] / Torque [Nm]')


%% RLS Parameter Identification
% Use eqs 18-21 from Section V-B. 
% --------------------------------------

model_P_k_last = 
model_Y_k = X_k * model_P_k_last.transpose





















