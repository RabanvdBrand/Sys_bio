%% COMPUTER PRACTICAL  

%% IMPLEMENTATION BASIC CANCER MODEL 
% Set the time of the simulation
t0 = 0; tf = 20;  %(day)

% Set the initial condition
% C, D
x0 = [49.0497;0.171; 2; 0.9; 0.0019];

% Solve the ODE system
%opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[T,X] = ode15s(@Project1ODE,[t0 tf],x0,[]); %Project1ODE is the function that contains the ODE model 
%T saves all time steps; 
% X saves all solutions (first column C, second column D)

%figure 2 plots
figure;
subplot(1,5,1);
plot(T,X(:,1));
xlabel('time (day)');
ylabel('cancer cell volume (mm^3)');

subplot(1,5,2);
plot(T,X(:,2));
xlabel('time (day)');
ylabel('chemotherapeutic drug concentration D (mg/ml');

subplot(1,5,3);
plot(T,X(:,3));
xlabel('time (day)');
ylabel('endothelial cell density(mm^3)');

subplot(1,5,4);
plot(T,X(:,4));
xlabel('time (day)');
ylabel('VEGF concentration (mg/ml)');

subplot(1,5,5);
plot(T,X(:,5));
xlabel('time (day)');
ylabel('Anti VEGF (mg/ml)');