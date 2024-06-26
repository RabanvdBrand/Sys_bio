function dYdt = Project1ODE(t,Y,p)
% ODE cancer model, including cancer cells (C), chemotherapy drug (D),
% endothelial cells (E), VEGF (V), anti-VEGF (A)
%
% input: 
% - t:  time
% - Y:  vector of state variables (C, D)
% - dYdt:   column vector of derivatives

%% Parameters
C_m = 10^4; % Tumour carrying capacity
d_c = 0.103; % Natural death rate tumour cells (daily)
b = 0.1685; % death rate of cancer cells due to chemotherapy
P_ce = 0.4579; % proliferation rate of tumour cells
d_d = 0.1825; % clearance of therapeutic agent
b_k = 1.0839*10^-6; % clearance of chemotherapeutic agents due to binding
E_thres = 4.5; % endothelial threshold
P_e = 0.03; % proliferation rate of endothelial cells
E_m = 10; % carrying capacity endothelial cells
d_e = 0.05; % natural death rate of endothelial cells
d_v = 0.1; % decay rate of VEGF
d_b = 0.05; % clearance of anti-VEGF drug
P_v0 = 10^-3; % baseline secretion by tumour cells
k = 1; % rate of VEGF  induced proliferation
U_b = 10; % removal of VEGF through binding
U_bk = 0.15; % clearance of the anti VEGF drug

k_c = 0.06;

%% Variables 
C = Y(1);
D = Y(2);
E = Y(3);
V = Y(4);
A = Y(5);
if E < E_thres
    P_v0 = 20^-3;
end
%  * E * (1-E/E_m)
%
%% ODEs 
dAdt = -d_b*A-U_bk*V*A;
dVdt = P_v0 - d_v * V - U_b*V*A;
dEdt = k * V + P_e * E * (1-(E/E_m))- d_e * E;
dCdt = ((P_ce * E)/(E_m + E)) * (1-(C/C_m)) - d_c * C - b * C * D/(k_c + D);
dDdt = -d_d*D-b_k*C*D/(k_c+D);

dYdt = [dCdt; dDdt; dEdt; dVdt;  dAdt];
end


