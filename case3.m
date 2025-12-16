%% ROBUST H-INFINITY SOF CONTROL FOR TCP/AQM ROUTERS
%  Reproduction of Case 3: Input Delay (h) = State Delay (d)
%  Based on Paper: "Robust H-inf Static Output Feedback Control..."
%  UPDATED: 
%    1. Performs actual Iterative LMI Optimization for F1
%    2. Plots arranged into 2 Figures with 2 Subplots each
clear; clc; close all;

%% ---------------------------------------------------------
%  1. INITIALIZING PARAMETERS & LINEARIZATION
%  ---------------------------------------------------------
disp('1. INITIALIZING PARAMETERS...');

% Network Parameters (Section 4)
N = 1200;                   % Number of TCP flows
C0_bps = 100e6;             % Link bandwidth (100 Mbps)
PacketSize = 500 * 8;       % Packet size in bits (4000 bits)
C0 = C0_bps / PacketSize;   % Link capacity (25,000 packets/sec)
Tp_base = 0.2;              % Base propagation delay (sec)

% Target Operating Points
q0 = 500;                   % Target queue length (packets)

% Derived Equilibrium Values
% Note: W0 must satisfy W0 = C0 * RTT / N to be in equilibrium
tau0 = q0/C0 + Tp_base;     % Equilibrium RTT
W0_calc = (C0 * tau0) / N;  % Equilibrium Window 
p0 = 2 / (W0_calc^2);       % Equilibrium Prob

disp(['   Equilibrium W0: ', num2str(W0_calc)]);
disp(['   Equilibrium p0: ', num2str(p0)]);

% Linearized State Space Matrices (Equation 3)
% State x = [delta_W; delta_q]
a11 = -N / (tau0^2 * C0);
a12 = -1 / (tau0^2 * C0);
a21 = N / tau0;
a22 = -1 / tau0;
A = [a11, a12; a21, a22];

% Ad (Delayed State Matrix)
Ad = [-N / (tau0^2 * C0), 1 / (tau0^2 * C0); 
       0,                 0];
   
% Input Matrix B
b1 = -tau0 * C0^2 / (2 * N^2);
B = [b1; 0];

% Disturbance Matrix D
d1 = (tau0 - Tp_base)/(tau0^2 * C0);
d2 = -(tau0 - Tp_base)/(tau0^2 * C0);
d3 = -Tp_base/tau0;
D_mat = [d1, d2; d3, 0];

% Output Matrices
L = [0, 1];    % Controlled output z (Queue deviation)
C_out = [0, 1]; % Measured output y (Queue deviation)

%% ---------------------------------------------------------
%  2. CONTROLLER DESIGN (ITERATIVE LMI OPTIMIZATION)
%  ---------------------------------------------------------
disp(' ');
disp('2. CONTROLLER DESIGN (Solving BMI via Iteration)...');

if exist('yalmip','file') ~= 2
    warning('YALMIP is not installed. Using paper value F1 = 1e-5.');
    F1_val = 1.00e-5;
else
    % Solver Settings 
    ops = sdpsettings('solver','sdpt3','verbose',0); 

    % --- Initialization ---
    F1_val = 1.0e-6; 
    gamma_best = 1e9;
    tol = 1e-4;
    max_iter = 15;

    disp(['   Starting Iterative LMI (Max Iter: ' num2str(max_iter) ')']);
    
    for iter = 1:max_iter
        % STEP A: Fix F1, Solve for P, Q, Gamma
        P  = sdpvar(2,2,'symmetric');
        Q1 = sdpvar(2,2,'symmetric');
        g2 = sdpvar(1,1); 
        
        M_delayed = Ad + B * F1_val * C_out;
        
        % LMI: [Phi P*M P*D; M'*P -Q 0; D'*P 0 -g2*I] < 0
        Phi_11 = P*A + A'*P + Q1 + L'*L;
        LMI_StepA = [Phi_11,       P*M_delayed,    P*D_mat;
                     M_delayed'*P, -Q1,            zeros(2,2);
                     D_mat'*P,     zeros(2,2),     -g2*eye(2)] <= -1e-6*eye(6);
                 
        ConstraintsA = [LMI_StepA, P >= 1e-6*eye(2), Q1 >= 1e-6*eye(2), g2 >= 0];
        solA = optimize(ConstraintsA, g2, ops);
        
        if solA.problem ~= 0, break; end
        
        P_curr = value(P);
        
        % STEP B: Fix P, Solve for F1, Q, Gamma
        F1_var = sdpvar(1,1);
        Q1_B   = sdpvar(2,2,'symmetric');
        g2_B   = sdpvar(1,1);
        
        Term_12 = P_curr * Ad + P_curr * B * F1_var * C_out;
        Phi_11_B = P_curr*A + A'*P_curr + Q1_B + L'*L;
        
        LMI_StepB = [Phi_11_B,     Term_12,        P_curr*D_mat;
                     Term_12',    -Q1_B,           zeros(2,2);
                     D_mat'*P_curr,zeros(2,2),     -g2_B*eye(2)] <= -1e-6*eye(6);
                 
        ConstraintsB = [LMI_StepB, Q1_B >= 1e-6*eye(2), g2_B >= 0];
        solB = optimize(ConstraintsB, g2_B, ops);
        
        if solB.problem == 0
            F1_val = value(F1_var);
            g_new = sqrt(value(g2_B));
            if abs(g_new - gamma_best) < tol, break; end
            gamma_best = g_new;
        end
    end
    disp(['   OPTIMIZED F1 GAIN: ', num2str(F1_val)]);
end

%% ---------------------------------------------------------
%  3. NONLINEAR SIMULATION
%  ---------------------------------------------------------
disp(' ');
disp('3. RUNNING SIMULATION...');
T_sim = 10;             
dt = 0.0005;            
steps = floor(T_sim / dt);
time = (0:steps-1) * dt;

% Arrays
W = zeros(1, steps);
q = zeros(1, steps);
p = zeros(1, steps);
rtt_hist = zeros(1, steps);

% Initial Conditions
W(1) = 0.1; % Non-zero start to allow fluid model evolution
q(1) = 0;   
p(1) = p0;  

for k = 1:steps-1
    % A. Time Varying Delay (Random Noise)
    xi = rand - 0.5; 
    Tp_noisy = Tp_base + 0.04 * xi;
    current_tau = max(0.001, q(k)/C0 + Tp_noisy);
    rtt_hist(k) = current_tau;
    
    % B. Delayed Variable Lookup
    delay_idx = round(current_tau / dt);
    prev_k = k - delay_idx;
    
    if prev_k < 1
        W_delayed = W(1);
        p_delayed = p(1);
        tau_delayed = tau0; 
    else
        W_delayed = W(prev_k);
        p_delayed = p(prev_k);
        tau_delayed = rtt_hist(prev_k);
    end
    
    % C. CONTROLLER (Static Output Feedback) 
    % u(t) = F1*y(t) -> p(t) = p0 + F1*(q(t)-q0)
    y_current = q(k) - q0;
    u_signal = F1_val * y_current;
    
    p_control = p0 + u_signal;
    p(k+1) = max(0, min(1, p_control)); % Saturation 0 <= p <= 1
    
    % D. NONLINEAR FLUID DYNAMICS 
    if W(k) < 1e-6
        term_dec = 0;
    else
        term_dec = (W(k)/2) * (W_delayed / tau_delayed) * p_delayed;
    end
    term_inc = 1 / current_tau;
    
    dW = term_inc - term_dec;
    rate_in = (N * W(k)) / current_tau;
    dq = rate_in - C0;
    
    % Non-negative Queue Constraint
    if q(k) <= 0 && dq < 0
        dq = 0;
    end
    
    % Euler Integration
    W(k+1) = W(k) + dW * dt;
    q(k+1) = q(k) + dq * dt;
    
    % Buffer Saturation
    if q(k+1) > 800; q(k+1) = 800; end
    if q(k+1) < 0; q(k+1) = 0; end
    W(k+1) = max(0, W(k+1));
end

%% ---------------------------------------------------------
%  4. PLOTTING
%  ---------------------------------------------------------

% --- FIGURE 1: System States (Queue & Window) ---
figure('Name', 'Figure 1: System States', 'Color', 'w');

% Subplot 1.1: Queue Length
subplot(2, 1, 1);
plot(time, q, 'b', 'LineWidth', 1.5); hold on;
yline(q0, 'r--', 'Target (500)', 'LineWidth', 1.5);
title(['Queue Length (Optimized F_1 = ' num2str(F1_val, '%.2e') ')']);
ylabel('Queue (packets)');
xlabel('Time (sec)');
grid on;
ylim([0 900]); % Matches paper scale
legend('Proposed SOF', 'Target', 'Location', 'Best');

% Subplot 1.2: Window Size
subplot(2, 1, 2);
plot(time, W, 'b', 'LineWidth', 1.5); hold on;
yline(W0_calc, 'r--', 'Equilibrium', 'LineWidth', 1.5);
title('Average TCP Window Size');
ylabel('Window (packets)');
xlabel('Time (sec)');
grid on;
xlim([0 5]);% Matches paper scale
legend('Proposed SOF', 'Equilibrium', 'Location', 'Best');


% --- FIGURE 2: Control & Environment (Input & RTT) ---
figure('Name', 'Figure 2: Control & RTT', 'Color', 'w');

% Subplot 2.1: Control Input p
subplot(2, 1, 1);
plot(time, p, 'b', 'LineWidth', 1.5); hold on;
yline(p0, 'k:', 'Equilibrium p_0', 'LineWidth', 1.5);
title('Control Input p (Mark Probability)');
ylabel('Probability p');
xlabel('Time (sec)');
grid on;
ylim([0 0.25]); % Matches paper scale
legend('Proposed SOF', 'Equilibrium', 'Location', 'Best');

% Subplot 2.2: Round Trip Time (RTT)
subplot(2, 1, 2);
plot(time, rtt_hist, 'b', 'LineWidth', 0.5);
title('Round Trip Time (with Random Noise)');
ylabel('RTT (sec)'); 
xlabel('Time (sec)');
grid on;
ylim([0.16 0.25]); % Matches paper scale