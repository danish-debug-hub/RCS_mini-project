%% =========================================================
%  PI-AQM vs H∞ Static Output Feedback Comparison
%  with properly tuned PI gains and disturbance rejection analysis
% =========================================================

clear; clc; close all;

%% --- Network parameters ---
N = 1200;
C0 = (100e6)/(500*8);  % Capacity in packets/sec
Tp0 = 0.2;             % Propagation delay (s)
q0 = 500;              % Target queue length (packets)

tau0 = q0/C0 + Tp0;    % Round-trip time at equilibrium
W0 = (C0*tau0)/N;      % Window size at equilibrium
p0 = 2/(W0^2);         % Equilibrium drop probability

% H∞ SOF gain
F1 = 1e-5;

%% --- Compute PI gains using model-based tuning ---
% Linearized TCP/AQM model around equilibrium:
% G(s) = (R*C^3)/(2N^2) * e^{-sR} / (s + 2N/(R^2*C))

R0 = tau0;  % Round-trip time at equilibrium
C = C0;     % Link capacity

% Plant transfer function parameters
K = (R0 * C^3) / (2 * N^2);
a = 2*N / (R0^2 * C);

% Method 1: Ziegler-Nichols-like tuning for AQM
% Based on the dominant time constant and gain
% Gain margin approach
phase_margin_target = 60; % degrees
wc_desired = 2*pi/(R0);   % Cross-over frequency around 1/R0

% Calculate initial proportional gain based on plant gain at wc
mag_at_wc = abs(K/(a + 1j*wc_desired)); % Approximate magnitude
Kp_initial = 0.6 / mag_at_wc;  % Conservative gain

% Integral gain: Ki = Kp * (a_crossover/10) for good phase margin
Ki_initial = Kp_initial * wc_desired / 10;

% Method 2: Pole placement approach (more systematic)
% Desired closed-loop poles for PI control
zeta = 0.7;          % Damping ratio
wn = 2/R0;           % Natural frequency (~2/R0 for reasonable response)

% For PI controller C(s) = Kp + Ki/s = Kp*(1 + 1/(Ti*s))
% Closed-loop characteristic equation: 1 + K*Kp*(1 + 1/(Ti*s))/(s+a) = 0
% Solve for Kp and Ki
Ti = 2*zeta/wn - a/(K*wn^2);  % Reset time
if Ti > 0
    Kp_pp = (2*zeta*wn - a)/K;
    Ki_pp = Kp_pp/Ti;
else
    % Fall back to simpler tuning if pole placement fails
    Kp_pp = Kp_initial;
    Ki_pp = Ki_initial;
end

% Choose the more conservative gains for stability
Kp = min(Kp_initial, Kp_pp) * 0.8;  % Add safety margin
Ki = min(Ki_initial, Ki_pp) * 0.8;

% Alternative: Use internal model control (IMC) tuning for first-order with delay
lambda = R0;  % Closed-loop time constant (set equal to RTT for reasonable response)
Kp_imc = (2*lambda + R0*a)/(K*R0);
Ki_imc = Kp_imc/(lambda + R0/2);

% Use IMC tuning as it handles delay better
if Kp_imc > 0 && Ki_imc > 0
    Kp = Kp_imc * 0.7;  % Conservative scaling
    Ki = Ki_imc * 0.7;
end

% Ensure gains are reasonable
Kp = max(1e-10, min(Kp, 1e-2));  % Bound between reasonable values
Ki = max(1e-10, min(Ki, 1e-2));

fprintf('Computed PI gains:\n');
fprintf('  Kp = %e\n', Kp);
fprintf('  Ki = %e\n', Ki);
fprintf('  Equilibrium p0 = %f\n', p0);
fprintf('  RTT (tau0) = %f s\n', tau0);
fprintf('  Window size (W0) = %f packets\n', W0);

%% --- Simulation parameters ---
T = 20;               % Increased simulation time for disturbance analysis
dt = 5e-4;            % Time step (s)
steps = floor(T/dt);
time = (0:steps-1)*dt;

%% --- Disturbance parameters ---
% Define disturbance scenarios
disturbance_type = 'capacity_change';  % Options: 'load_change', 'capacity_change'

switch disturbance_type
    case 'load_change'
        % Step change in number of TCP connections
        disturbance_time = 5;  % When disturbance occurs (s)
        disturbance_duration = 2;  % Duration (s)
        N_disturbance = 1400;  % Increased number of connections (50% increase)
        disturbance_magnitude = (N_disturbance - N)/N;  % Relative change
        
    case 'capacity_change'
        % Step change in link capacity
        disturbance_time = 5;
        disturbance_duration = 5;
        C0_disturbance = C0 * 0.7;  % 30% capacity reduction
        disturbance_magnitude = (C0_disturbance - C0)/C0;

end

%% --- Initial states ---
W_h = zeros(1,steps); q_h = zeros(1,steps); p_h = zeros(1,steps);
W_p = zeros(1,steps); q_p = zeros(1,steps); p_p = zeros(1,steps);

W_h(1)=0.1; q_h(1)=0; p_h(1)=p0;
W_p(1)=0.1; q_p(1)=0; p_p(1)=p0;

int_err = 0;

%% --- Anti-windup for PI controller ---
% To prevent integrator windup when p is saturated
p_max = 1.0;  % Maximum drop probability
p_min = 0.0;  % Minimum drop probability
integral_max = (p_max - p0)/Ki;  % Max integral contribution
integral_min = (p_min - p0)/Ki;  % Min integral contribution

%% --- Disturbance signal ---
disturbance = zeros(1, steps);
dist_start_idx = round(disturbance_time/dt);
dist_end_idx = round((disturbance_time + disturbance_duration)/dt);

switch disturbance_type
    case 'load_change'
        % Apply disturbance by changing number of connections
        N_current = N * ones(1, steps);
        N_current(dist_start_idx:dist_end_idx) = N_disturbance;
        
    case 'capacity_change'
        % Apply disturbance by changing capacity
        C_current = C0 * ones(1, steps);
        C_current(dist_start_idx:dist_end_idx) = C0_disturbance;
        
    case 'impulse'
        % Add impulse to queue
        impulse_end_idx = round((disturbance_time + impulse_duration)/dt);
        disturbance(dist_start_idx:impulse_end_idx) = impulse_magnitude/impulse_duration;
end

%% --- Simulation loop ---
for k = 1:steps-1
    % Apply current disturbance parameters
    switch disturbance_type
        case 'load_change'
            N_curr = N_current(k);
            C_curr = C0;
        case 'capacity_change'
            N_curr = N;
            C_curr = C_current(k);
        otherwise
            N_curr = N;
            C_curr = C0;
    end
    
    tau = max(0.01, q_h(k)/C_curr + Tp0);
    d_idx = max(1, k - round(tau/dt));

    % ---- H∞ SOF ----
    p_h(k+1) = max(p_min, min(p_max, p0 + F1*(q_h(k)-q0)));

    dW = (1/tau) - (W_h(k)/2)*(W_h(d_idx)/tau)*p_h(d_idx);
    dq = (N_curr*W_h(k))/tau - C_curr;
    
    % Add impulse disturbance if applicable
    if strcmp(disturbance_type, 'impulse')
        dq = dq + disturbance(k);
    end

    if q_h(k)<=0 && dq<0, dq=0; end
    W_h(k+1)=max(0,W_h(k)+dW*dt);
    q_h(k+1)=min(800,max(0,q_h(k)+dq*dt));

    % ---- PI-AQM ----
    tau_p = max(0.01, q_p(k)/C_curr + Tp0);
    d_idx_p = max(1, k - round(tau_p/dt));
    
    err = q_p(k) - q0;
    
    % Anti-windup: only integrate when not saturated
    p_temp = p0 + Kp*err + Ki*int_err;
    if (p_temp <= p_max && p_temp >= p_min)
        int_err = int_err + err*dt;
    end
    
    % Bound integral term
    int_err = max(integral_min, min(integral_max, int_err));
    
    p_p(k+1) = max(p_min, min(p_max, p0 + Kp*err + Ki*int_err));

    dW = (1/tau_p) - (W_p(k)/2)*(W_p(d_idx_p)/tau_p)*p_p(d_idx_p);
    dq_p = (N_curr*W_p(k))/tau_p - C_curr;
    
    % Add impulse disturbance if applicable
    if strcmp(disturbance_type, 'impulse')
        dq_p = dq_p + disturbance(k);
    end

    if q_p(k)<=0 && dq_p<0, dq_p=0; end
    W_p(k+1)=max(0,W_p(k)+dW*dt);
    q_p(k+1)=min(800,max(0,q_p(k)+dq_p*dt));
end

%% --- Performance metrics ---
% Calculate ISE (Integral of Squared Error)
ise_h = sum((q_h - q0).^2) * dt;
ise_p = sum((q_p - q0).^2) * dt;

% Calculate IAE (Integral of Absolute Error)
iae_h = sum(abs(q_h - q0)) * dt;
iae_p = sum(abs(q_p - q0)) * dt;

% Calculate settling time (within 5% of target)
settling_threshold = 0.05 * q0;
settle_idx_h = find(abs(q_h - q0) <= settling_threshold, 1, 'first');
settle_idx_p = find(abs(q_p - q0) <= settling_threshold, 1, 'first');
settle_h = settle_idx_h * dt;
settle_p = settle_idx_p * dt;

% Disturbance rejection metrics (after disturbance starts)
dist_idx = dist_start_idx;
post_dist_steps = steps - dist_idx;
ise_h_dist = sum((q_h(dist_idx:end) - q0).^2) * dt / post_dist_steps;
ise_p_dist = sum((q_p(dist_idx:end) - q0).^2) * dt / post_dist_steps;

% Peak deviation during disturbance
max_dev_h = max(abs(q_h(dist_idx:dist_end_idx) - q0));
max_dev_p = max(abs(q_p(dist_idx:dist_end_idx) - q0));

% Recovery time (time to return to within 5% after disturbance ends)
recovery_start = dist_end_idx;
recovery_idx_h = find(abs(q_h(recovery_start:end) - q0) <= settling_threshold, 1, 'first');
recovery_idx_p = find(abs(q_p(recovery_start:end) - q0) <= settling_threshold, 1, 'first');
recovery_time_h = recovery_idx_h * dt;
recovery_time_p = recovery_idx_p * dt;

fprintf('\nPerformance Comparison:\n');
fprintf('  H∞ SOF: ISE = %.2f, IAE = %.2f, Settling time = %.3f s\n', ...
    ise_h, iae_h, settle_h);
fprintf('  PI-AQM: ISE = %.2f, IAE = %.2f, Settling time = %.3f s\n', ...
    ise_p, iae_p, settle_p);

fprintf('\nDisturbance Rejection Performance (%s):\n', disturbance_type);
fprintf('  H∞ SOF: Peak deviation = %.2f packets, Recovery time = %.3f s, Normalized ISE = %.4f\n', ...
    max_dev_h, recovery_time_h, ise_h_dist);
fprintf('  PI-AQM: Peak deviation = %.2f packets, Recovery time = %.3f s, Normalized ISE = %.4f\n', ...
    max_dev_p, recovery_time_p, ise_p_dist);

%% --- Main comparison plot ---
figure('Color','w', 'Position', [100, 100, 1200, 800]);

% Queue length comparison
subplot(2,1,1);
plot(time,q_h,'b','LineWidth',1.5); hold on;
plot(time,q_p,'r--','LineWidth',1.5);
yline(q0,'k:','Target','LineWidth',1.5);

% Mark disturbance region
if disturbance_duration > 0
    fill([disturbance_time, disturbance_time+disturbance_duration, ...
          disturbance_time+disturbance_duration, disturbance_time], ...
         [0, 0, 800, 800], [0.9 0.9 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    text(disturbance_time + disturbance_duration/2, 50, ...
        sprintf('Disturbance\n%.1f%%', disturbance_magnitude*100), ...
        'HorizontalAlignment', 'center', 'BackgroundColor', 'white');
end

legend('H∞ SOF','PI-AQM','Target','Location','Best');
xlabel('Time (s)');
ylabel('Queue Length (packets)');
title(sprintf('Queue Length Comparison\nPI gains: Kp=%.2e, Ki=%.2e', Kp, Ki));
grid on;
ylim([0, 800]);



% Window size comparison
subplot(2,1,2);
plot(time,W_h,'b','LineWidth',1.5); hold on;
plot(time,W_p,'r--','LineWidth',1.5);
yline(W0,'k:','Equilibrium W0','LineWidth',1.5);

legend('H∞ SOF','PI-AQM','Equilibrium','Location','Best');
xlabel('Time (s)');
ylabel('Window Size (packets)');
title('Window Size Comparison');
grid on;




%% --- Disturbance Rejection Comparison Plot ---
figure('Color','w', 'Position', [100, 100, 1000, 600]);

% Focus on disturbance region with zoom
zoom_start = max(0, disturbance_time - 2);
zoom_end = min(T, disturbance_time + disturbance_duration + 8);
zoom_idx = time >= zoom_start & time <= zoom_end;

subplot(2,1,1);
plot(time(zoom_idx), q_h(zoom_idx), 'b-', 'LineWidth', 2); hold on;
plot(time(zoom_idx), q_p(zoom_idx), 'r--', 'LineWidth', 2);
yline(q0, 'k:', 'Target', 'LineWidth', 1.5);

% Mark disturbance region
fill([disturbance_time, disturbance_time+disturbance_duration, ...
      disturbance_time+disturbance_duration, disturbance_time], ...
     [0, 0, 800, 800], [0.9 0.9 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

legend('H∞ SOF', 'PI-AQM', 'Target', 'Location', 'Best');
xlabel('Time (s)');
ylabel('Queue Length (packets)');
title(sprintf('Disturbance Rejection: %s\n(Zoomed View)', strrep(disturbance_type, '_', ' ')));
grid on;
xlim([zoom_start, zoom_end]);

% Drop probability during disturbance
subplot(2,1,2);
plot(time(zoom_idx), p_h(zoom_idx), 'b-', 'LineWidth', 2); hold on;
plot(time(zoom_idx), p_p(zoom_idx), 'r--', 'LineWidth', 2);
yline(p0, 'k:', 'Equilibrium', 'LineWidth', 1.5);

legend('H∞ SOF', 'PI-AQM', 'Equilibrium', 'Location', 'Best');
xlabel('Time (s)');
ylabel('Drop Probability');
title('Controller Response to Disturbance');
grid on;
xlim([zoom_start, zoom_end]);

