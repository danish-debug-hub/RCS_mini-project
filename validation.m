%% =========================================================
%  MONTE CARLO ROBUSTNESS ANALYSIS
%  Robust H∞ Static Output Feedback (TCP/AQM)
% =========================================================

clear; clc; close all;

%% --- Network parameters (same as paper) ---
N = 1200;
C0 = (100e6)/(500*8);
Tp0 = 0.2;
q0 = 500;

tau0 = q0/C0 + Tp0;
W0 = (C0*tau0)/N;
p0 = 2/(W0^2);

% H∞ SOF gain (use optimized or paper value)
F1 = 1e-5;

%% --- Simulation parameters ---
T = 10;
dt = 5e-4;
steps = floor(T/dt);

MC = 50;   % Monte Carlo runs
q_peak = zeros(1,MC);
q_ss = zeros(1,MC);

disp('Running Monte Carlo simulations...');

for m = 1:MC
    % Random initial conditions
    W = zeros(1,steps);
    q = zeros(1,steps);
    p = zeros(1,steps);

    W(1) = 0.2*rand;
    q(1) = 200*rand;
    p(1) = p0;

    rtt = tau0*ones(1,steps);

    for k = 1:steps-1
        % Random RTT perturbation (bounded)
        Tp = Tp0 + 0.03*(rand-0.5);
        tau = max(0.01, q(k)/C0 + Tp);
        rtt(k) = tau;

        d_idx = max(1, k - round(tau/dt));

        % H∞ static output feedback
        p(k+1) = max(0,min(1, p0 + F1*(q(k)-q0)));

        % Nonlinear TCP dynamics
        dW = (1/tau) - (W(k)/2)*(W(d_idx)/tau)*p(d_idx);
        dq = (N*W(k))/tau - C0;

        if q(k)<=0 && dq<0
            dq = 0;
        end

        W(k+1) = max(0, W(k) + dW*dt);
        q(k+1) = min(800, max(0, q(k) + dq*dt));
    end

    q_peak(m) = max(q);
    q_ss(m) = mean(q(end-300:end));
end

%% --- Results visualization ---
figure('Color','w');
subplot(2,1,1);
histogram(q_peak,15);
xlabel('Peak Queue (packets)');
ylabel('Frequency');
title('Monte Carlo Peak Queue Length');
grid on;

subplot(2,1,2);
histogram(q_ss,15);
xlabel('Steady-State Queue (packets)');
ylabel('Frequency');
title('Monte Carlo Steady-State Queue');
grid on;

disp('Monte Carlo robustness analysis completed.');
