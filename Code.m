%% ======================================================================
%  Mathematical Modelling of Alzheimer's Disease: 
%  A Nonlinear Dynamical System Model of Protein Accumulation and Neuroinflammation
%  Jadean Govender (27836894) and Nala Steward (27854639)
%
%  VARIABLES
%    A = amyloid-beta
%    T = tau
%    I = inflammation
%    D = dead neurons
%
%  EQUATIONS (time is in days)
%    dA/dt = pA + kA1*D - cA*A
%    dT/dt = pT - cT*T
%    dI/dt = kA2*h(A) + kT*h(T) - cI*I
%    dD/dt = pD + kA*h(A) + kTt*h(T) + kI*h(I) - cD*D
%
%  h is a Hill function. From 0 to 1:
%    h(x) = x^q / (K^q + x^q)
%
%% ======================================================================

clear; clc; close all;


%% PARAMETERS USED IN THE MODEL
%  All parameters are stored in one structure called p.

p.pA  = 1.00;   % baseline production of amyloid-beta
p.pT  = 0.10;   % baseline production of tau 
p.pD  = 0.01;   % baseline death of neurons

p.cA  = 1.92;   % clearance rate of amyloid-beta
p.cT  = 0.07;   % clearance rate of tau
p.cI  = 0.70;   % clearance rate of inflammation
p.cD  = 0.10;   % clearance rate of dead neurons

p.kA1 = 1.30;   % amyloid released by each dead neuron
p.kA2 = 1.30;   % how strongly amyloid drives inflammation (same as kA1)
p.kT = 0.05;   % how strongly tau drives inflammation

p.kA = 0.00;   % how strongly amyloid kills neurons (0 = off)
p.kTt = 0.00;  % how strongly tau kills neurons (0 = off)
p.kI = 1.00;   % how strongly inflammation kills neurons

p.KA  = 6.0;    % Hill half-saturation constant for amyloid
p.KT  = 0.5;    % Hill half-saturation constant for tau
p.KI  = 0.5;    % Hill half-saturation constant for inflammation

p.n   = 2;      % Hill coefficient for amyloid and tau
p.m   = 2;      % Hill coefficient for inflammation

% Settings for the simulations
tEnd = 365;     % length of the normal simulation (days)
                %checks that the model behaves as expected

tTip = 2000;    % length of the simulations in the threshold test (days)
                %looks for the tipping point

Imax = 5;       % biggest starting inflammation

% Settings for the sensitivity analysis
relStep = 0.01; % change each parameter by 1 percent

Tstar  = p.pT / p.cT;                                   % steady state of tau


%% EQUILIBRIUM POINTS AND STABILITY OF THE MODEL
%  where all four derivatives are zero

eqs = find_equilibria(p); % finds equalibrium points at baseline parameters -> matrix

if isempty(eqs)
    error('No equilibrium was found.');
end

for j = 1:size(eqs, 1)
    y_eq = eqs(j, :)';            
    J = ad_jac(y_eq, p);          % Jacobian matrix
    ev = eig(J);                  % eigenvalues of the Jacobian
    label = stability_label(ev);  % 'stable', 'saddle' or 'unstable'

    fprintf('\nEquilibrium %d (%s)\n', j, label);
    fprintf('  A = %.4f  T = %.4f  I = %.4f  D = %.4f\n', y_eq);
    fprintf('  Real parts of eigenvalues: ');
    fprintf('%.4f  ', sort(real(ev)));
    fprintf('\n');
end


%% SIMULATION OVER TIME 
%  Start with everything at zero and let the model run.
%   ode15s is used because the model is stiff (some things are fast, some slow).

y0 = [0; 0; 0; 0];                               % start values of A, T, I, D
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);  % accuracy settings

[t, y] = ode15s(@(t, y) ad_rhs(t, y, p), [0 tEnd], y0, opts);

ax.XColor = 'black'; 
ax.YColor = 'black';

fig1 = figure('Name', 'Time series', 'Color', 'black');
names4 = {'Amyloid-\beta A(t)', 'Tau T(t)', 'Inflammation I(t)', 'Dead neurons D(t)'};

for k = 1:4
    subplot(2, 2, k);                            % 2 by 2 grid of plots
    plot(t, y(:, k), 'LineWidth', 2);            % solution of variable k
    hold on;
    % dashed line = the lowest (healthy) equilibrium
    plot([0 tEnd], [eqs(1, k) eqs(1, k)], 'k--');
    xlabel('Time (days)');
    ylabel(names4{k});
    title(names4{k});
    grid on;
end


%% BIFURCATION DIAGRAM AND HYSTERESIS 
%  Blue = stable, red = unstable.

% Bifurication:

kvals = linspace(0, 3, 301);        % kA1 values to test
nEq = zeros(1, length(kvals));      % number of equilibria for each kA1

fig2 = figure('Name', 'Bifurcation diagram', 'Color', 'black');
hold on;

for i = 1:length(kvals)
    q = p;                          
    q.kA1 = kvals(i);               % change kA1
    q.kA2 = kvals(i);               % keep kA2 equal to kA1

    e = find_equilibria(q);         
    nEq(i) = size(e, 1);

    for j = 1:size(e, 1)
        ev = eig(ad_jac(e(j, :)', q));
        if all(real(ev) < 0)
            plot(kvals(i), e(j, 1), 'b.', 'MarkerSize', 6);   % stable
        else
            plot(kvals(i), e(j, 1), 'r.', 'MarkerSize', 6);   % unstable
        end
    end
end

% Hysteresis: 
% start from the steady state of the previous step. If the two paths are different, the system remembers its history.

kC = linspace(0, 3, 81);            
up = zeros(1, length(kC));          % A when kA1 goes up
down = zeros(1, length(kC));        % A when kA1 goes down
yc = [0; Tstar; 0; 0];              % start from the healthy state

for i = 1:length(kC)                % going up
    q = p;
    q.kA1 = kC(i);
    q.kA2 = kC(i);
    yc = settle(yc, q);
    up(i) = yc(1);
end

for i = length(kC):-1:1             % going down
    q = p;
    q.kA1 = kC(i);
    q.kA2 = kC(i);
    yc = settle(yc, q);
    down(i) = yc(1);
end

plot(kC, up, 'ko-', 'MarkerFaceColor', 'g', 'MarkerSize', 4);
plot(kC, down, 'kx-', 'MarkerSize', 6);
xlabel('k_{A1}');
ylabel('Equilibrium amyloid-\beta A^*');
title('Bifurcation diagram in k_{A1}');
grid on;

% Dummy points so the legend only has one entry per curve
h1 = plot(nan, nan, 'b.', 'MarkerSize', 14);
h2 = plot(nan, nan, 'r.', 'MarkerSize', 14);
h3 = plot(nan, nan, 'ko-', 'MarkerFaceColor', 'g');
h4 = plot(nan, nan, 'kx-');
legend([h1 h2 h3 h4], {'stable', 'unstable', 'kA1 increasing', 'kA1 decreasing'}, ...
       'Location', 'northwest');

% The tipping points are where the number of equilibria changes
fprintf('\n=== Tipping points (kA1) ===\n');
changes = find(diff(nEq) ~= 0);     % places where the number changes
if isempty(changes)
    fprintf('The number of equilibria never changes, so no threshold here.\n');
else
    for c = 1:length(changes)
        i = changes(c);
        fprintf('Number of equilibria goes from %d to %d between kA1 = %.3f and %.3f\n', ...
                nEq(i), nEq(i+1), kvals(i), kvals(i+1));
    end
end


%%  MAP OF THE NUMBER OF EQUILIBRIA 
%   1 equilibrium  = only one outcome (monostable)
%   3 equilibria   = two possible outcomes (bistable)

nMap = 40;                              % grid is nMap by nMap
kA1_list = linspace(0.5, 3, nMap);      % kA1 values (x axis)
KA_list = linspace(1, 15, nMap);        % KA values (y axis)
countMap = zeros(nMap, nMap);           % matrix to store the counts

for i = 1:nMap                          % i = row = KA
    for j = 1:nMap                      % j = column = kA1
        q = p;
        q.kA1 = kA1_list(j);
        q.kA2 = kA1_list(j);            % kA2 follows kA1
        q.KA = KA_list(i);
        e = find_equilibria(q);
        countMap(i, j) = size(e, 1);
    end
end

fig3 = figure('Name', 'Equilibria map', 'Color', 'black');
imagesc(kA1_list, KA_list, countMap);   
set(gca, 'YDir', 'normal');             
colorbar;
hold on;
plot(p.kA1, p.KA, 'wp', 'MarkerFaceColor', 'w', 'MarkerSize', 10);   % baseline
xlabel('k_{A1}');
ylabel('K_A');
title('Number of equilibria (white star = baseline)');

fprintf('\n=== Two-parameter map ===\n');
fprintf('Points with 3 equilibria: %d out of %d\n', sum(countMap(:) == 3), numel(countMap));

% 3 equilibrium points

%% FIND THE SADDLE (THE THRESHOLD) 
%  If there are 3 equilibria:
%   lowest: healthy state
%   middle: saddle (the threshold) 
%   highest: chronic (diseased) state.

haveBistable = false;           % becomes true if the system is bistable

eqb = find_equilibria(p);       % equilibria at baseline
fprintf('\n=== Baseline point: %d equilibria ===\n', size(eqb, 1));

labels = cell(size(eqb, 1), 1);
for j = 1:size(eqb, 1)
    ev = eig(ad_jac(eqb(j, :)', p));
    labels{j} = stability_label(ev);
    fprintf('  Eq %d: A=%.4f T=%.4f I=%.4f D=%.4f  %s\n', ...
            j, eqb(j, 1), eqb(j, 2), eqb(j, 3), eqb(j, 4), labels{j});
end

saddleIdx = find(strcmp(labels, 'saddle'));    % which equilibrium is the saddle

if size(eqb, 1) == 3 && length(saddleIdx) == 1
    haveBistable = true;
else
    fprintf('Not bistable here, so Sections 6 and 7 are skipped.\n');
end

if haveBistable
    lowEq = eqb(1, :)';             % healthy state [A; T; I; D]
    sadEq = eqb(saddleIdx, :)';     % saddle = threshold
    highEq = eqb(3, :)';            % chronic state

    fprintf('\nHealthy:  A = %.4f, I = %.4f, D = %.4f\n', lowEq(1), lowEq(3), lowEq(4));
    fprintf('Saddle:   A = %.4f, I = %.4f, D = %.4f\n', sadEq(1), sadEq(3), sadEq(4));
    fprintf('Chronic:  A = %.4f, I = %.4f, D = %.4f\n', highEq(1), highEq(3), highEq(4));
end


%% TEST THE THRESHOLD BY SIMULATION 
%  Test it:
%   start at the healthy state
%   change ONE variable (I or A) to a bigger starting value
%   run for a long time
%   if A ends above the saddle value of A, the system "tipped"
% Bisection is used to find the smallest starting value that tips.

if haveBistable
    A_sad = sadEq(1);
    I_sad = sadEq(3);

    I0_crit = critical_start(p, lowEq, A_sad, 3, Imax, tTip);
    A0_crit = critical_start(p, lowEq, A_sad, 1, 1.5*highEq(1), tTip);

    fprintf('\n=== Thresholds found by simulation ===\n');
    fprintf('I of the saddle (guess):        %.4f\n', I_sad);
    if isnan(I0_crit)
        fprintf('Critical I(0): not found.\n');
    else
        fprintf('Critical I(0) from simulation:  %.4f\n', I0_crit);
    end
    if isnan(A0_crit)
        fprintf('Critical A(0): not found.\n');
    else
        fprintf('Critical A(0) from simulation:  %.4f (A of saddle = %.4f)\n', A0_crit, A_sad);
    end

    % ----- Plot 1: different starting inflammation -----
    if ~isnan(I0_crit)
        I0list = [1.1*I_sad, 0.9*I0_crit, 1.1*I0_crit];   % three start values

        fig4 = figure('Name', 'Threshold test in I', 'Color', 'w');
        for r = 1:length(I0list)
            start = lowEq;              % begin at the healthy state
            start(3) = I0list(r);       % but change I
            [tt, yy] = ode15s(@(t, y) ad_rhs(t, y, p), [0 tTip], start, opts);

            subplot(2, 1, 1); hold on;
            plot(tt, yy(:, 3), 'LineWidth', 1.6, 'DisplayName', sprintf('I(0) = %.3f', I0list(r)));
            subplot(2, 1, 2); hold on;
            plot(tt, yy(:, 1), 'LineWidth', 1.6, 'DisplayName', sprintf('I(0) = %.3f', I0list(r)));
        end
        subplot(2, 1, 1);
        xlabel('Time (days)'); ylabel('I'); grid on; legend('Location', 'northeast');
        title('Healthy start with extra inflammation');
        subplot(2, 1, 2);
        plot([0 tTip], [A_sad A_sad], 'k--', 'DisplayName', 'A of saddle');
        xlabel('Time (days)'); ylabel('A'); grid on;
    end

    % ----- Plot 2: different starting amyloid -----
    if isnan(A0_crit)
        A0list = A_sad * [0.5 1 2 4];       
    else
        A0list = A0_crit * [0.5 0.9 1.1 1.5];
    end

    fig5 = figure('Name', 'Tipping', 'Color', 'black');
    hold on;
    for r = 1:length(A0list)
        start = lowEq;
        start(1) = A0list(r);               
        [tt, yy] = ode15s(@(t, y) ad_rhs(t, y, p), [0 tTip], start, opts);
        plot(tt, yy(:, 1), 'LineWidth', 1.6, 'DisplayName', sprintf('A(0) = %.3f', A0list(r)));
    end
    plot([0 tTip], [A_sad A_sad], 'k--', 'DisplayName', 'A of saddle');
    xlabel('Time (days)');
    ylabel('Amyloid-\beta A(t)');
    title('Different starting values of A');
    legend('Location', 'east');
    grid on;
end


%% SENSITIVITY ANALYSIS -----------------------------------------------
% How much does the healthy A* change if one parameter changes a little?
%    S = (relative change in A*) / (relative change in the parameter)
%
% We change each parameter by +1% and -1%. n and m are left out because
% they are exponents. Parameters equal to 0 give S = 0.

parNames = {'pA','pT','pD','cA','cT','cI','cD','kA1','kA2','kT', ...
            'kA','kTt','kI','KA','KT','KI'};
nPar = length(parNames);
S = zeros(nPar, 1);             % sensitivities
A_base = eqs(1, 1);             % healthy A at baseline

for i = 1:nPar
    name = parNames{i};

    qp = p;
    qp.(name) = p.(name) * (1 + relStep);   % parameter 1% bigger
    qm = p;
    qm.(name) = p.(name) * (1 - relStep);   % parameter 1% smaller

    ep = find_equilibria(qp);
    em = find_equilibria(qm);

    if isempty(ep) || isempty(em)
        S(i) = 0;                           % nothing found
    else
        % take the equilibrium closest to the baseline healthy A
        [~, ip] = min(abs(ep(:, 1) - A_base));
        [~, im] = min(abs(em(:, 1) - A_base));
        A_plus = ep(ip, 1);
        A_minus = em(im, 1);

        S(i) = ((A_plus - A_minus) / A_base) / (2 * relStep);
    end
end

[~, order] = sort(abs(S), 'descend');       % biggest effect first

fig6 = figure('Name', 'Sensitivity', 'Color', 'black');
barh(S(order));
set(gca, 'YTick', 1:nPar, 'YTickLabel', parNames(order), 'YDir', 'reverse');
xlabel('Sensitivity of A^*');
title('Local sensitivity analysis');
grid on;


%% ======================================================================
%  FUNCTIONS 
%% ======================================================================

function h = hill(x, K, q)      % Hill function.
    h = x.^q ./ (K^q + x.^q);   % x ranges between 0 and 1
end
                  
function dh = dhill(x, K, q)   % Derivative of the Hill function
    dh = q * K^q * x.^(q-1) ./ (K^q + x.^q).^2;
end

function dy = ad_rhs(~, y, p)  % y is a column [A; T; I; D], time isn't used.
    A = y(1);
    T = y(2);
    I = y(3);
    D = y(4);

    dA = p.pA + p.kA1*D - p.cA*A;
    dT = p.pT - p.cT*T;
    dI = p.kA2*hill(A, p.KA, p.n) + p.kT*hill(T, p.KT, p.n) - p.cI*I;
    dD = p.pD + p.kA*hill(A, p.KA, p.n) + p.kTt*hill(T, p.KT, p.n) ...
         + p.kI*hill(I, p.KI, p.m) - p.cD*D;

    dy = [dA; dT; dI; dD];
end


function J = ad_jac(y, p) % Jacobian
                          % J(i,j) is how much equation i changes when variable j changes.
    A = y(1);
    T = y(2);
    I = y(3);

    J = zeros(4, 4);        

    % Row 1: dA/dt
    J(1, 1) = -p.cA;
    J(1, 4) = p.kA1;

    % Row 2: dT/dt
    J(2, 2) = -p.cT;

    % Row 3: dI/dt
    J(3, 1) = p.kA2 * dhill(A, p.KA, p.n);
    J(3, 2) = p.kT * dhill(T, p.KT, p.n);
    J(3, 3) = -p.cI;

    % Row 4: dD/dt
    J(4, 1) = p.kA * dhill(A, p.KA, p.n);
    J(4, 2) = p.kTt * dhill(T, p.KT, p.n);
    J(4, 3) = p.kI * dhill(I, p.KI, p.m);
    J(4, 4) = -p.cD;
end


function label = stability_label(ev) % stability labelling, based on eigenvalues

    nPositive = sum(real(ev) > 0);
    if nPositive == 0
        label = 'stable';
    elseif nPositive == 1
        label = 'saddle';
    else
        label = 'unstable';
    end
end


function [I, D] = get_I_D(A, p)
    % At equilibrium I and D can be written using A only.
    % (T is always pT/cT at equilibrium.)
    Tst = p.pT / p.cT;
    I = (p.kA2*hill(A, p.KA, p.n) + p.kT*hill(Tst, p.KT, p.n)) / p.cI;
    D = (p.pD + p.kA*hill(A, p.KA, p.n) + p.kTt*hill(Tst, p.KT, p.n) ...
         + p.kI*hill(I, p.KI, p.m)) / p.cD;
end


function g = reduced_g(A, p)
    % Puts I and D into the equation for dA/dt, which leaves ONE equation
    % in A. The equilibria are the values of A where g(A) = 0.
    [~, D] = get_I_D(A, p);
    g = p.pA + p.kA1*D - p.cA*A;
end


function eqs = find_equilibria(p)
    % Finds all equilibria. 
    % Each row of the output is [A T I D].

    Tst = p.pT / p.cT;

    % A can never be bigger than Amax, so we only look between 0 and Amax
    Dmax = (p.pD + p.kA + p.kTt + p.kI) / p.cD;
    Amax = (p.pA + p.kA1*Dmax) / p.cA;

    Agrid = linspace(0, 1.01*Amax, 20001);   
    g = reduced_g(Agrid, p);                 % g at every A value

    idx = find(g(1:end-1) .* g(2:end) < 0); % where does g, undergo a sign change

    eqs = zeros(length(idx), 4);                       % find roots in between sign change locations
    for j = 1:length(idx)
        a = Agrid(idx(j));                   
        b = Agrid(idx(j) + 1);               
        A_root = fzero(@(A) reduced_g(A, p), [a b]);
        [I, D] = get_I_D(A_root, p);
        eqs(j, :) = [A_root, Tst, I, D];
    end
end


function yEnd = settle(y0, p)
    % Runs the model until it stops changing much.
    % Runs for 1000 days at a time, for up to 20 times.
    o = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    y = y0;
    for piece = 1:20
        [~, yy] = ode15s(@(t, y) ad_rhs(t, y, p), [0 1000], y, o);
        y = yy(end, :)';                          % Save the last state
        if max(abs(ad_rhs(0, y, p))) < 1e-6       
            break;
        end
    end
    yEnd = y;
end


function chronic = ends_chronic(p, y0, A_sad, tTip)
    % Runs the model for tTip days.
    % Checks if the final amyloid value is above the saddle value.
    o = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    [~, yy] = ode15s(@(t, y) ad_rhs(t, y, p), [0 tTip], y0, o);
    chronic = yy(end, 1) > A_sad;
end


function vcrit = critical_start(p, lowEq, A_sad, idx, vmax, tTip)
    % Finds the starting value that causes the system to tip.
    % idx = 1 for amyloid and idx = 3 for inflammation.
    % The other variables start at their healthy equilibrium values.
    % Gives NaN if even vmax does not tip.

    ytest = lowEq;
    ytest(idx) = vmax;
    if ~ends_chronic(p, ytest, A_sad, tTip)
        vcrit = NaN;                 % system doesn't tip at any point
        return;
    end

    v_low = lowEq(idx);              % this value does not tip
    v_high = vmax;                   % this value tips
    for it = 1:30
        mid = (v_low + v_high) / 2;  % try the middle
        ytest = lowEq;
        ytest(idx) = mid;            % Tip occurs, so lower the upper value
        if ends_chronic(p, ytest, A_sad, tTip)
            v_high = mid;            % No tip, so raise the lower value
        else
            v_low = mid;             % does not tip, so the threshold is higher
        end
    end
    vcrit = (v_low + v_high) / 2;
end