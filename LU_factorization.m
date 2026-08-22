clear; clc;
mpc = case30;

%% 1. Load and Parse MATPOWER Data
baseMVA = mpc.baseMVA;
nbus = size(mpc.bus,1);
nbranch = size(mpc.branch,1);

bus_id = mpc.bus(:,1);            
bus_type = mpc.bus(:,2);         % 1=PQ(24), 2=PV(5), 3=Slack(1) for case_30
Pd = mpc.bus(:,3) / baseMVA;     %real power demand (MW)
Qd = mpc.bus(:,4) / baseMVA;     %reactive power demand (MVAr)
Gs = mpc.bus(:,5) / baseMVA;     %shunt conductance (MW demanded at V = 1.0 p.u.)
Bs = mpc.bus(:,6) / baseMVA;     %shunt susceptance (MVAr injected at V = 1.0 p.u.)
V = mpc.bus(:,8);                % Initial Voltage Magnitude Guess
delta = mpc.bus(:,9) * pi / 180; %Initial Voltage angle converted to radians

%Calculate Pg, Qg using generator data
Pg = zeros(nbus,1);
Qg = zeros(nbus,1);
for i = 1:size(mpc.gen,1)
    idx = find(bus_id == mpc.gen(i,1));
    Pg(idx) = Pg(idx) + mpc.gen(i,2) / baseMVA;
    Qg(idx) = Qg(idx) + mpc.gen(i,3) / baseMVA;
    %update voltage of PV/slack bus
    if (bus_type(idx) == 2 || bus_type(idx) == 3)
        V(idx) = mpc.gen(i,6);
    end
end

P_spec = Pg - Pd;
Q_spec = Qg - Qd;

pq = find(bus_type == 1);
pv = find(bus_type == 2);
slack = find(bus_type == 3);
non_slack = [pv; pq];

num_non_slack = length(non_slack);
num_pq = length(pq);

%% 2. Formulate Ybus matrix
Ybus = zeros(nbus, nbus);

for k = 1:nbranch
    from = find(bus_id == mpc.branch(k, 1));
    to = find(bus_id == mpc.branch(k, 2));
    r = mpc.branch(k, 3);
    x = mpc.branch(k, 4);
    b_shunt = mpc.branch(k, 5);
    ratio = mpc.branch(k, 9);
    
    if ratio == 0
        ratio = 1.0;
    end
    
    y = 1 / (r + 1i * x);
    
    Ybus(from, to) = Ybus(from, to) - y / ratio;
    Ybus(to, from) = Ybus(to, from) - y / ratio;   
    Ybus(from, from) = Ybus(from, from) + (y + 1i * b_shunt / 2) / (ratio^2);
    Ybus(to, to) = Ybus(to, to) + y + 1i * b_shunt / 2;
end

for i = 1:nbus
    Ybus(i, i) = Ybus(i, i) + Gs(i) + 1i * Bs(i);
end

G = real(Ybus);
B = imag(Ybus);

%% 3. Newton Raphson iterations
tol = 1e-4;
max_iter = 100;
iter = 0;
converged = false;

while iter < max_iter
    iter = iter + 1;
    
    P_calc = zeros(nbus, 1);
    Q_calc = zeros(nbus, 1);
    
    for i = 1:nbus
        for j = 1:nbus
            theta_ij = delta(i) - delta(j);
            P_calc(i) = P_calc(i) + V(i)*V(j)*(G(i,j)*cos(theta_ij) + B(i,j)*sin(theta_ij));
            Q_calc(i) = Q_calc(i) + V(i)*V(j)*(G(i,j)*sin(theta_ij) - B(i,j)*cos(theta_ij));
        end
    end
 
    dP = P_spec(non_slack) - P_calc(non_slack);
    dQ = Q_spec(pq) - Q_calc(pq);
    mismatch = [dP; dQ];
    
    fprintf('Iteration %d: Max Mismatch = %f\n', iter, max(abs(mismatch)));
    
    if max(abs(mismatch)) < tol
        converged = true;
        break;
    end
    
    J11 = zeros(num_non_slack, num_non_slack);
    J12 = zeros(num_non_slack, num_pq);
    J21 = zeros(num_pq, num_non_slack);
    J22 = zeros(num_pq, num_pq);
    
    for i = 1:num_non_slack
        n = non_slack(i);
        
        % J11: dP / dDelta
        for j = 1:num_non_slack
            m = non_slack(j);
            if n == m
                J11(i,j) = -Q_calc(n) - (V(n)^2)*B(n,n);
            else
                J11(i,j) = V(n)*V(m)*(G(n,m)*sin(delta(n)-delta(m)) - B(n,m)*cos(delta(n)-delta(m)));
            end
        end
        
        % J12: dP / dV
        for j = 1:num_pq
            m = pq(j);
            if n == m
                J12(i,j) = (P_calc(n) / V(n)) + G(n,n)*V(n);
            else
                J12(i,j) = V(n)*(G(n,m)*cos(delta(n)-delta(m)) + B(n,m)*sin(delta(n)-delta(m)));
            end
        end
    end
    
    for i = 1:num_pq
        n = pq(i);
        
        % J21: dQ / dDelta
        for j = 1:num_non_slack
            m = non_slack(j);
            if n == m
                J21(i,j) = P_calc(n) - (V(n)^2)*G(n,n);
            else
                J21(i,j) = -V(n)*V(m)*(G(n,m)*cos(delta(n)-delta(m)) + B(n,m)*sin(delta(n)-delta(m)));
            end
        end
        
        % J22: dQ / dV
        for j = 1:num_pq
            m = pq(j);
            if n == m
                J22(i,j) = (Q_calc(n) / V(n)) - B(n,n)*V(n);
            else
                J22(i,j) = V(n)*(G(n,m)*sin(delta(n)-delta(m)) - B(n,m)*cos(delta(n)-delta(m)));
            end
        end
    end
    
    J = [J11, J12; J21, J22];
    J(abs(J) < 1e-10) = 0;

    % Solve using LU factorization
    [dx, L, U] = calculate_lu(J, mismatch);
    
    % Update State Variables
    delta(non_slack) = delta(non_slack) + dx(1:num_non_slack);
    V(pq) = V(pq) + dx(num_non_slack+1 : end);
end


%% 4. Display Results
if converged
    fprintf('\nSystem Converged in %d iterations\n', iter);
    fprintf('Bus\t\t V (pu)\t\t Angle (deg)\n');
    for i = 1:nbus
        fprintf('%d\t\t %.4f\t\t %.4f\n', bus_id(i), V(i), delta(i)*180/pi);
    end
else
    fprintf('\nFailed to converge after %d iterations.\n', max_iter);
end

%% 5. LU factorization function
function [x, L, U] = calculate_lu(J, mismatch)
    [n, ~] = size(J);
    L = eye(n); 
    U = zeros(n);

    % LU Decomposition
    for i = 1:n
        % Calculate the row of U
        for j = i:n
            sum_U = 0;
            for k = 1:i-1
                sum_U = sum_U + L(i,k) * U(k,j);
            end
            U(i, j) = J(i, j) - sum_U;
        end
        
        % Calculate the column of L
        for j = i+1:n
            sum_L = 0;
            for k = 1:i-1
                sum_L = sum_L + L(j, k) * U(k, i);
            end
            L(j, i) = (J(j, i) - sum_L) / U(i, i);
        end
    end
    
    % Forward Substitution (Ly = b)
    y = zeros(n, 1);
    for i = 1:n
        sum_y = 0;
        for j = 1:i-1
            sum_y = sum_y + L(i, j) * y(j);
        end
        y(i) = mismatch(i) - sum_y;
    end
    
    % Backward Substitution (Ux = y)
    x = zeros(n, 1);
    for i = n:-1:1
        sum_x = 0;
        for j = i+1:n
            sum_x = sum_x + U(i, j) * x(j);
        end
        x(i) = (y(i) - sum_x) / U(i, i);
    end
end
