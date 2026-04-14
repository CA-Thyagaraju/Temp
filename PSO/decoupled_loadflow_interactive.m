
function decoupled_loadflow_interactive()
% ---------- INPUT ----------
y = input('Enter line admittance matrix y: ');

NB = size(y,1);        % number of buses

Ybus = zeros(NB, NB);

% ---------- FORM YBUS ----------
for i = 1:NB
    for j = 1:NB
        
        if i == j
            % Diagonal element
            Ybus(i,i) = sum(y(i,:));
        else
            % Off-diagonal element
            Ybus(i,j) = -y(i,j);
        end
        
    end
end

% ---------- POWER DATA ----------
S = input('Enter power vector S (Pg-Pl + j(Qg-Ql)): ');

% ---------- SLACK BUS ----------
slack = input('Slack bus number: ');
Vslack = input('Slack bus voltage: ');

% ---------- PV BUS DATA ----------
pv_input = input('Enter PV bus data as [bus_number  voltage]: ');

if isempty(pv_input)
    pv = [];
else
    pv = pv_input(:,1).';
end

% ---------- INITIALIZATION ----------
V = ones(NB,1);
delta = zeros(NB,1);

% Slack bus voltage
V(slack) = Vslack;

% PV bus voltage magnitudes
for r = 1:size(pv_input,1)
    V(pv_input(r,1)) = pv_input(r,2);
end

% Separate G and B matrices
G = real(Ybus);
B = imag(Ybus);

%TOLERANCE AND ITERATIONS SETTING
tol = 1e-4;
max_iter = 50;

% ---------- ITERATION ----------
for iter = 1:max_iter

    P = zeros(NB,1);
    Q = zeros(NB,1);

    % ---- Compute Pi and Qi ----
    for i = 1:NB
        for k = 1:NB
            P(i) = P(i) + V(i)*V(k)*(G(i,k)*cos(delta(i)-delta(k)) + B(i,k)*sin(delta(i)-delta(k)));

            Q(i) = Q(i) + V(i)*V(k)*(G(i,k)*sin(delta(i)-delta(k)) - B(i,k)*cos(delta(i)-delta(k)));
        end
    end

    % ---- Active power mismatch ----
    non_slack = setdiff(1:NB, slack);
    dP = real(S(non_slack)) - P(non_slack);


% ---------- FORM H MATRIX ----------
ntheta = length(non_slack);
H = zeros(ntheta);

for ii = 1:ntheta
    i = non_slack(ii);

    for kk = 1:ntheta
        k = non_slack(kk);

        if i == k
            H(ii,kk) = -Q(i) - B(i,i)*V(i)^2;
        else
            H(ii,kk) = V(i)*V(k)*(G(i,k)*sin(delta(i)-delta(k)) - B(i,k)*cos(delta(i)-delta(k)));
        end
    end
end

% ---------- SOLVE FOR Δδ ----------
dtheta = H \ dP;

% ---------- UPDATE ANGLES ----------
delta(non_slack) = delta(non_slack) + dtheta;

% ---------- RECOMPUTE Q AFTER δ UPDATE ----------
Q = zeros(NB,1);
for i = 1:NB
    for k = 1:NB
        Q(i) = Q(i) + V(i)*V(k)*(G(i,k)*sin(delta(i)-delta(k)) ...
               - B(i,k)*cos(delta(i)-delta(k)));
    end
end


% ---------- PQ BUS SET ----------
pq = setdiff(non_slack, pv);

% ---------- REACTIVE POWER MISMATCH ----------
dQ = imag(S(pq)) - Q(pq);

% ---------- CONVERGENCE CHECK ----------
if max(abs([dP; dQ])) < tol
    break;
end

% ---------- FORM L MATRIX ----------
if isempty(pq)
    continue;
end
npq = length(pq);
L = zeros(npq);

for ii = 1:npq
    i = pq(ii);

    for kk = 1:npq
        k = pq(kk);

        if i == k
            L(ii,kk) = Q(i) - B(i,i)*V(i)^2;
        else
            L(ii,kk) = V(i)*V(k)*(G(i,k)*sin(delta(i)-delta(k)) - B(i,k)*cos(delta(i)-delta(k)));
        end
    end
end

% ---------- SOLVE FOR ΔV ----------
dV = L \ dQ;

% ---------- UPDATE VOLTAGES ----------
for ii = 1:npq
    i = pq(ii);
    V(i) = V(i) + dV(ii);
end

end
% ---------- RESULTS ----------
fprintf('\nConverged in %d iterations\n\n', iter);

for i = 1:NB
    fprintf('Bus %d: |V| = %.4f pu , Angle = %.2f deg\n', i, V(i), rad2deg(delta(i)));
end
