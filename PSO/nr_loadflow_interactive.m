function nr_loadflow_interactive()

% ---------- INPUTS ----------
Ybus = input('Ybus = ');
n = size(Ybus,1);              % number of buses

S = input('Apparent Power: ');
S = -S;                        % load convention (P - jQ)

Vslack = input('Slack bus voltage: ');

tol = 0.0001;
max_iter = 100;

% ---------- INITIALIZATION ----------
V = ones(n,1);                 % initial guess
V(1) = Vslack;                 % slack bus fixed

% ---------- NEWTON–RAPHSON ----------
for iter = 1:max_iter

    V_old = V;

    % ---- Calculate P and Q ----
    P = zeros(n,1);
    Q = zeros(n,1);

    for i = 1:n
        for k = 1:n
            G = real(Ybus(i,k));
            B = imag(Ybus(i,k));
            P(i) = P(i) + abs(V(i))*abs(V(k)) * ...
                   (G*cos(angle(V(i))-angle(V(k))) + ...
                    B*sin(angle(V(i))-angle(V(k))));
            Q(i) = Q(i) + abs(V(i))*abs(V(k)) * ...
                   (G*sin(angle(V(i))-angle(V(k))) - ...
                    B*cos(angle(V(i))-angle(V(k))));
        end
    end

    % ---- Mismatch (skip slack) ----
    dP = real(S(2:n)) - P(2:n);
    dQ = imag(S(2:n)) - Q(2:n);

    if max([abs(dP); abs(dQ)]) < tol
        break;
    end

    % ---- Jacobian ----
    npq = n - 1;

    J1 = zeros(npq,npq);
    J2 = zeros(npq,npq);
    J3 = zeros(npq,npq);
    J4 = zeros(npq,npq);

    for i = 2:n
        for k = 2:n
            ii = i-1; kk = k-1;
            G = real(Ybus(i,k));
            B = imag(Ybus(i,k));
            ang = angle(V(i)) - angle(V(k));

            if i == k
                J1(ii,kk) = -Q(i) - imag(Ybus(i,i))*abs(V(i))^2;
                J2(ii,kk) =  P(i)/abs(V(i)) + real(Ybus(i,i))*abs(V(i));
                J3(ii,kk) =  P(i) - real(Ybus(i,i))*abs(V(i))^2;
                J4(ii,kk) =  Q(i)/abs(V(i)) - imag(Ybus(i,i))*abs(V(i));
            else
                J1(ii,kk) = abs(V(i))*abs(V(k)) * ...
                             (G*sin(ang) - B*cos(ang));
                J2(ii,kk) = abs(V(i)) * ...
                            (G*cos(ang) + B*sin(ang));
                J3(ii,kk) = -abs(V(i))*abs(V(k)) * ...
                            (G*cos(ang) + B*sin(ang));
                J4(ii,kk) = abs(V(i)) * ...
                            (G*sin(ang) - B*cos(ang));
            end
        end
    end

    % ---- Solve corrections ----
    J = [J1 J2; J3 J4];
    dx = J \ [dP; dQ];

    dtheta = dx(1:npq);
    dV = dx(npq+1:end);

    % ---- Update voltages ----
    for i = 2:n
        V(i) = (abs(V(i)) + dV(i-1)) * ...
               exp(1j*(angle(V(i)) + dtheta(i-1)));
    end

    if max(abs(V - V_old)) < tol
        break;
    end
end

% ---------- RESULTS ----------
fprintf('Converged in %d iterations\n\n', iter);
for i = 1:n
    fprintf('Bus %d: |V| = %.4f pu , Angle = %.2f deg\n', ...
        i, abs(V(i)), rad2deg(angle(V(i))));
end

end
