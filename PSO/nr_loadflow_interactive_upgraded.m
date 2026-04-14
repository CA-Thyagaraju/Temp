function nr_loadflow_interactive()

Ybus = input('Ybus = ');
n = size(Ybus,1);

S = input('Apparent Power: ');

slack = input('Slack bus number: ');
Vslack = input('Slack bus voltage: ');
S(slack) = 0;

tol = 0.0001;
max_iter = 100;

% ---------- INITIALIZATION ----------
V = ones(n,1);
V(slack) = Vslack;

% PV bus voltages (Nx2 matrix input)
pv_input = input('PV voltage matrix [bus_number  complex_voltage]: ');

if isempty(pv_input)
    pv = [];
else
    pv = pv_input(:,1).';   % extract PV bus numbers
    
    for r = 1:size(pv_input,1)
        bus = pv_input(r,1);
        V(bus) = pv_input(r,2);
    end
end

% ---------- NR ITERATION ----------
for iter = 1:max_iter

    P = zeros(n,1);
    Q = zeros(n,1);

    for i = 1:n
        for k = 1:n
            G = real(Ybus(i,k));
            B = imag(Ybus(i,k));
            P(i) = P(i) + abs(V(i))*abs(V(k)) * (G*cos(angle(V(i))-angle(V(k))) + B*sin(angle(V(i))-angle(V(k))));
            Q(i) = Q(i) + abs(V(i))*abs(V(k)) * (G*sin(angle(V(i))-angle(V(k))) - B*cos(angle(V(i))-angle(V(k))));
        end
    end

    non_slack = setdiff(1:n, slack);
    pq = setdiff(non_slack, pv);

    dP = real(S(non_slack)) - P(non_slack);
    dQ = imag(S(pq)) - Q(pq);

    if max([abs(dP); abs(dQ)]) < tol
        break;
    end

    ntheta = length(non_slack);
    npq = length(pq);

    J1 = zeros(ntheta,ntheta);
    J2 = zeros(ntheta,npq);
    J3 = zeros(npq,ntheta);
    J4 = zeros(npq,npq);

    for ii = 1:ntheta
        i = non_slack(ii);
        for kk = 1:ntheta
            k = non_slack(kk);
            G = real(Ybus(i,k));
            B = imag(Ybus(i,k));
            ang = angle(V(i)) - angle(V(k));
            if i == k
                J1(ii,kk) = -Q(i) - imag(Ybus(i,i))*abs(V(i))^2;
            else
                J1(ii,kk) = abs(V(i))*abs(V(k)) * (G*sin(ang) - B*cos(ang));
            end
        end
    end

    for ii = 1:ntheta
        i = non_slack(ii);
        for kk = 1:npq
            k = pq(kk);
            G = real(Ybus(i,k));
            B = imag(Ybus(i,k));
            ang = angle(V(i)) - angle(V(k));
            if i == k
                J2(ii,kk) = P(i) + real(Ybus(i,i))*abs(V(i))^2;
            else
                J2(ii,kk) = abs(V(i)) * (G*cos(ang) + B*sin(ang));
            end
        end
    end

    for ii = 1:npq
        i = pq(ii);
        for kk = 1:ntheta
            k = non_slack(kk);
            G = real(Ybus(i,k));
            B = imag(Ybus(i,k));
            ang = angle(V(i)) - angle(V(k));
            if i == k
                J3(ii,kk) = P(i) - real(Ybus(i,i))*abs(V(i))^2;
            else
                J3(ii,kk) = -abs(V(i))*abs(V(k)) * (G*cos(ang) + B*sin(ang));
            end
        end
    end

    for ii = 1:npq
        i = pq(ii);
        for kk = 1:npq
            k = pq(kk);
            G = real(Ybus(i,k));
            B = imag(Ybus(i,k));
            ang = angle(V(i)) - angle(V(k));
            if i == k
                J4(ii,kk) = Q(i) - imag(Ybus(i,i))*abs(V(i))^2;
            else
                J4(ii,kk) = abs(V(i)) * (G*sin(ang) - B*cos(ang));
            end
        end
    end

    J = [J1 J2; J3 J4];
    dx = J \ [dP; dQ];

    dtheta = dx(1:ntheta);
    dV = dx(ntheta+1:end);

    for ii = 1:ntheta
        i = non_slack(ii);
        V(i) = abs(V(i)) * exp(1j*(angle(V(i)) + dtheta(ii)));
    end

    for ii = 1:npq
        i = pq(ii);
        V(i) = abs(V(i))*(1 + dV(ii)) * exp(1j*angle(V(i)));
    end
end

fprintf('Converged in %d iterations\n\n', iter);
for i = 1:n
    fprintf('Bus %d: |V| = %.4f pu , Angle = %.2f deg\n', i, abs(V(i)), rad2deg(angle(V(i))));
end

end
