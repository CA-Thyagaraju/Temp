function gs_loadflow_interactive()

Ybus = input('Ybus = ');
n = size(Ybus,1);

S = input('Apparant Power: ');
S = -(S);

Vslack = input('Slack bus voltage: ');

alpha = input('Acceleration factor: ');

tol = 0.0001;
max_iter = 100;

V = ones(n,1);
V(1) = Vslack;

for iter = 1:max_iter
    V_old = V;
    
    for i = 2:n
        sumYV = 0;
        for k = 1:n
            if k ~= i
                sumYV = sumYV + Ybus(i,k)*V(k);
            end
        end
        
        V_new = (1/Ybus(i,i)) * (conj(S(i))/conj(V(i)) - sumYV);
        V(i) = V(i) + alpha*(V_new - V(i));
    end
    
    if max(abs(V - V_old)) < tol
        break;
    end
end

fprintf('Converged in %d iterations\n\n', iter);

for i = 1:n
    fprintf('Bus %d: |V| = %.4f pu , Angle = %.2f deg\n', i, abs(V(i)), rad2deg(angle(V(i))));
end

% Slack bus power
sumI = 0;
for k = 1:n
    sumI = sumI + Ybus(1,k)*V(k);
end
S_slack = V(1)*conj(sumI);

fprintf('\nSlack Bus Power:\n');
fprintf('P = %.4f pu , Q = %.4f pu\n', real(S_slack), imag(S_slack));

% Line flows
y = input('\nEnter line admittance matrix y (same size as Ybus): ');

fprintf('\nLine Flows:\n');
for i = 1:n
    for k = i+1:n
        if y(i,k) ~= 0
            Sik = V(i)*conj((V(i)-V(k))*y(i,k));
            Ski = V(k)*conj((V(k)-V(i))*y(i,k));
            fprintf('S%d%d = %.4f + j%.4f\n', i,k, real(Sik), imag(Sik));
            fprintf('S%d%d = %.4f + j%.4f\n\n', k,i, real(Ski), imag(Ski));
        end
    end
end

end
