function Ybus = formYbus(y)
n = size(y,1);          
Ybus = zeros(n,n);

for i = 1:n
    for j = 1:n
        if i == j
            % Diagonal element: self + sum of connected line admittances
            Ybus(i,i) = y(i,i) + sum(y(i,[1:i-1, i+1:n]));
        else
            % Off-diagonal element
            Ybus(i,j) = -y(i,j);
        end
    end
end
end
