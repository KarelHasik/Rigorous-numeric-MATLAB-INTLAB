function l = theta6R(m)
% Rigorous upper bound for theta_6(m)
longprecision(32);
a = -intval(37)/intval(24);

Am = a * m ./ (1 + m) - 1 + (1 + m) .* log(1 + m) ./ m;  % interval
l = Am;   % initialize with  Am

for j = 1:6
    l = USigmaRR(l, m);  % rigorous update with interval input
end
end
