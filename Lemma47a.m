% Initialize
longprecision(32);             % Set INTLAB to high precision
a = intval('-1.5');            % Replace with appropriate value for `a`
epsval = intval('1e-16');      % Small positive number (for robustness)

M(1) = intval('0.0093') - epsval;
j = 1;

z(1) = 1 + 0.25 * a^2 / LDaR(M(1), -M(1), a);
M(2) = 1 / (z(1) + sqrt(z(1)^2 - 1)) - epsval;

while (sup(M(j)) < inf(M(j+1))) && (sup(M(j+1)) < 0.26)
    j = j + 1;
    z(j) = 1 + 0.25 * a^2 / LDaR(M(j), -M(j), a);
    M(j+1) = 1 / (z(j) + sqrt(z(j)^2 - 1)) - epsval;
end

% Display results: intval ans = 0.23434557682237, j = 62, M(j) = 0.23435, M(j+1) = 0.23435
M(j)
disp(['j = ', num2str(j)])
disp(['M(j) = ',  num2str(sup(M(j)))])
disp(['M(j+1) = ', num2str(sup(M(j+1)))])

%intval ans = 0.23434557682237, j = 62, M(j) = 0.23435, M(j+1) = 0.23435