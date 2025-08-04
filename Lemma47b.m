% Initialization
longprecision(32);                 % Set INTLAB to high-precision mode
a = -intval(3)./intval(2);                % a = -3/2 as an INTLAB interval
epsval = intval('1e-16');           % Use eps consistent with interval arithmetic

M(1) = intval('0.21') - epsval;    % Initial value for M(1)
j = 1;

% Compute z(1) and M(2) using interval-safe expressions
z(1) = 1 + 0.25 * a^2 / LDaR(intval('0.27'), -M(1), a);
M(2) = 1 / (z(1) + sqrt(z(1)^2 - 1)) - epsval;

% Loop using rigorous bounds for interval comparisons
while (sup(M(j)) < inf(M(j+1))) && (sup(M(j+1)) < 0.26)
    j = j + 1;
    z(j) = 1 + 0.25 * a^2 / LDaR(intval('0.27'), -M(j), a);
    M(j+1) = 1 / (z(j) + sqrt(z(j)^2 - 1)) - epsval;
end

% Display final results j = 12, M(j) = 0.25155, M(j+1) = 0.25155
disp(['j = ', num2str(j)])
disp(['M(j) = ', num2str(sup(M(j)))])
disp(['M(j+1) = ', num2str(sup(M(j+1)))])

%j = 12, M(j) = 0.25156, M(j+1) = 0.25156