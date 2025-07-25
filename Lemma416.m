longprecision(32);  % Use high-precision interval arithmetic

% Step 1: Initial setup
MM(1) = intval('0.0093');      % Start with rigorous interval
M(1) = inf(MM(1));             % Use lower bound
MM(2) = LSigmaRR(M(1), -M(1));  % First function evaluation
M(2) = inf(MM(2));             % Lower bound of result

% Step 2: Loop
j = 1;
while true
    j = j + 1;
    MM(j+1) = LSigmaR(M(j), -M(j));
    M(j+1) = inf(MM(j+1));
    
    % Check stopping condition *after* defining M(j+1)
    if ~(M(j) < M(j+1) && M(j) < 0.22)
        break
    end
end

% Step 3: Display results
MM(j)
disp(['j     = ', num2str(j)])
disp(['M(j)  = ', num2str(M(j))])
disp(['M(j-1)= ', num2str(M(j-1))])

%intval ans = 0.26269143925511, j     = 97, M(j)  = 0.2627, M(j-1)= 0.20989