tic
% Define subinterval [-0.0385, -0.0285] partitioned into 1000 subintervals
n = 1000;                       % number of subdivisions in m
h = intval('0.00002');        % step size in M-direction
delta = 0.00001;             % width of each m-subinterval

W = cell(1,n);                % use cell array to handle variable N(j)

for j = 1:n
    % Define interval endpoints as infsup for rigor
    m_left = -0.0285 - j * delta;
    m_right = -0.0285 - (j-1) * delta;
    m2_j = infsup(m_left, m_left);    % lower endpoint as interval
    m1_j = infsup(m_right, m_right);  % upper endpoint as interval
    
    mY_j = hull(m1_j, m2_j);          % rigorously enclosing interval
    a2_j = -m1_j;                     % interval a2(j)
    
    % Use formula with interval safety
    A_j = -(intval(37)./intval(24))*m2_j/(1 + m2_j) - 1 + (1 + m2_j)*log(1 + m2_j)/m2_j;
    
    % Determine number of steps N(j)
    N_j = floor(sup(A_j - a2_j) *5*1e4) + 1;
    
    % Preallocate W{j}
    W{j} = zeros(1, N_j);
    
    % Loop over M-subdivision for fixed mY(j)
    for k = 1:N_j
        M_k = hull(a2_j + (k-1)*h, a2_j + k*h);  % rigorously enclose subinterval
        W{j}(k) = sup(dL(M_k, mY_j));           % compute sup of dL over sub-box
    end
end

% Compute max across all cells, max is equal to 0.6013
% Elapsed time is 220 seconds
maxval = max(cellfun(@max, W));
disp(maxval)
toc

%0.601280338274214, Elapsed time is 221.814207 seconds.