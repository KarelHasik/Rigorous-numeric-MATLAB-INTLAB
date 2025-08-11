function z = m5RR(M)
% Computes a rigorous lower bound for m_5(M) using interval iteration
% using downward rounding and point interval enclosure

longprecision(32);

if ~isa(M, 'intval'); lo=pred(M,2); hi=succ(M,2); M = infsup(lo,hi); end

l = -M;  % Initial guess as interval
j = 1;

tic

% Save current rounding mode
rndold = getround;

while (sup(l) <= inf(LR(M, l))) && (j <= 5)
    setround(-1);                      % Downward rounding for safety
    l = intval(inf(LR(M, l)));         % Point interval that underestimates true value
    setround(rndold);                 % Restore rounding
    j = j + 1;
end

z = l;
toc
end
