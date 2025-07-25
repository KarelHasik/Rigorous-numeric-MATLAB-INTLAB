longprecision(32);
tstart = tic;

rndold = getround;  % Save rounding mode

M_0 = intval('0.377');
m_0 = m5RR(M_0);

% Compute rigorous enclosure of sup(theta6R(m_0))
setround(1);                         % upward rounding
M1_sup = sup(theta6R(m_0));          % safely rounded up
M_1 = infsup(M1_sup, M1_sup);        % point interval that contains it
setround(rndold);                    % restore rounding mode

j = 1;

while (M_1 < M_0) && (inf(M_1) > 0.0094)
    j = j + 1;
    M_0 = M_1;
    m_0 = m5RR(M_0);                   % m5R returns interval

    setround(1);
    M1_sup = sup(theta6R(m_0));
    M_1 = infsup(M1_sup, M1_sup);    % construct safe enclosure
    setround(rndold);
end

% Output results for [0.1005,0.377] j =63, intval M_0 = 0.10067100489005, intval M_1 = 0.09991953654082, 
% intval m_0 = -0.08951396624275, telapsed =19.332572500000002
% [0.05,0.1007] j =209, intval M_0 = 0.05005586376854,intval M_1 = 0.04996841011585, 
% intval m_0 = -0.04718651638436, telapsed = 57.207350583333338
%[0.04,0.05006] j =166, intval M_0 = 0.04001741563917, intval M_1 =
%0.03997469071316, intval m_0 = -0.03816908020074, telapsed = 46.824151791666673
%[0.03,0.04002] j =382, intval M_0 = 0.03000548386822, intval M_1 =
%0.02998898015919, intval m_0 = -0.02895814763688, telapsed =1.212718312083334e+02
%[0.02,0.03001] j =1273, intval M_0 = 0.02000072061956, intval M_1 = 0.01999678449971
%intval m_0 = -0.01953171704053, telapsed =3.693981672916668e+02
%[0.0094,0.377] j =24336, intval M_0 = 0.00940004991117, intval M_1 = 
% 0.00939998983810, intval m_0 = -0.00929559077686, telapsed =7.052058126916667e+03
j
M_0
M_1
m_0
telapsed = toc(tstart)
