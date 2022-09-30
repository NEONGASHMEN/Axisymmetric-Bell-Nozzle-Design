function M = mach_frm_pr(pr,gamma)
temp_gamma = (1-gamma)/gamma;
M_sq = (2/(gamma-1))*((pr)^(temp_gamma) - 1);
M = sqrt(M_sq);                                 %%isentropic flow rltn
end