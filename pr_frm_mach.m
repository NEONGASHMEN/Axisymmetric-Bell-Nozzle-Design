function pr = pr_frm_mach(M,gamma)
gamma_term1 = 0.5*(gamma-1);
gamma_term2 = gamma/(1-gamma);
pr = (1 + gamma_term1*(M^2))^gamma_term2;
end