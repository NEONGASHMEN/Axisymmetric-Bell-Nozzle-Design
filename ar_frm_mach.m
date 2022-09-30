function ar = ar_frm_mach(M,gamma)
expTerm = 0.5*(gamma+1)/(gamma-1);
ar = (0.5*(gamma+1))^(-expTerm)*(1 + 0.5*(gamma-1)*M*M)^(expTerm)/M;;
end