function nu = pm(M,gamma)
temp_gamma = sqrt((gamma+1)/(gamma-1));
temp_M = sqrt(M*M - 1);
nu = (temp_gamma*atan(inv(temp_gamma)*temp_M)) - atan(temp_M);
end