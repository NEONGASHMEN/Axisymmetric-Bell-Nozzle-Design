function result = pd(dis_mode,x,L,Pt,Pe,q)

if dis_mode == "Parabolic"
    result = exp((L^(-2))*((x-L)^2)*(log(Pt)-log(Pe)) + log(Pe));
elseif dis_mode == "Cubic"
    big_term1 = (q*L + 2*(log(Pt)-log(Pe)))*(L^(-3));
    big_term2 = (2*q*L + 3*(log(Pt)-log(Pe)))*(L^(-2));
    result = exp(((x^3)*big_term1) -  ((x^2)*big_term2) + q*x + log(Pt));
elseif dis_mode == "Cosine"
    log_term1 = log(Pt) - log(Pe);
    log_term2 = log(Pt) + log(Pe);
    result = 0.5*log_term1*cos(pi*x/L) + 0.5*log_term2;
else
    disp("Please provide the Axis Pressure Distribution scheme !");
end