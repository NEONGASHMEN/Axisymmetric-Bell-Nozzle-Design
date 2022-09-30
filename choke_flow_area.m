function area = choke_flow_area(mdot,Pc,Tc,gamma,R)

gamma_term = ((gamma+1)/2)^((gamma+1)/(2*(1-gamma)));
area = mdot*sqrt(R*Tc/gamma)/(gamma_term*Pc);
end
