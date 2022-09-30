%%Algorithm was developed thanks to paper published by Kyril Palaveev &
%%Dr. Torsten Schenkel. They've also made a similar program in Python;
%%https://github.com/KyrilPalaveev/Axisymmetric-Method-of-Characteristics-Nozzle-Designer-with-GUI

clear
clc
close all

%%IP
gamma = 1.165;      %Ratio of sp. heats
Pe = 101325;        %Optimise pressure
F = 5000;           %Thrust
Pc = 2270000;       %Chamber pressure
Tc = 1200;          %Chamber temperature
R = 367;            %Gas constant
g = 9.81;           %Acceleration due to gravity
L = 0.12;           %Characteristic net length
rsltn = 200;        %Net resolution
pds = "Parabolic";  %Pressure Distribution Scheme along nozzle axis
q = 0;              %q factor for Cubic Distribution
itr = 4;            %Iterations per net node
theta_max = 5;      %Divergence angle in degrees
plot_net = "No";    %Plot Kernel mesh

%%Calculating crucial params
pr = Pe/Pc;
Me = isentropic('var',"M",'pr',pr,'gamma',gamma);
tr = isentropic('var',"tr",'pr',pr,'gamma',gamma);
Te = tr*Tc;
a_e = sqrt(gamma*R*Te);
Ve = Me*a_e;
mdot = F/Ve;
Isp = F/(mdot*g);
At = choke_flow_area(mdot,Pc,Tc,gamma,R);
Rt = sqrt(At/pi);
Pt = Pc*isentropic('var',"pr",'gamma',gamma,'M',1);
Cv = R/(gamma-1);
Cp = R + Cv;
theta_max = theta_max*pi/180;
Ae = At*ar_frm_mach(Me,gamma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------------------MOC--------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Finding the pts inside the characteristic net
pts_no = rsltn + 1;
for i = rsltn:-1:2
    pts_no = pts_no + i;
end

%Create a matrix to store net params
params = struct('x',nan,'y',nan,'M',nan,'V',nan,'nu',nan,'mu',nan,'theta',nan,'Kn',nan,'Kp',nan,'P',nan,'T',nan);
net = repmat(params,1,pts_no);
net(1).x = 0;
net(1).y = 0;
net(1).M = 1.00000166;
net(1).P = Pc*isentropic('var',"pr",'M',net(1).M,'gamma',gamma);
net(1).T = Tc*isentropic('var',"tr",'gamma',gamma,'M',net(1).M);
net(1).V = (net(1).M)*sqrt(gamma*R*(net(1).T));
net(1).mu = mach_angle(net(1).M);
net(1).nu = pm(net(1).M,gamma);
net(1).theta = 0;
net(1).Kn = net(1).theta + net(1).nu;
net(1).Kp = net(1).theta - net(1).nu;

%Symmetry Line params
nxt_axis_pt = rsltn+1;
axis_pt = 1;
prev_x = 0;
axis_nodes = zeros(1,rsltn);
wall_nodes = zeros(1,rsltn);
axis_nodes(1) = 1;
wall_nodes(1) = 1+rsltn;
for i = 1:(rsltn-1)
    axis_pt = axis_pt + nxt_axis_pt;
    net(axis_pt).x = prev_x + (L/(rsltn-1));
    net(axis_pt).y = 0;
    net(axis_pt).P = pd(pds,net(axis_pt).x,L,Pt,Pe,q);
    net(axis_pt).M = isentropic('var',"M",'gamma',gamma,'pr',(net(axis_pt).P)/Pc);
    net(axis_pt).T = Tc*isentropic('var',"tr",'gamma',gamma,'M',net(axis_pt).M);
    net(axis_pt).V = (net(axis_pt).M)*sqrt(gamma*R*(net(axis_pt).T));
    net(axis_pt).mu = mach_angle(net(axis_pt).M);
    net(axis_pt).nu = pm(net(axis_pt).M,gamma);
    net(axis_pt).theta = 0;
    net(axis_pt).Kn = net(axis_pt).theta + net(axis_pt).nu;
    net(axis_pt).Kp = net(axis_pt).theta - net(axis_pt).nu;
    axis_nodes(i+1) = axis_pt;
    wall_nodes(i+1) = axis_pt + nxt_axis_pt - 2;
    prev_x = net(axis_pt).x;
    nxt_axis_pt = nxt_axis_pt - 1;
end
k = rsltn - 1;
row_nodes = cell(rsltn-1,1);
for i = 1:(rsltn -1)
    l = rsltn + 2;
    row_nodes{i,1} = zeros(1,k);
    for j = 1:k
        if (i == 1) && (j == 1)
            row_nodes{1,1}(1) = 2;
        elseif j == 1
            row_nodes{i,1}(j) = row_nodes{i-1,1}(j) + 1;
        else
            row_nodes{i,1}(j) = row_nodes{i,1}(j-1) + l;
        end
         l = l-1;
    end
    k = k-1;
end

%Interior nodes params
for h = 1:length(row_nodes)
    for i = row_nodes{h,1}

        if ~(ismember(i,axis_nodes)||ismember(i,wall_nodes))

            LR = net(i-1);
            for j = 1:length(axis_nodes)
                if i > axis_nodes(j)
                    parent_idx = j;
                end
            end
            RR = net(axis_nodes(parent_idx+1)+(i-axis_nodes(parent_idx)-1));

            for j = 1:itr
                if j == 1
                    slp_lr = tan(LR.theta + LR.mu);
                    slp_rr = tan(RR.theta - RR.mu);
                    term4x = RR.x*slp_rr - LR.x*slp_lr;
                    x_x = (LR.y - RR.y + term4x)/(slp_rr - slp_lr);
                    y_x = LR.y + (x_x - LR.x)*slp_lr;
                    term4V = 1/(cot(LR.mu)/LR.V + cot(RR.mu)/RR.V);
                    if (LR.y && RR.y) == 0
                        A = tan(LR.mu)*sin(LR.mu)*sin(LR.theta)/(y_x*cos(LR.theta+LR.mu));
                        B = tan(RR.mu)*sin(RR.mu)*sin(RR.theta)/(y_x*cos(RR.theta-RR.mu));
                    else
                        A = tan(LR.mu)*sin(LR.mu)*sin(LR.theta)/(LR.y*cos(LR.theta+LR.mu));
                        B = tan(RR.mu)*sin(RR.mu)*sin(RR.theta)/(RR.y*cos(RR.theta-RR.mu));
                    end
                    V_x = term4V*(cot(LR.mu)*(1+A*(x_x-LR.x))+cot(RR.mu)*(1+B*(x_x-RR.x))+RR.theta-LR.theta);
                    term4theta = cot(LR.mu)*(V_x-LR.V)/(LR.V);
                    theta_x = LR.theta + term4theta - A*(x_x-LR.x);
                    T_x = Tc - V_x*V_x/(2*Cp);
                    M_x = V_x/sqrt(gamma*R*T_x);
                    mu_x = mach_angle(M_x);
                    nu_x = pm(M_x,gamma);
                else
                    tanterm1 = 0.5*(tan(RR.theta-RR.mu) + tan(theta_x-mu_x));
                    tanterm2 = 0.5*(tan(LR.theta+LR.mu) + tan(theta_x+mu_x));
                    x_x2 = (LR.y - RR.y + RR.x*tanterm1 - LR.x*tanterm2)/(tanterm1-tanterm2);
                    term4y = 0.5*(tan(LR.theta+LR.mu) + tan(theta_x+mu_x));
                    y_x2 = LR.y + term4y*(x_x - LR.x);
                    if (LR.y && RR.y) == 0
                        LR_A = tan(LR.mu)*sin(LR.mu)*sin(LR.theta)/(y_x*cos(LR.theta+LR.mu));
                        RR_B = tan(RR.mu)*sin(RR.mu)*sin(RR.theta)/(y_x*cos(RR.theta-RR.mu));
                    else
                        LR_A = tan(LR.mu)*sin(LR.mu)*sin(LR.theta)/(LR.y*cos(LR.theta+LR.mu));
                        RR_B = tan(RR.mu)*sin(RR.mu)*sin(RR.theta)/(RR.y*cos(RR.theta-RR.mu));
                    end
                    A_x = tan(mu_x)*sin(mu_x)*sin(theta_x)/(y_x*cos(theta_x+mu_x));
                    B_x = tan(mu_x)*sin(mu_x)*sin(theta_x)/(y_x*cos(theta_x-mu_x));
                    term4N = 2/(tan(LR.mu) + tan(mu_x));
                    N = term4N*((2*LR.V/(LR.V+V_x)) + 0.5*(LR_A+A_x)*(x_x-LR.x));
                    term4C = 2/(tan(RR.mu) + tan(mu_x));
                    C = term4C*((2*RR.V/(RR.V+V_x)) + 0.5*(RR_B+B_x)*(x_x-RR.x));
                    term1dtr = (1/(LR.V+V_x))*(1/(tan(LR.mu)+tan(mu_x)));
                    term2dtr = (1/(RR.V+V_x))*(1/(tan(RR.mu)+tan(mu_x)));
                    V_x2 = (N + C - LR.theta + RR.theta)/(4*(term1dtr+term2dtr));
                    term1_4theta = 2*(V_x2 - LR.V)/(V_x2 + LR.V);
                    term2_4theta = 2/(tan(LR.mu) + tan(mu_x));
                    theta_x2 = LR.theta + term2_4theta*(term1_4theta - 0.5*(LR_A+A_x)*(x_x-LR.x));
                    T_x2 = Tc - V_x2*V_x2/(2*Cp);
                    M_x2 = V_x2/sqrt(gamma*R*T_x2);
                    mu_x2 = mach_angle(M_x2);
                    nu_x2 = pm(M_x2,gamma);
                    x_x = x_x2;
                    y_x = y_x2;
                    V_x = V_x2;
                    theta_x = theta_x2;
                    M_x = M_x2;
                    T_x = T_x2;
                    nu_x = nu_x2;
                    mu_x = mu_x2;
                end         
            end

            net(i).x = x_x;
            net(i).y = y_x;
            net(i).V = V_x;
            net(i).theta = theta_x;
            net(i).T = T_x;
            net(i).M = M_x;
            net(i).nu = nu_x;
            net(i).mu = mu_x;
        end
    end
end

%Wall contour params
j = 1;
for i = wall_nodes
    net(i).theta = net(i-1).theta;
    net(i).mu = net(i-1).mu;
    net(i).M = net(i-1).M;
    net(i).V = net(i-1).V;
    net(i).nu = net(i-1).nu;
    net(i).T = net(i-1).T;
    if j == 1
        term1 = tan(net(i).theta + net(i).mu);
        term2 = tan(theta_max);
        numrtr = net(i-1).y - net(i-1).x*term1 - Rt;
        dnmtr = term2 - term1;
        net(i).x = numrtr/dnmtr;
        net(i).y = net(i-1).y + term1*(net(i).x - net(i-1).x);
    elseif j == 2
        term1 = tan(net(i).theta + net(i).mu);
        term2 = tan(0.5*(theta_max + net(wall_nodes(j-1)).theta));
        numrtr = net(i-1).y - net(i-1).x*term1 + net(wall_nodes(j-1)).x*term2 - net(wall_nodes(j-1)).y;
        dnmtr = term2 - term1;
        net(i).x = numrtr/dnmtr;
        net(i).y = net(i-1).y + term1*(net(i).x - net(i-1).x);
    else  
        term1 = tan(net(i).theta + net(i).mu);
        term2 = tan(0.5*(net(wall_nodes(j-2)).theta + net(wall_nodes(j-1)).theta));
        numrtr = net(i-1).y - net(i-1).x*term1 + net(wall_nodes(j-1)).x*term2 - net(wall_nodes(j-1)).y;
        dnmtr = term2 - term1;
        net(i).x = numrtr/dnmtr;
        net(i).y = net(i-1).y + term1*(net(i).x - net(i-1).x);
    end
    j = j + 1;
end

%Print Stuff
Ae_numerical = pi*net(end).y*net(end).y;
disp(['Exit Area: ' num2str(Ae_numerical*1e6) ' mm2']);
disp(['Isentropic Exit Area: ' num2str(Ae*1e6) ' mm2']);
disp("  ");
XP_numerical = Ae_numerical/At;
XP_ideal = Ae/At;
disp(['XP ratio: ' num2str(XP_numerical)]);
disp(['XP ratio (Ideal): ' num2str(XP_ideal)]);
disp(['Percentage Error: ' num2str((XP_ideal - XP_numerical)/XP_ideal)]);
disp(" ");
disp(['Length of Nozzle: ' num2str(net(end).x*1e3) ' mm']);

%Export Stuff
crdnts = zeros(length(wall_nodes)+1,3);
j = 2;
crdnts(1,1) = 0;
crdnts(1,2) = Rt;
for i = wall_nodes
    crdnts(j,1) = net(i).x;
    crdnts(j,2) = net(i).y;
    j = j + 1;
end
writematrix(crdnts*1e3,'contour.dat');


%Plot Contour
fig1 = figure();
axis([0 0.25 0 0.25]);
for i = 1:length(wall_nodes)
    if i == 1
        line([0 net(wall_nodes(i)).x],[Rt net(wall_nodes(i)).y]);
    else
        line([net(wall_nodes(i-1)).x net(wall_nodes(i)).x],[net(wall_nodes(i-1)).y net(wall_nodes(i)).y]);
        hold on;
    end
end

%Plot net
if (plot_net == "Yes") || (plot_net == "yes")
    if rsltn > 49
        disp("No of Characteristic Lines are high. Plotting the mesh will take time. Please be patient !");
    end
    j = rsltn - 1;
    for i = 1:length(row_nodes)
        for j = 1:length(row_nodes{i})
            if i == 1
                x = net(row_nodes{i}(j)).x;
                lx = net(axis_nodes(j)).x;
                rx = net(axis_nodes(j+1)).x;
                y = net(row_nodes{i}(j)).y;
                ly = net(axis_nodes(j)).y;
                ry = net(axis_nodes(j+1)).y;
                line([lx x],[ly y],'Color','k');
                line([rx x],[ry y],'Color','k');
            else
                x = net(row_nodes{i}(j)).x;
                lx = net(row_nodes{i-1}(j)).x;
                rx = net(row_nodes{i-1}(j+1)).x;
                y = net(row_nodes{i}(j)).y;
                ly = net(row_nodes{i-1}(j)).y;
                ry = net(row_nodes{i-1}(j+1)).y;
                line([lx x],[ly y],'Color','k');
                line([rx x],[ry y],'Color','k');
            end
        end
    end
end








