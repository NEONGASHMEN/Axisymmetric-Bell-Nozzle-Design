function result = isentropic(args)
arguments
    args.M = '-';
    args.gamma = '-';
    args.pr = '-';
    args.tr = '-';
    args.var;
end

%%vpasolve() gives +0.00i << disrupts calcs
if args.gamma == '-'
    disp("Please provide a value for gamma !")
else
    if args.var == "M"
        if args.pr ~= '-'
            result = mach_frm_pr(args.pr,args.gamma);
        elseif args.tr ~= '-'
            args.pr = args.tr^(args.gamma/(args.gamma-1));
            result = mach_frm_pr(args.pr,args.gamma);
        else
            disp("Please enter pr or tr !");
        end
    end

    if args.var == "pr"
        if args.M ~= '-'
            result = pr_frm_mach(args.M,args.gamma);
            %syms x
            %result = double(vpasolve(mach_frm_pr(x,args.gamma)==args.M,x));
        elseif args.tr ~= '-'
            result = args.tr^(args.gamma/(args.gamma-1));
        else
            disp("Please enter M or tr !");
        end
    end

    if args.var == "tr"
        if args.pr ~= '-'
            result = args.pr^((args.gamma-1)/args.gamma);
        elseif args.M ~= '-'
            args.pr = pr_frm_mach(args.M,args.gamma);
            %syms x
            %args.pr = double(vpasolve(mach_frm_pr(x,args.gamma)==args.M,x));
            result = args.pr^((args.gamma-1)/args.gamma);
        else
            disp("Please enter M or pr !");
        end
    end
end


% if args.var == "M"
%     result = mach_frm_pr(args.pr,args.gamma);
% end
% if args.var == "pr"
%     syms x
%     result = double(vpasolve(mach_frm_pr(x,args.gamma)==args.M,x));
% end
% if args.var == "gamma"
%     syms x
%     %result = double(vpasolve(mach_frm_pr(args.pr,x)==args.M,x));
%     disp("Not possible to find the value of gamma using vpasolve");
% end
% if args.var == "tr"
%     result = args.pr^((gamma-1)/gamma);
end
    