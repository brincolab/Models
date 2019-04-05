function [joptim,curr_e,curr_i] = balance_j_fic(Tmax,dt,I0,Jexte,Jexti,w,JN,C,we,...
    gamma,sigma,taog,taon,wgaine,wgaini,receptors,g_e,g_i,Ie,Ii,ce,ci,delta,burnout)
%
% Finds optimal J for feedback inhibitory control that keeps the average
% firing rates of units close to 3 Hz
% Same parameters than dmf_deco2018, but J is to be estimated and user
% should plug the value of deltaJ, the step to modify J.
% Initial parameters
%
N=size(C,1); % number of units
J=ones(N,1); % To be estimated
delta=delta*ones(N,1); % deltaJ for optimization

for k=1:50000 % optimization loop
    [curr_e,curr_i] = dmf_deco18(Tmax,dt,I0,Jexte,Jexti,w,JN,C,J,we,...
        gamma,sigma,taog,taon,wgaine,wgaini,receptors,g_e,g_i,Ie,Ii,ce,ci); % running model
    curr_e_n = curr_e-125/310; % normalized
    currm=mean(curr_e_n(burnout:end,:),1); % unit average current, removes first 1000 points
    
    flag=0;
    
    for n=1:N % units loop. This is not the bottleneck aprox 0.05s for all units
        if abs(currm(n)+0.026)>0.005 % if averaged current - 125/310 exceeds -0.026 +- 0.005
            if currm(n)<-0.026 % if unit average current is below 0.026
                J(n)=J(n)-delta(n); % downregulte J by delta
                delta(n)=delta(n)-0.001; % reduce delta
                if delta(n)<0.001
                    delta(n)=0.001; % controls delta from going below 0.001
                end
            else % if unit average current is above 0.026
                J(n)=J(n)+delta(n); % upregulate J
            end
        else
            flag=flag+1;
        end
    end
    if flag==N % all units optimized
        break;
    end    
end
joptim = J;