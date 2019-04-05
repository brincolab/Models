function [curr_e,curr_i] = dmf_deco18(Tmax,dt,I0,Jexte,Jexti,w,JN,C,J,we,...
    gamma,sigma,taog,taon,wgaine,wgaini,receptors,g_e,g_i,Ie,Ii,ce,ci)
% Implements Deco 2018 Curr Biol Dynamic Mean Field Model
%
% INPUTS
% Tmax = simulation time points
% dt = integration step [ms]
% I0 = external current [nA]
% Jexte = external->E coupling
% Jexte = external->I coupling
% w = local excitatory recurrence
% JN = excitatory synaptic coupling [nA]
% C = connectivity matrix, should fit length of receptors
% J = feedback inhibitory control parameter to be estimated
% we = Global Coupling parameter
% gamma = Kinetic Parameter of Excitation
% sigma = amplitude of noise [nA]
% taog = GABA time constant
% taon = NMDA time constant
% wgaine = excitatory gain modulation
% wgaini = excitatory gain modulation
% receptors = receptor density per region
% g_e = excitatory conductance
% g_i = same than above but inhibitory
% Ie = excitatory threshold for nonlineariy
% Ii = same than above but inhibitory
% ce = excitatory non linear shape parameter
% ci = same than above but inhibitory
%
% OUTPUTS
% curr_e = currents of excitatory pools per region
% curr_i = currents of inhibitory pools per region

% Initial parameters
N=size(C,1); % number of units

% Model simulations
curr_e=zeros(Tmax,N);
curr_i = curr_e;
sn=0.001*ones(N,1);
sg=0.001*ones(N,1);
for t=1:Tmax  % integration in time
    xn=I0*Jexte+w*JN*sn+we*JN*C*sn-J.*sg;% excitatory currents
    xg=I0*Jexti+JN*sn-sg;% inhibitory currents    
    rn = curr2rate(xn,wgaine,g_e,Ie,ce,receptors); % excitatory population rate function        
    rg = curr2rate(xg,wgaini,g_i,Ii,ci,receptors);% inhibitory population rate function    
    sn=sn+dt*(-sn/taon+(1-sn)*gamma.*rn./1000.)+sqrt(dt)*sigma*randn(N,1);
    sn(sn>1) = 1;
    sn(sn<0) = 0;
    sg=sg+dt*(-sg/taog+rg./1000.)+sqrt(dt)*sigma*randn(N,1);
    sg(sg>1) = 1;
    sg(sg<0) = 0;
    % storing variables
    curr_e(t,:) = xn;
    curr_i(t,:) = xg;
    
end