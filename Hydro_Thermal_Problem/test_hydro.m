
%addpath 'C:\Users\vince\Dropbox\Softwares\Dual_SDDP'  
addpath 'C:\Users\vince\Dropbox\Softwares\Dual_SDDP'
addpath 'C:\Program Files\Mosek\9.0'
addpath 'C:\Program Files\Mosek\9.0\toolbox\r2015a'

mat=dlmread('C:\Users\vince\Dropbox\Softwares\Dual_SDDP\data\deficit.txt');
%mat=dlmread('C:\Users\vince\Dropbox\Softwares\Dual_SDDP\data\deficit.txt');
%deficit_obj(i) is the cost of deficit in subsystem i
deficit_obj=mat(:,1);
%upper bound on fictitious plant production in subsystem i is deficit_depth(i)*demand(i) 
deficit_depth=mat(:,2);

%demands(t,i) is demand for stage t subsystem i

demands=dlmread('C:\Users\vince\Dropbox\Softwares\Dual_SDDP\data\demand.txt');

%number of stages
T=12;

%Number of scenarios per stage
M=50;

%upper_bound_exchanges(i,j) is the upper bound on the exchange between i and j
upper_bound_exchanges=dlmread('C:\Users\vince\Dropbox\Softwares\Dual_SDDP\data\exchange.txt');

%cost_exchange(i,j) is the exchange cost between i and j
cost_exchange=dlmread('C:\Users\vince\Dropbox\Softwares\Dual_SDDP\data\exchange_cost.txt');

%upper_bound_reservoirs(i) is the upper bound on reservoir of subsystem i level
upper_bound_reservoirs=dlmread('C:\Users\vince\Dropbox\Softwares\Dual_SDDP\data\upper_bound_reservoir.txt'); 
%initial_reservoir_level(i) is the initial level of reservoir i
initial_reservoir_level=dlmread('C:\Users\vince\Dropbox\Softwares\Dual_SDDP\data\initial_reservoir_level.txt');
initial_reservoir_level=0.6*initial_reservoir_level;
%initial_inflows(i) is inflows for first stage for reservoir i
initial_inflows=dlmread('C:\Users\vince\Dropbox\Softwares\Dual_SDDP\data\initial_inflows.txt');
%upper_hydro(i) is the upper bound on hydro production in subsystem i
upper_hydro=dlmread('C:\Users\vince\Dropbox\Softwares\Dual_SDDP\data\upper_hydro.txt');
%scenario_inflows_j(t,s) is inflow for scenario s stage t+1 subsystem j
scenario_inflows_1=dlmread('C:\Users\vince\Dropbox\Softwares\Dual_SDDP\data\scenario_inflows_0.txt');
scenario_inflows_2=dlmread('C:\Users\vince\Dropbox\Softwares\Dual_SDDP\data\scenario_inflows_1.txt');
scenario_inflows_3=dlmread('C:\Users\vince\Dropbox\Softwares\Dual_SDDP\data\scenario_inflows_2.txt');
scenario_inflows_4=dlmread('C:\Users\vince\Dropbox\Softwares\Dual_SDDP\data\scenario_inflows_3.txt');


thermal_0=dlmread('C:\Users\vince\Dropbox\Softwares\Dual_SDDP\data\thermal_0.txt');
thermal_1=dlmread('C:\Users\vince\Dropbox\Softwares\Dual_SDDP\data\thermal_1.txt');
thermal_2=dlmread('C:\Users\vince\Dropbox\Softwares\Dual_SDDP\data\thermal_2.txt');
thermal_3=dlmread('C:\Users\vince\Dropbox\Softwares\Dual_SDDP\data\thermal_3.txt');
%lower_bounds_thermal_j(i) is lower bound on thermal production of thermal
%plant i in subsystem j
lower_bounds_thermal_0=thermal_0(:,1);
lower_bounds_thermal_1=thermal_1(:,1);
lower_bounds_thermal_2=thermal_2(:,1);
lower_bounds_thermal_3=thermal_3(:,1);
lower_bounds_thermal=cell(1,4);
lower_bounds_thermal{1,1}=lower_bounds_thermal_0;
lower_bounds_thermal{1,2}=lower_bounds_thermal_1;
lower_bounds_thermal{1,3}=lower_bounds_thermal_2;
lower_bounds_thermal{1,4}=lower_bounds_thermal_3;

%upper_bounds_thermal_j(i) is lower bound on thermal production of thermal
%plant i in subsystem j
upper_bounds_thermal_0=thermal_0(:,2);
upper_bounds_thermal_1=thermal_1(:,2);
upper_bounds_thermal_2=thermal_2(:,2);
upper_bounds_thermal_3=thermal_3(:,2);
upper_bounds_thermal=cell(1,4);
upper_bounds_thermal{1,1}=upper_bounds_thermal_0;
upper_bounds_thermal{1,2}=upper_bounds_thermal_1;
upper_bounds_thermal{1,3}=upper_bounds_thermal_2;
upper_bounds_thermal{1,4}=upper_bounds_thermal_3;


%upper_bounds_thermal_j(i) is lower bound on thermal production of thermal
%plant i in subsystem j
cost_thermal_0=thermal_0(:,3);
cost_thermal_1=thermal_1(:,3);
cost_thermal_2=thermal_2(:,3);
cost_thermal_3=thermal_3(:,3);
cost_thermal=cell(1,4);
cost_thermal{1,1}=cost_thermal_0;
cost_thermal{1,2}=cost_thermal_1;
cost_thermal{1,3}=cost_thermal_2;
cost_thermal{1,4}=cost_thermal_3;

nb_thermal=[length(cost_thermal_0);length(cost_thermal_1);length(cost_thermal_2);length(cost_thermal_3)];

size_u=(2*4+sum(nb_thermal)+10+16)*ones(T,1);
size_x=4*ones(1,T);
size_b=9*ones(T,1);
size_f=zeros(T,1);
upper_u=cell(1,T);
lower_u=cell(1,T);
for t=1:T
    upper_u{1,t}=[];
    lower_u{1,t}=[];
    for i=1:4
        upper_u{1,t}=[upper_u{1,t};upper_hydro(i);inf];
        lower_u{1,t}=[lower_u{1,t};0;0];
        for j=1:nb_thermal(i)
            upper_u{1,t}=[upper_u{1,t};upper_bounds_thermal{1,i}(j)];
            lower_u{1,t}=[lower_u{1,t};lower_bounds_thermal{1,i}(j)];
        end
    end
    upper_u{1,t}=[upper_u{1,t};upper_bound_exchanges(1,2);upper_bound_exchanges(1,3);upper_bound_exchanges(1,5)];
    upper_u{1,t}=[upper_u{1,t};upper_bound_exchanges(2,1)];
    upper_u{1,t}=[upper_u{1,t};upper_bound_exchanges(3,1);upper_bound_exchanges(3,5)];
    upper_u{1,t}=[upper_u{1,t};upper_bound_exchanges(4,5)];
    upper_u{1,t}=[upper_u{1,t};upper_bound_exchanges(5,1);upper_bound_exchanges(5,3);upper_bound_exchanges(5,4)];
    for i=1:4
        for j=1:4
            upper_u{1,t}=[upper_u{1,t};deficit_depth(j)*demands(t,i)];
        end
    end
    lower_u{1,t}=[lower_u{1,t};zeros(26,1)];
end

lower_x=cell(1,T);
upper_x=cell(1,T);
for t=1:T
    lower_x{1,t}=zeros(4,1);
    for i=1:4
        upper_x{1,t}=[upper_x{1,t};upper_bound_reservoirs(i)];
    end
end

gamma=0.9906;
costc=cell(1,T);
costd=cell(1,T);
for t=1:T
    costc{1,t}=zeros(4,1);
    costd{1,t}=[];
    for i=1:4
        costd{1,t}=[costd{1,t};0;0.005];
        costd{1,t}=[costd{1,t};cost_thermal{1,i}];    
    end
    costd{1,t}=[costd{1,t};cost_exchange(1,2);cost_exchange(1,3);cost_exchange(1,5)];
    costd{1,t}=[costd{1,t};cost_exchange(2,1)];
    costd{1,t}=[costd{1,t};cost_exchange(3,1);cost_exchange(3,5)];
    costd{1,t}=[costd{1,t};cost_exchange(4,5)];
    costd{1,t}=[costd{1,t};cost_exchange(5,1);cost_exchange(5,3);cost_exchange(5,4)];
    for i=1:4
        for j=1:4
            costd{1,t}=[costd{1,t};deficit_depth(j)*demands(t,i)];
        end
    end
    costd{1,t}=(gamma^(t-1))*costd{1,t};
end

bs=cell(1,T);
bs{1,1}=cell(1,1);
fs=cell(1,T);
fs{1,1}=cell(1,1);
for t=2:T
    fs{1,t}=cell(1,M);
end
bs{1,1}{1,1}=[initial_inflows+initial_reservoir_level;demands(1,:)';0];
for t=2:T
    bs{1,t}=cell(1,M);
    for j=1:M
        bs{1,t}{1,j}=[scenario_inflows_1(t-1,j);scenario_inflows_2(t-1,j);scenario_inflows_3(t-1,j);scenario_inflows_4(t-1,j);demands(t,:)';0];
    end
end

probabilities=cell(1,T-1);
for t=1:T-1
    probabilities{1,t}=(1/50)*ones(1,M);
end

nb_scenarios_rhs=[1;50*ones(T-1,1)];
nb_iter_max=5000;
lower_pi=cell(1,T);
upper_pi=cell(1,T);
for t=1:T
    lower_pi{1,t}=-(10^6)*ones(9,1);
    upper_pi{1,t}=(10^6)*ones(9,1);
end

a_subi=cell(1,T);
a_subj=cell(1,T);
a_valij=cell(1,T);
for t=1:T
    a_subi{1,t}=[1,2,3,4];
    a_subj{1,t}=[1,2,3,4];
    a_valij{1,t}=[1,1,1,1];
end

b_subi=cell(1,T-1);
b_subj=cell(1,T-1);
b_valij=cell(1,T-1);
for t=1:T-1
    b_subi{1,t}=[1,2,3,4];
    b_subj{1,t}=[1,2,3,4];
    b_valij{1,t}=[-1,-1,-1,-1];
end

c_subi=cell(1,T);
c_subj=cell(1,T);
c_valij=cell(1,T);
nb_v=2+nb_thermal;
nb_var=sum(nb_v);

for t=1:T
    for i=1:4
        c_subi{1,t}=[c_subi{1,t},i,i];
        c_subj{1,t}=[c_subj{1,t},2*(i-1)+sum(nb_thermal(1:i-1))+1,2*(i-1)+sum(nb_thermal(1:i-1))+2];
        c_valij{1,t}=[c_valij{1,t},1,1];
    end
    
    %Subsystem 1
    c_subi{1,t}=[c_subi{1,t},5];
    c_subj{1,t}=[c_subj{1,t},1];
    c_valij{1,t}=[c_valij{1,t},1];
    
    c_subi{1,t}=[c_subi{1,t},5*ones(1,nb_thermal(1))];
    c_subj{1,t}=[c_subj{1,t},[3:1:2+nb_thermal(1)]];
    c_valij{1,t}=[c_valij{1,t},ones(1,nb_thermal(1))];
    
    c_subi{1,t}=[c_subi{1,t},5*ones(1,4)];
    c_subj{1,t}=[c_subj{1,t},[nb_var+11:1:nb_var+14]];
    c_valij{1,t}=[c_valij{1,t},ones(1,4)];
    
    c_subi{1,t}=[c_subi{1,t},5*ones(1,3)];
    c_subj{1,t}=[c_subj{1,t},nb_var+4,nb_var+5,nb_var+8];
    c_valij{1,t}=[c_valij{1,t},ones(1,3)];
    
    c_subi{1,t}=[c_subi{1,t},5*ones(1,3)];
    c_subj{1,t}=[c_subj{1,t},nb_var+1,nb_var+2,nb_var+3];
    c_valij{1,t}=[c_valij{1,t},-ones(1,3)];
    
    %Subsystem 2
    i=2;
    c_subi{1,t}=[c_subi{1,t},4+i];
    c_subj{1,t}=[c_subj{1,t},sum(nb_v(1:i-1))+1];
    c_valij{1,t}=[c_valij{1,t},1];
    
    c_subi{1,t}=[c_subi{1,t},(4+i)*ones(1,nb_thermal(i))];
    c_subj{1,t}=[c_subj{1,t},[sum(nb_v(1:i-1))+3:1:sum(nb_v(1:i-1))+2+nb_thermal(i)]];
    c_valij{1,t}=[c_valij{1,t},ones(1,nb_thermal(i))];
    
    c_subi{1,t}=[c_subi{1,t},(4+i)*ones(1,4)];
    c_subj{1,t}=[c_subj{1,t},[nb_var+10+4*(i-1)+1:nb_var+10+4*i]];
    c_valij{1,t}=[c_valij{1,t},ones(1,4)];
    
    c_subi{1,t}=[c_subi{1,t},4+i];
    c_subj{1,t}=[c_subj{1,t},nb_var+1];
    c_valij{1,t}=[c_valij{1,t},1];
    
    c_subi{1,t}=[c_subi{1,t},4+i];
    c_subj{1,t}=[c_subj{1,t},nb_var+4];
    c_valij{1,t}=[c_valij{1,t},-1];
    
    %Subsystem 3
    i=3;
    c_subi{1,t}=[c_subi{1,t},4+i];
    c_subj{1,t}=[c_subj{1,t},sum(nb_v(1:i-1))+1];
    c_valij{1,t}=[c_valij{1,t},1];
    
    c_subi{1,t}=[c_subi{1,t},(4+i)*ones(1,nb_thermal(i))];
    c_subj{1,t}=[c_subj{1,t},[sum(nb_v(1:i-1))+3:1:sum(nb_v(1:i-1))+2+nb_thermal(i)]];
    c_valij{1,t}=[c_valij{1,t},ones(1,nb_thermal(i))];
    
    c_subi{1,t}=[c_subi{1,t},(4+i)*ones(1,4)];
    c_subj{1,t}=[c_subj{1,t},[nb_var+10+4*(i-1)+1:nb_var+10+4*i]];
    c_valij{1,t}=[c_valij{1,t},ones(1,4)];
    
    c_subi{1,t}=[c_subi{1,t},4+i,4+i];
    c_subj{1,t}=[c_subj{1,t},nb_var+2,nb_var+9];
    c_valij{1,t}=[c_valij{1,t},1,1];
    
    c_subi{1,t}=[c_subi{1,t},4+i,4+i];
    c_subj{1,t}=[c_subj{1,t},nb_var+5,nb_var+6];
    c_valij{1,t}=[c_valij{1,t},-1,-1];
    
    %Subsystem 4
    i=4;
    c_subi{1,t}=[c_subi{1,t},4+i];
    c_subj{1,t}=[c_subj{1,t},sum(nb_v(1:i-1))+1];
    c_valij{1,t}=[c_valij{1,t},1];
    
    c_subi{1,t}=[c_subi{1,t},(4+i)*ones(1,nb_thermal(i))];
    c_subj{1,t}=[c_subj{1,t},[sum(nb_v(1:i-1))+3:1:sum(nb_v(1:i-1))+2+nb_thermal(i)]];
    c_valij{1,t}=[c_valij{1,t},ones(1,nb_thermal(i))];
    
    c_subi{1,t}=[c_subi{1,t},(4+i)*ones(1,4)];
    c_subj{1,t}=[c_subj{1,t},[nb_var+10+4*(i-1)+1:nb_var+10+4*i]];
    c_valij{1,t}=[c_valij{1,t},ones(1,4)];
    
    c_subi{1,t}=[c_subi{1,t},4+i];
    c_subj{1,t}=[c_subj{1,t},nb_var+10];
    c_valij{1,t}=[c_valij{1,t},1];
    
    c_subi{1,t}=[c_subi{1,t},4+i];
    c_subj{1,t}=[c_subj{1,t},nb_var+7];
    c_valij{1,t}=[c_valij{1,t},-1];
    
    %Subsystem 5
    
    i=5;
    c_subi{1,t}=[c_subi{1,t},4+i,4+i,4+i];
    c_subj{1,t}=[c_subj{1,t},nb_var+3,nb_var+6,nb_var+7];
    c_valij{1,t}=[c_valij{1,t},1,1,1];
    
    c_subi{1,t}=[c_subi{1,t},4+i,4+i,4+i];
    c_subj{1,t}=[c_subj{1,t},nb_var+8,nb_var+9,nb_var+10];
    c_valij{1,t}=[c_valij{1,t},-1,-1,-1];

    d_subi{1,t}=[];
    d_subj{1,t}=[];
    d_valij{1,t}=[];
    e_subi{1,t}=[];
    e_subj{1,t}=[];
    e_valij{1,t}=[];
    f_subi{1,t}=[];
    f_subj{1,t}=[];
    f_valij{1,t}=[];

end


itermax=5000;
penaltiesQ=(10^(9))*ones(T-1,itermax);
penaltiesR=(10^(9))*ones(T-1,itermax);
large_bound=10^6;
nb_iter_max=200;
[lower_bounds,upper_bounds,time]=primal_sddp_control_rhs_eq(T,size_x,size_u,size_b,size_f,upper_x,lower_x,upper_u,lower_u,costc,costd,bs,probabilities,nb_scenarios_rhs,1000,a_subi,a_subj,a_valij,b_subi,b_subj,b_valij,c_subi,c_subj,c_valij,0,1.96,0.01);
big_m=10^9;
nb_iter_max=1000;
[upper_bounds,time,Errors]=dual_sddp_control_rhs_eq(T,size_x,size_u,size_b,upper_x,lower_x,upper_u,lower_u,costc,costd,bs,probabilities,nb_scenarios_rhs,nb_iter_max,lower_pi,upper_pi,a_subi,a_subj,a_valij,b_subi,b_subj,b_valij,c_subi,c_subj,c_valij,big_m,penaltiesQ,penaltiesR);


