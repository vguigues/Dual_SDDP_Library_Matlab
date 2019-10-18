

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%T: number of stages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%size_x(t): number of x variables for stage t, t=1,..,T.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%size_u(t): number of u variables for stage t, t=1,..,T.---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%size_b(t): number of components in b(t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%upper_x: upper_x{1,t} is the upper bound on x(t) for stage t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%lower_x: lower_x{1,t} is the lower bound on x(t) for stage t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%upper_u: upper_u{1,t} is the upper bound on u(t) for stage t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%lower_u: lower_u{1,t} is the lower bound on u(t) for stage t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%costc{1,t}: cost vector c for stage t.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%costd{1,t}: cost vector d for stage t.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bs{1,t}{1,j}: rhs b_t for stage t scenario j.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The technology matrix A_{t} for stage t is given by
%a_subi{1,t}, a_subj{1,t}, a_valij{1,t}:
%A_{t}(a_subi{1,t}(k),a_subj{1,t}(k))=a_valij{1,t}(k).
%Similarly for matrices B_{t}, C_{t}, D_{t}, E_{t}, F_{t} for stage t>=2. For instance B_t is 
%given by b_subi{1,t-1}, b_subj{1,t-1}, b_valij{1,t-1}:
%B_{t}(b_subi{1,t-1}(k),b_subj{1,t-1}(k))= b_valij{1,t-1}(k).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%probabilities{1,t-1}(j) is the probability of scenario j for stage t=2,..,T.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nb_scenarios_rhs(t): number of scenarios for stage t=1,..,T.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nb_iter_max: maximal number of iterations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%lower_pi: lower_pi{1,t} lower bound on pi(t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%upper_pi: upper_pi{1,t} upper_bound on pi(t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%upper_bounds(k) is the upper bound for iteration k.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%time(k) is the CPU time needed to solve iteration k of the problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [upper_bounds,time,Errors]=dual_sddp_control_rhs_eq(T,size_x,size_u,size_b,upper_x,lower_x,upper_u,lower_u,costc,costd,bs,probabilities,nb_scenarios_rhs,nb_iter_max,lower_pi,upper_pi,a_subi,a_subj,a_valij,b_subi,b_subj,b_valij,c_subi,c_subj,c_valij,big_m,penaltiesQ,penaltiesR)

upper_bounds=[];
time=[];
Errors=cell(1,T-1);

dynamic_subi_af=cell(1,T);
dynamic_subj_af=cell(1,T);
dynamic_valij_af=cell(1,T);

dynamic_blx=cell(1,T);
dynamic_bux=cell(1,T);

dynamic_buc=cell(1,T);
dynamic_blc=cell(1,T);

dynamic_c=cell(1,T);

%Initializations 

% Stage 1

dynamic_subi_af{1,1}=[c_subj{1,1}];
dynamic_subj_af{1,1}=[c_subi{1,1}];
dynamic_valij_af{1,1}=[c_valij{1,1}];

dynamic_subi_af{1,1}=[dynamic_subi_af{1,1},[1:size_u(1)],[1:size_u(1)]];
dynamic_subj_af{1,1}=[dynamic_subj_af{1,1},[size_b(1)+2*size_x(1)+1:size_b(1)+2*size_x(1)+2*size_u(1)]];
dynamic_valij_af{1,1}=[dynamic_valij_af{1,1},ones(1,2*size_u(1))];

dynamic_subi_af{1,1}=[dynamic_subi_af{1,1},size_u(1)+1];
dynamic_subj_af{1,1}=[dynamic_subj_af{1,1},size_b(1)+2*(size_x(1)+size_u(1))+1];
dynamic_valij_af{1,1}=[dynamic_valij_af{1,1},1];

dynamic_buc{1,1}=[costd{1,1};big_m];
dynamic_blc{1,1}=[costd{1,1};-inf];

dynamic_blx{1,1}=[lower_pi{1,1}];
dynamic_bux{1,1}=[upper_pi{1,1}];

dynamic_c{1,1}=[bs{1,1}{1,1}];
for i=1:size_x(1)
    if (upper_x{1,1}(i)==inf)
       dynamic_c{1,1}=[dynamic_c{1,1};0];
       dynamic_blx{1,1}=[dynamic_blx{1,1};0];
       dynamic_bux{1,1}=[dynamic_bux{1,1};0];
    else
       dynamic_c{1,1}=[dynamic_c{1,1};upper_x{1,1}(i)];
       dynamic_blx{1,1}=[dynamic_blx{1,1};-inf];
       dynamic_bux{1,1}=[dynamic_bux{1,1};0];
    end
end

for i=1:size_x(1)
    if (lower_x{1,1}(i)==-inf)
       dynamic_c{1,1}=[dynamic_c{1,1};0];
       dynamic_blx{1,1}=[dynamic_blx{1,1};0];
       dynamic_bux{1,1}=[dynamic_bux{1,1};0];
    else
       dynamic_c{1,1}=[dynamic_c{1,1};lower_x{1,1}(i)];
       dynamic_blx{1,1}=[dynamic_blx{1,1};0];
       dynamic_bux{1,1}=[dynamic_bux{1,1};inf];
    end
end

for i=1:size_u(1)
    if (upper_u{1,1}(i)==inf)
       dynamic_c{1,1}=[dynamic_c{1,1};0];
       dynamic_blx{1,1}=[dynamic_blx{1,1};0];
       dynamic_bux{1,1}=[dynamic_bux{1,1};0];
    else
       dynamic_c{1,1}=[dynamic_c{1,1};upper_u{1,1}(i)];
       dynamic_blx{1,1}=[dynamic_blx{1,1};-inf];
       dynamic_bux{1,1}=[dynamic_bux{1,1};0];
    end
end

for i=1:size_u(1)
    if (lower_u{1,1}(i)==-inf)
       dynamic_c{1,1}=[dynamic_c{1,1};0];
       dynamic_blx{1,1}=[dynamic_blx{1,1};0];
       dynamic_bux{1,1}=[dynamic_bux{1,1};0];
    else
       dynamic_c{1,1}=[dynamic_c{1,1};lower_u{1,1}(i)];
       dynamic_blx{1,1}=[dynamic_blx{1,1};0];
       dynamic_bux{1,1}=[dynamic_bux{1,1};inf];
    end
end

dynamic_blx{1,1}=[dynamic_blx{1,1};-inf];
dynamic_bux{1,1}=[dynamic_bux{1,1};inf];
dynamic_c{1,1}=[dynamic_c{1,1};1];

% Stage T

dynamic_subi_af{1,T}=[];
dynamic_subj_af{1,T}=[];
dynamic_valij_af{1,T}=[];

for i=1:nb_scenarios_rhs(T)
    dynamic_subi_af{1,T}=[dynamic_subi_af{1,T},a_subj{1,T}+(i-1)*size_x(T),c_subj{1,T}+(i-1)*size_u(T)+nb_scenarios_rhs(T)*size_x(T)];
    dynamic_subj_af{1,T}=[dynamic_subj_af{1,T},a_subi{1,T}+(i-1)*(size_b(T)+2*(size_x(T)+size_u(T))),c_subi{1,T}+(i-1)*(size_b(T)+2*(size_x(T)+size_u(T)))];
    dynamic_valij_af{1,T}=[dynamic_valij_af{1,T},a_valij{1,T},c_valij{1,T}];
end

for i=1:nb_scenarios_rhs(T)
    dynamic_subi_af{1,T}=[dynamic_subi_af{1,T},b_subj{1,T-1}+nb_scenarios_rhs(T)*(size_x(T)+size_u(T))];
    dynamic_subj_af{1,T}=[dynamic_subj_af{1,T},b_subi{1,T-1}+(i-1)*(size_b(T)+2*(size_x(T)+size_u(T)))];
    dynamic_valij_af{1,T}=[dynamic_valij_af{1,T},probabilities{1,T-1}(i)*b_valij{1,T-1}];
end

for i=1:nb_scenarios_rhs(T)
    dynamic_subi_af{1,T}=[dynamic_subi_af{1,T},[(i-1)*size_x(T)+1:i*size_x(T)]];
    dynamic_subj_af{1,T}=[dynamic_subj_af{1,T},[(i-1)*(size_b(T)+2*(size_x(T)+size_u(T)))+size_b(T)+1:(i-1)*(size_b(T)+2*(size_x(T)+size_u(T)))+size_b(T)+size_x(T)]];
    dynamic_valij_af{1,T}=[dynamic_valij_af{1,T},ones(1,size_x(T))];
    
    dynamic_subi_af{1,T}=[dynamic_subi_af{1,T},[(i-1)*size_x(T)+1:i*size_x(T)]];
    dynamic_subj_af{1,T}=[dynamic_subj_af{1,T},[(i-1)*(size_b(T)+2*(size_x(T)+size_u(T)))+size_b(T)+size_x(T)+1:(i-1)*(size_b(T)+2*(size_x(T)+size_u(T)))+size_b(T)+2*size_x(T)]];
    dynamic_valij_af{1,T}=[dynamic_valij_af{1,T},ones(1,size_x(T))];
    
    dynamic_subi_af{1,T}=[dynamic_subi_af{1,T},[nb_scenarios_rhs(T)*size_x(T)+(i-1)*size_u(T)+1:nb_scenarios_rhs(T)*size_x(T)+i*size_u(T)]];
    dynamic_subj_af{1,T}=[dynamic_subj_af{1,T},[(i-1)*(size_b(T)+2*(size_x(T)+size_u(T)))+size_b(T)+2*size_x(T)+1:(i-1)*(size_b(T)+2*(size_x(T)+size_u(T)))+size_b(T)+2*size_x(T)+size_u(T)]];
    dynamic_valij_af{1,T}=[dynamic_valij_af{1,T},ones(1,size_u(T))];
    
    dynamic_subi_af{1,T}=[dynamic_subi_af{1,T},[nb_scenarios_rhs(T)*size_x(T)+(i-1)*size_u(T)+1:nb_scenarios_rhs(T)*size_x(T)+i*size_u(T)]];
    dynamic_subj_af{1,T}=[dynamic_subj_af{1,T},[(i-1)*(size_b(T)+2*(size_x(T)+size_u(T)))+size_b(T)+2*size_x(T)+1+size_u(T):i*(size_b(T)+2*(size_x(T)+size_u(T)))]];
    dynamic_valij_af{1,T}=[dynamic_valij_af{1,T},ones(1,size_u(T))];    
end

dynamic_subi_af{1,T}=[dynamic_subi_af{1,T},[nb_scenarios_rhs(T)*(size_x(T)+size_u(T))+1:nb_scenarios_rhs(T)*(size_x(T)+size_u(T))+size_x(T-1)]];
dynamic_subj_af{1,T}=[dynamic_subj_af{1,T},[nb_scenarios_rhs(T)*(size_b(T)+2*(size_x(T)+size_u(T)))+1:nb_scenarios_rhs(T)*(size_b(T)+2*(size_x(T)+size_u(T)))+size_x(T-1)]];
dynamic_valij_af{1,T}=[dynamic_valij_af{1,T},ones(1,size_x(T-1))];

dynamic_subi_af{1,T}=[dynamic_subi_af{1,T},[nb_scenarios_rhs(T)*(size_x(T)+size_u(T))+1:nb_scenarios_rhs(T)*(size_x(T)+size_u(T))+size_x(T-1)]];
dynamic_subj_af{1,T}=[dynamic_subj_af{1,T},[nb_scenarios_rhs(T)*(size_b(T)+2*(size_x(T)+size_u(T)))+size_x(T-1)+1:nb_scenarios_rhs(T)*(size_b(T)+2*(size_x(T)+size_u(T)))+2*size_x(T-1)]];
dynamic_valij_af{1,T}=[dynamic_valij_af{1,T},-ones(1,size_x(T-1))];

dynamic_buc{1,T}=[repmat(costc{1,T},nb_scenarios_rhs(T),1);repmat(costd{1,T},nb_scenarios_rhs(T),1);zeros(size_x(T-1),1)];
dynamic_blc{1,T}=dynamic_buc{1,T};

dynamic_c{1,T}=[];
dynamic_blx{1,T}=[];
dynamic_bux{1,T}=[];

for i=1:nb_scenarios_rhs(T)
    dynamic_blx{1,T}=[dynamic_blx{1,T};lower_pi{1,T}];
    dynamic_bux{1,T}=[dynamic_bux{1,T};upper_pi{1,T}];
    dynamic_c{1,T}=[dynamic_c{1,T};probabilities{1,T-1}(i)*bs{1,T}{1,i}];    
    for j=1:size_x(T)
        if (upper_x{1,T}(j)==inf)
           dynamic_c{1,T}=[dynamic_c{1,T};0];
           dynamic_blx{1,T}=[dynamic_blx{1,T};0];
           dynamic_bux{1,T}=[dynamic_bux{1,T};0];
        else
           dynamic_c{1,T}=[dynamic_c{1,T};probabilities{1,T-1}(i)*upper_x{1,T}(j)];
           dynamic_blx{1,T}=[dynamic_blx{1,T};-inf];
           dynamic_bux{1,T}=[dynamic_bux{1,T};0];
        end
    end

    for j=1:size_x(T)
        if (lower_x{1,T}(j)==-inf)
           dynamic_c{1,T}=[dynamic_c{1,T};0];
           dynamic_blx{1,T}=[dynamic_blx{1,T};0];
           dynamic_bux{1,T}=[dynamic_bux{1,T};0];
        else
           dynamic_c{1,T}=[dynamic_c{1,T};probabilities{1,T-1}(i)*lower_x{1,T}(j)];
           dynamic_blx{1,T}=[dynamic_blx{1,T};0];
           dynamic_bux{1,T}=[dynamic_bux{1,T};inf];
        end
    end

    for j=1:size_u(T)
        if (upper_u{1,T}(j)==inf)
           dynamic_c{1,T}=[dynamic_c{1,T};0];
           dynamic_blx{1,T}=[dynamic_blx{1,T};0];
           dynamic_bux{1,T}=[dynamic_bux{1,T};0];
        else
           dynamic_c{1,T}=[dynamic_c{1,T};probabilities{1,T-1}(i)*upper_u{1,T}(j)];
           dynamic_blx{1,T}=[dynamic_blx{1,T};-inf];
           dynamic_bux{1,T}=[dynamic_bux{1,T};0];
        end
    end

    for j=1:size_u(T)
        if (lower_u{1,T}(j)==-inf)
           dynamic_c{1,T}=[dynamic_c{1,T};0];
           dynamic_blx{1,T}=[dynamic_blx{1,T};0];
           dynamic_bux{1,T}=[dynamic_bux{1,T};0];
        else
           dynamic_c{1,T}=[dynamic_c{1,T};probabilities{1,T-1}(i)*lower_u{1,1}(j)];
           dynamic_blx{1,T}=[dynamic_blx{1,T};0];
           dynamic_bux{1,T}=[dynamic_bux{1,T};inf];
        end
    end
end

dynamic_blx{1,T}=[dynamic_blx{1,T};zeros(2*size_x(T-1),1)];
dynamic_bux{1,T}=[dynamic_bux{1,T};inf*ones(2*size_x(T-1),1)];
dynamic_c{1,T}=[dynamic_c{1,T};-penaltiesQ(T-1,1)*ones(size_x(T-1),1);-penaltiesR(T-1,1)*ones(size_x(T-1),1)];

% Stages t=2,...,T-1

for t=2:T-1
    dynamic_subi_af{1,t}=[];
    dynamic_subj_af{1,t}=[];
    dynamic_valij_af{1,t}=[];
    dimC=size_b(t)+2*(size_x(t)+size_u(t));
    for i=1:nb_scenarios_rhs(t)
        dynamic_subi_af{1,t}=[dynamic_subi_af{1,t},(i-1)*size_u(t)+c_subj{1,t}];
        dynamic_subj_af{1,t}=[dynamic_subj_af{1,t},(i-1)*dimC+c_subi{1,t}];
        dynamic_valij_af{1,t}=[dynamic_valij_af{1,t},c_valij{1,t}];
        
        
        dynamic_subi_af{1,t}=[dynamic_subi_af{1,t},[(i-1)*size_u(t)+1:i*size_u(t)]];
        dynamic_subj_af{1,t}=[dynamic_subj_af{1,t},[(i-1)*dimC+size_b(t)+2*size_x(t)+1:(i-1)*dimC+size_b(t)+2*size_x(t)+size_u(t)]];
        dynamic_valij_af{1,t}=[dynamic_valij_af{1,t},ones(1,size_u(t))];
        
        dynamic_subi_af{1,t}=[dynamic_subi_af{1,t},[(i-1)*size_u(t)+1:i*size_u(t)]];
        dynamic_subj_af{1,t}=[dynamic_subj_af{1,t},[(i-1)*dimC+size_b(t)+2*size_x(t)+size_u(t)+1:(i-1)*dimC+size_b(t)+2*(size_x(t)+size_u(t))]];
        dynamic_valij_af{1,t}=[dynamic_valij_af{1,t},ones(1,size_u(t))];
        
        dynamic_subi_af{1,t}=[dynamic_subi_af{1,t},nb_scenarios_rhs(t)*size_u(t)+b_subj{1,t}];
        dynamic_subj_af{1,t}=[dynamic_subj_af{1,t},(i-1)*dimC+b_subi{1,t}];
        dynamic_valij_af{1,t}=[dynamic_valij_af{1,t},probabilities{1,t-1}(i)*b_valij{1,t}];
                
        dynamic_subi_af{1,t}=[dynamic_subi_af{1,t},[nb_scenarios_rhs(t)*size_u(t)+1:nb_scenarios_rhs(t)*size_u(t)+size_x(t-1)]];
        dynamic_subj_af{1,t}=[dynamic_subj_af{1,t},[nb_scenarios_rhs(t)*(dimC+1)+1:nb_scenarios_rhs(t)*(dimC+1)+size_x(t-1)]];
        dynamic_valij_af{1,t}=[dynamic_valij_af{1,t},ones(1,size_x(t-1))];
        
        dynamic_subi_af{1,t}=[dynamic_subi_af{1,t},[nb_scenarios_rhs(t)*size_u(t)+1:nb_scenarios_rhs(t)*size_u(t)+size_x(t-1)]];
        dynamic_subj_af{1,t}=[dynamic_subj_af{1,t},[nb_scenarios_rhs(t)*(dimC+1)+size_x(t-1)+1:nb_scenarios_rhs(t)*(dimC+1)+2*size_x(t-1)]];
        dynamic_valij_af{1,t}=[dynamic_valij_af{1,t},-ones(1,size_x(t-1))];
        
        dynamic_subi_af{1,t}=[dynamic_subi_af{1,t},nb_scenarios_rhs(t)*size_u(t)+size_x(t-1)+i];
        dynamic_subj_af{1,t}=[dynamic_subj_af{1,t},nb_scenarios_rhs(t)*dimC+i];
        dynamic_valij_af{1,t}=[dynamic_valij_af{1,t},1];
    end
    dynamic_buc{1,t}=[repmat(costd{1,t},nb_scenarios_rhs(t),1);zeros(size_x(t-1),1);big_m*ones(nb_scenarios_rhs(t),1)];
    dynamic_blc{1,t}=[repmat(costd{1,t},nb_scenarios_rhs(t),1);zeros(size_x(t-1),1);-inf*ones(nb_scenarios_rhs(t),1)];
    
    dynamic_blx{1,t}=[];
    dynamic_bux{1,t}=[];
    
    dynamic_c{1,t}=[];
    
    for i=1:nb_scenarios_rhs(t)
        dynamic_blx{1,t}=[dynamic_blx{1,t};lower_pi{1,t}];
        dynamic_bux{1,t}=[dynamic_bux{1,t};upper_pi{1,t}];
        dynamic_c{1,t}=[dynamic_c{1,t};probabilities{1,t-1}(i)*bs{1,t}{1,i}];
        for j=1:size_x(t)
            if (upper_x{1,t}(j)==inf)
                dynamic_c{1,t}=[dynamic_c{1,t};0];
                dynamic_blx{1,t}=[dynamic_blx{1,t};0];
                dynamic_bux{1,t}=[dynamic_bux{1,t};0];
            else
                dynamic_c{1,t}=[dynamic_c{1,t};probabilities{1,t-1}(i)*upper_x{1,t}(j)];
                dynamic_blx{1,t}=[dynamic_blx{1,t};-inf];
                dynamic_bux{1,t}=[dynamic_bux{1,t};0];
            end
        end
        
        for j=1:size_x(t)
            if (lower_x{1,t}(j)==-inf)
                dynamic_c{1,t}=[dynamic_c{1,t};0];
                dynamic_blx{1,t}=[dynamic_blx{1,t};0];
                dynamic_bux{1,t}=[dynamic_bux{1,t};0];
            else
                dynamic_c{1,t}=[dynamic_c{1,t};probabilities{1,t-1}(i)*lower_x{1,t}(j)];
                dynamic_blx{1,t}=[dynamic_blx{1,t};0];
                dynamic_bux{1,t}=[dynamic_bux{1,t};inf];
            end
        end
        
        for j=1:size_u(t)
            if (upper_u{1,t}(j)==inf)
                dynamic_c{1,t}=[dynamic_c{1,t};0];
                dynamic_blx{1,t}=[dynamic_blx{1,t};0];
                dynamic_bux{1,t}=[dynamic_bux{1,t};0];
            else
                dynamic_c{1,t}=[dynamic_c{1,t};probabilities{1,t-1}(i)*upper_u{1,t}(j)];
                dynamic_blx{1,t}=[dynamic_blx{1,t};-inf];
                dynamic_bux{1,t}=[dynamic_bux{1,t};0];
            end
        end
        
        for j=1:size_u(t)
            if (lower_u{1,t}(j)==-inf)
                dynamic_c{1,t}=[dynamic_c{1,t};0];
                dynamic_blx{1,t}=[dynamic_blx{1,t};0];
                dynamic_bux{1,t}=[dynamic_bux{1,t};0];
            else
                dynamic_c{1,t}=[dynamic_c{1,t};probabilities{1,t-1}(i)*lower_u{1,t}(j)];
                dynamic_blx{1,t}=[dynamic_blx{1,t};0];
                dynamic_bux{1,t}=[dynamic_bux{1,t};inf];
            end
        end
    end
    for i=1:nb_scenarios_rhs(t)
        dynamic_c{1,t}=[dynamic_c{1,t};probabilities{1,t-1}(i)];
    end
    dynamic_blx{1,t}=[dynamic_blx{1,t};-inf*ones(nb_scenarios_rhs(t),1);zeros(2*size_x(t-1),1)];
    dynamic_bux{1,t}=[dynamic_bux{1,t};inf*ones(nb_scenarios_rhs(t),1);inf*ones(2*size_x(t-1),1)];
    dynamic_c{1,t}=[dynamic_c{1,t};-penaltiesQ(t-1,1)*ones(size_x(t-1),1);-penaltiesR(t-1,1)*ones(size_x(t-1),1)];
end

Cum_Probas=cell(1,T-1);
for t=1:T-1
    Cum_Probas{1,t}=[0,cumsum(probabilities{1,t})];
end

iter=1;

while (iter<=nb_iter_max)
    iter
    pis=cell(1,T-1);
    nus1=cell(1,T-1);
    nus2=cell(1,T-1);
    tic
    clear prob; 
    prob.c=dynamic_c{1,1};
    prob.blx=dynamic_blx{1,1};
    prob.bux=dynamic_bux{1,1};
    prob.blc=dynamic_blc{1,1};
    prob.buc=dynamic_buc{1,1};  
    prob.a = sparse(dynamic_subi_af{1,1},dynamic_subj_af{1,1},dynamic_valij_af{1,1},size_u(1)+iter,size_b(1)+1+2*(size_x(1)+size_u(1)));
    [~,res]=mosekopt('maximize echo(0)',prob);
    sol=res.sol.bas.xx;
    solsta=strcat('MSK_SOL_STA_', res.sol.bas.solsta);
    if (strcmp(solsta,'MSK_SOL_STA_PRIMAL_INFEASIBLE_CER')==1)
        disp('Unfeasible primal problem');
        pause
    elseif (strcmp(solsta,'MSK_SOL_STA_DUAL_INFEASIBLE_CER')==1)
        disp('Primal infinite optimal value');
        pause
    else
        pis{1,1}=sol(1:size_b(1));
        nus1{1,1}=sol(size_b(1)+1:size_b(1)+size_x(1));
        nus2{1,1}=sol(size_b(1)+size_x(1)+1:size_b(1)+2*size_x(1));
    end
    
    for t=2:T-1
         clear prob;
         dimCs=nb_scenarios_rhs(t)*(size_b(t)+2*(size_x(t)+size_u(t))+1);
         dynamic_c{1,t}(dimCs+1:dimCs+size_x(t-1))=-penaltiesQ(t-1,iter)*ones(size_x(t-1),1);
         dynamic_c{1,t}(dimCs+size_x(t-1)+1:dimCs+2*size_x(t-1))=-penaltiesR(t-1,iter)*ones(size_x(t-1),1);
         prob.c=dynamic_c{1,t};
         prob.blx=dynamic_blx{1,t};
         prob.bux=dynamic_bux{1,t};
         aux=costc{1,t-1}-sparse(a_subj{1,t-1},a_subi{1,t-1},a_valij{1,t-1},size_x(t-1),size_b(t-1))*pis{1,t-1}-nus1{1,t-1}-nus2{1,t-1};
         dynamic_blc{1,t}(nb_scenarios_rhs(t)*size_u(t)+1:nb_scenarios_rhs(t)*size_u(t)+size_x(t-1))=aux;
         dynamic_buc{1,t}(nb_scenarios_rhs(t)*size_u(t)+1:nb_scenarios_rhs(t)*size_u(t)+size_x(t-1))=aux;
         prob.blc=dynamic_blc{1,t};
         prob.buc=dynamic_buc{1,t};
         prob.a =sparse(dynamic_subi_af{1,t},dynamic_subj_af{1,t},dynamic_valij_af{1,t},nb_scenarios_rhs(t)*(iter+size_u(t))+size_x(t-1),nb_scenarios_rhs(t)*(1+size_b(t)+2*(size_x(t)+size_u(t)))+2*(size_x(t-1)));
         [~,res]=mosekopt('maximize echo(0)',prob);
        sol=res.sol.bas.xx;
        solsta=strcat('MSK_SOL_STA_', res.sol.bas.solsta);
        if (strcmp(solsta,'MSK_SOL_STA_PRIMAL_INFEASIBLE_CER')==1)
            disp('Unfeasible primal problem');
            pause
        elseif (strcmp(solsta,'MSK_SOL_STA_DUAL_INFEASIBLE_CER')==1)
            disp('Primal infinite optimal value');
            pause
        else
            Alea_Uniform=rand;
            [~,Index] = histc(Alea_Uniform,Cum_Probas{1,t-1});
            if (Alea_Uniform==1)
                Index=nb_scenarios_rhs(t);
            end
            dimC=size_b(t)+2*(size_x(t)+size_u(t));
            pis{1,t}=sol((Index-1)*dimC+1:(Index-1)*dimC+size_b(t));
            nus1{1,t}=sol((Index-1)*dimC+size_b(t)+1:(Index-1)*dimC+size_b(t)+size_x(t));
            nus2{1,t}=sol((Index-1)*dimC+size_b(t)+size_x(t)+1:(Index-1)*dimC+size_b(t)+2*size_x(t));
        end
    end
    
    for t=T:-1:2
        if (t==T)
            clear prob;
            dimC=size_b(T)+2*(size_x(T)+size_u(T));
            dimCs=nb_scenarios_rhs(T)*dimC;
            
            dynamic_c{1,T}(dimCs+1:dimCs+size_x(T-1))=-penaltiesQ(T-1,iter)*ones(size_x(T-1),1);
            dynamic_c{1,T}(dimCs+size_x(T-1)+1:dimCs+2*size_x(T-1))=-penaltiesR(T-1,iter)*ones(size_x(T-1),1);
            prob.c=dynamic_c{1,T};
            
            prob.blx=dynamic_blx{1,T};
            prob.bux=dynamic_bux{1,T};
            
            aux=costc{1,T-1}-sparse(a_subj{1,T-1},a_subi{1,T-1},a_valij{1,T-1},size_x(T-1),size_b(T-1))*pis{1,T-1}-nus1{1,T-1}-nus2{1,T-1};
            dynamic_blc{1,T}(nb_scenarios_rhs(T)*(size_u(T)+size_x(T))+1:nb_scenarios_rhs(T)*(size_u(T)+size_x(T))+size_x(T-1))=aux;
            dynamic_buc{1,T}(nb_scenarios_rhs(T)*(size_u(T)+size_x(T))+1:nb_scenarios_rhs(T)*(size_u(T)+size_x(T))+size_x(T-1))=aux;     
            prob.blc=dynamic_blc{1,T};
            prob.buc=dynamic_buc{1,T};
            
            prob.a = sparse(dynamic_subi_af{1,T},dynamic_subj_af{1,T},dynamic_valij_af{1,T},nb_scenarios_rhs(T)*(size_x(T)+size_u(T))+size_x(T-1),dimCs+2*size_x(T-1));
            [~,res]=mosekopt('maximize echo(0)',prob);
            sol=res.sol.bas.xx;
            solsta=strcat('MSK_SOL_STA_', res.sol.bas.solsta);
            if (strcmp(solsta,'MSK_SOL_STA_PRIMAL_INFEASIBLE_CER')==1)
                disp('Unfeasible primal problem');
                pause
            elseif (strcmp(solsta,'MSK_SOL_STA_DUAL_INFEASIBLE_CER')==1)
                disp('Primal infinite optimal value');
                pause
            else
                Errors{1,T-1}=[Errors{1,T-1};norm(sol(dimCs+1:dimCs+2*size_x(T-1)),1)];
                dual1=res.sol.bas.slc;
                dual2=res.sol.bas.suc;
                dual3=res.sol.bas.slx;
                dual4=res.sol.bas.sux;
                daux=dual1(nb_scenarios_rhs(T)*(size_x(T)+size_u(T))+1:nb_scenarios_rhs(T)*(size_x(T)+size_u(T))+size_x(T-1))-dual2(nb_scenarios_rhs(T)*(size_x(T)+size_u(T))+1:nb_scenarios_rhs(T)*(size_x(T)+size_u(T))+size_x(T-1));
                slopepi=-sparse(a_subi{1,T-1},a_subj{1,T-1},a_valij{1,T-1},size_b(T-1),size_x(T-1))*daux;
                slopenu1=-daux;
                slopenu2=-daux;
                intercept=prob.buc'*(dual1-dual2)-slopepi'*pis{1,T-1}-slopenu1'*nus1{1,T-1}-slopenu2'*nus2{1,T-1};
                for i=1:nb_scenarios_rhs(T)
                    intercept=intercept+dual3((i-1)*dimC+1:(i-1)*dimC+size_b(T))'*prob.blx((i-1)*dimC+1:(i-1)*dimC+size_b(T))-dual4((i-1)*dimC+1:(i-1)*dimC+size_b(T))'*prob.bux((i-1)*dimC+1:(i-1)*dimC+size_b(T));
                end
            end
        else
            clear prob;
            dimC=size_b(t)+2*(size_x(t)+size_u(t));
            dimCs=nb_scenarios_rhs(t)*(dimC+1);
            dynamic_c{1,t}(dimCs+1:dimCs+size_x(t-1))=-penaltiesQ(t-1,iter)*ones(size_x(t-1),1);
            dynamic_c{1,t}(dimCs+size_x(t-1)+1:dimCs+2*size_x(t-1))=-penaltiesR(t-1,iter)*ones(size_x(t-1),1);
            prob.c=dynamic_c{1,t};
            prob.blx=dynamic_blx{1,t};
            prob.bux=dynamic_bux{1,t};
            aux=costc{1,t-1}-sparse(a_subj{1,t-1},a_subi{1,t-1},a_valij{1,t-1},size_x(t-1),size_b(t-1))*pis{1,t-1}-nus1{1,t-1}-nus2{1,t-1};
            dynamic_blc{1,t}(nb_scenarios_rhs(t)*size_u(t)+1:nb_scenarios_rhs(t)*size_u(t)+size_x(t-1))=aux;
            dynamic_buc{1,t}(nb_scenarios_rhs(t)*size_u(t)+1:nb_scenarios_rhs(t)*size_u(t)+size_x(t-1))=aux;
            prob.blc=dynamic_blc{1,t};
            prob.buc=dynamic_buc{1,t};
            prob.a = sparse(dynamic_subi_af{1,t},dynamic_subj_af{1,t},dynamic_valij_af{1,t},nb_scenarios_rhs(t)*(iter+size_u(t)+1)+size_x(t-1),nb_scenarios_rhs(t)*(1+size_b(t)+2*(size_x(t)+size_u(t)))+2*size_x(t-1));            
            [~,res]=mosekopt('maximize echo(0)',prob);
            sol=res.sol.bas.xx;
            solsta=strcat('MSK_SOL_STA_', res.sol.bas.solsta);
            if (strcmp(solsta,'MSK_SOL_STA_PRIMAL_INFEASIBLE_CER')==1)
                disp('Unfeasible primal problem');
                pause
            elseif (strcmp(solsta,'MSK_SOL_STA_DUAL_INFEASIBLE_CER')==1)
                disp('Primal infinite optimal value');
                pause
            else
                Errors{1,t-1}=[Errors{1,t-1};norm(sol(dimCs+1:dimCs+2*size_x(t-1)),1)];
                dual1=res.sol.bas.slc;
                dual2=res.sol.bas.suc;
                dual3=res.sol.bas.slx;
                dual4=res.sol.bas.sux;
                daux=dual1(nb_scenarios_rhs(t)*size_u(t)+1:nb_scenarios_rhs(t)*size_u(t)+size_x(t-1))-dual2(nb_scenarios_rhs(t)*size_u(t)+1:nb_scenarios_rhs(t)*size_u(t)+size_x(t-1));
                slopepi=-sparse(a_subi{1,t-1},a_subj{1,t-1},a_valij{1,t-1},size_b(t-1),size_x(t-1))*daux;
                slopenu1=-daux;
                slopenu2=-daux;
                intercept=prob.buc(1:nb_scenarios_rhs(t)*size_u(t)+size_x(t-1))'*(dual1(1:nb_scenarios_rhs(t)*size_u(t)+size_x(t-1))-dual2(1:nb_scenarios_rhs(t)*size_u(t)+size_x(t-1)))-slopepi'*pis{1,t-1}-slopenu1'*nus1{1,t-1}-slopenu2'*nus2{1,t-1};
                intercept=intercept-prob.buc(nb_scenarios_rhs(t)*size_u(t)+size_x(t-1)+1:nb_scenarios_rhs(t)*size_u(t)+size_x(t-1)+nb_scenarios_rhs(t)*(iter+1))'*dual2(nb_scenarios_rhs(t)*size_u(t)+size_x(t-1)+1:nb_scenarios_rhs(t)*size_u(t)+size_x(t-1)+nb_scenarios_rhs(t)*(iter+1));
                for i=1:nb_scenarios_rhs(t)
                    intercept=intercept+prob.blx((i-1)*dimC+1:(i-1)*dimC+size_b(t))'*dual3((i-1)*dimC+1:(i-1)*dimC+size_b(t))-prob.bux((i-1)*dimC+1:(i-1)*dimC+size_b(t))'*dual4((i-1)*dimC+1:(i-1)*dimC+size_b(t));
                end
            end
        end
        if (t==2)
            row=size_u(1)+iter+1;
            dynamic_subi_af{1,1}=[dynamic_subi_af{1,1},row*ones(1,size_b(1)),row*ones(1,size_x(1)),row*ones(1,size_x(1)),row];
            dynamic_subj_af{1,1}=[dynamic_subj_af{1,1},[1:size_b(1)],[size_b(1)+1:size_b(1)+size_x(1)],[size_b(1)+size_x(1)+1:size_b(1)+2*size_x(1)],size_b(1)+2*(size_x(1)+size_u(1))+1];
            dynamic_valij_af{1,1}=[dynamic_valij_af{1,1},-slopepi',-slopenu1',-slopenu2',1];
            dynamic_buc{1,1}=[dynamic_buc{1,1};intercept];
            dynamic_blc{1,1}=[dynamic_blc{1,1};-inf];
        else
            dimC=size_b(t-1)+2*(size_x(t-1)+size_u(t-1));
            dynamic_buc{1,t-1}=[dynamic_buc{1,t-1};intercept*ones(nb_scenarios_rhs(t-1),1)];
            dynamic_blc{1,t-1}=[dynamic_blc{1,t-1};-inf*ones(nb_scenarios_rhs(t-1),1)];
            for i=1:nb_scenarios_rhs(t-1)
                row=size_x(t-1)+nb_scenarios_rhs(t-1)*size_u(t-1)+nb_scenarios_rhs(t-1)*iter+i;
                dynamic_subi_af{1,t-1}=[dynamic_subi_af{1,t-1},row*ones(1,size_b(t-1)),row*ones(1,size_x(t-1)),row*ones(1,size_x(t-1)),row];
                dynamic_subj_af{1,t-1}=[dynamic_subj_af{1,t-1},[(i-1)*dimC+1:(i-1)*dimC+size_b(t-1)],[(i-1)*dimC+size_b(t-1)+1:(i-1)*dimC+size_b(t-1)+size_x(t-1)],[(i-1)*dimC+size_b(t-1)+size_x(t-1)+1:(i-1)*dimC+size_b(t-1)+2*size_x(t-1)],nb_scenarios_rhs(t-1)*dimC+i];
                dynamic_valij_af{1,t-1}=[dynamic_valij_af{1,t-1},-slopepi',-slopenu1',-slopenu2',1];
            end
        end
    end
    
    clear prob; 
    prob.c=dynamic_c{1,1};
    prob.blx=dynamic_blx{1,1};
    prob.bux=dynamic_bux{1,1};
    prob.blc=dynamic_blc{1,1};
    prob.buc=dynamic_buc{1,1};
    prob.a = sparse(dynamic_subi_af{1,1},dynamic_subj_af{1,1},dynamic_valij_af{1,1},size_u(1)+iter+1,size_b(1)+1+2*(size_x(1)+size_u(1)));
    [~,res]=mosekopt('maximize echo(0)',prob);
    sol=res.sol.bas.xx;
    solsta=strcat('MSK_SOL_STA_', res.sol.bas.solsta);
    if (strcmp(solsta,'MSK_SOL_STA_PRIMAL_INFEASIBLE_CER')==1)
        disp('Unfeasible primal problem');
        pause
    elseif (strcmp(solsta,'MSK_SOL_STA_DUAL_INFEASIBLE_CER')==1)
        disp('Primal infinite optimal value');
        pause
    else
        ub=sol'*prob.c;
        upper_bounds=[upper_bounds;ub];
    end 
    time=[time;toc];
    iter=iter+1;
end

