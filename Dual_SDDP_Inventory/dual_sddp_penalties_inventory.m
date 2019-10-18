 
function [upper_bounds,time,Errors]=dual_sddp_penalties_inventory(Cost_C,Init_Stock,Cost_b,Cost_h,Probabilities,Demand_First_Stage,Demand_Scenarios,T,M,nb_iter_max,penaltiesQ,penaltiesR,large_bound,big_m)

upper_bounds=[];
time=[];

Errors=cell(1,T-1);

dynamic_subi_a=cell(1,T);
dynamic_subj_a=cell(1,T);
dynamic_valij_a=cell(1,T);

dynamic_buc=cell(1,T);

dynamic_subi_a{1,1}=[1,1,1,1,2];
dynamic_subj_a{1,1}=[1,2,3,4,5];
dynamic_valij_a{1,1}=[1,-1,1,-1,1];
dynamic_buc{1,1}=[Cost_C(1);big_m];

for t=2:T-1
    dynamic_subi_a{1,t}=[];
    dynamic_subj_a{1,t}=[];
    dynamic_valij_a{1,t}=[];
    for i=1:M
        dynamic_subi_a{1,t}=[dynamic_subi_a{1,t},i,i,i,i];
        dynamic_subj_a{1,t}=[dynamic_subj_a{1,t},4*(i-1)+1,4*(i-1)+2,4*(i-1)+3,4*i];
        dynamic_valij_a{1,t}=[dynamic_valij_a{1,t},1,-1,1,-1];
    end
    for i=1:M
        dynamic_subi_a{1,t}=[dynamic_subi_a{1,t},M+1];
        dynamic_subj_a{1,t}=[dynamic_subj_a{1,t},4*i];
        dynamic_valij_a{1,t}=[dynamic_valij_a{1,t},Probabilities(t-1,i)];
    end
    dynamic_subi_a{1,t}=[dynamic_subi_a{1,t},M+1,M+1];
    dynamic_subj_a{1,t}=[dynamic_subj_a{1,t},5*M+1,5*M+2];
    dynamic_valij_a{1,t}=[dynamic_valij_a{1,t},1,-1];
    for i=1:M
        dynamic_subi_a{1,t}=[dynamic_subi_a{1,t},M+i+1];
        dynamic_subj_a{1,t}=[dynamic_subj_a{1,t},4*M+i];
        dynamic_valij_a{1,t}=[dynamic_valij_a{1,t},1];
    end
    dynamic_buc{1,t}=[Cost_C(t)*ones(M,1);0;big_m*ones(M,1)];
end

dynamic_subi_a{1,T}=[];
dynamic_subj_a{1,T}=[];
dynamic_valij_a{1,T}=[];
dynamic_buc{1,T}=[];
for i=1:M
    dynamic_subi_a{1,T}=[dynamic_subi_a{1,T},i,i,i];
    dynamic_subj_a{1,T}=[dynamic_subj_a{1,T},3*(i-1)+1,3*(i-1)+2,3*i];
    dynamic_valij_a{1,T}=[dynamic_valij_a{1,T},-1,1,-1];
    dynamic_buc{1,T}=[dynamic_buc{1,T};Cost_C(T)];
end
for i=1:M
    dynamic_subi_a{1,T}=[dynamic_subi_a{1,T},M+1];
    dynamic_subj_a{1,T}=[dynamic_subj_a{1,T},3*i];
    dynamic_valij_a{1,T}=[dynamic_valij_a{1,T},Probabilities(T-1,i)];
end
dynamic_subi_a{1,T}=[dynamic_subi_a{1,T},M+1,M+1];
dynamic_subj_a{1,T}=[dynamic_subj_a{1,T},3*M+1,3*M+2];
dynamic_valij_a{1,T}=[dynamic_valij_a{1,T},1,-1];


Cum_Probas=cell(1,T-1);
for t=1:T-1
    Cum_Probas{1,t}=[0,cumsum(Probabilities(t,:))];
end

iter=1;

while (iter<=nb_iter_max)
    iter
    pis=cell(1,T-1);
    tic
    clear prob;
    prob.c=[Demand_First_Stage;-Demand_First_Stage;Demand_First_Stage;-Init_Stock;1];
    prob.blx=[-large_bound;-Cost_b(1);-Cost_h(1);-inf;-inf];
    %prob.blx=[-large_bound;-Cost_b(1);-Cost_h(1);-large_bound;-inf];
    prob.bux=[large_bound;0;0;0;inf];
    %prob.blc=[Cost_C(1);-large_bound*ones(iter,1)];
    prob.blc=[Cost_C(1);-inf*ones(iter,1)];
    prob.buc=dynamic_buc{1,1};
    prob.a = sparse(dynamic_subi_a{1,1},dynamic_subj_a{1,1},dynamic_valij_a{1,1},1+iter,5);
    [~,res]=mosekopt('maximize echo(0)',prob);
    sol=res.sol.bas.xx;
%     sol(1)
%     sol(2)
%     sol(3)
%     sol(4)
%     sol(5)
%     sol'*prob.c - 23.660254
    solsta=strcat('MSK_SOL_STA_', res.sol.bas.solsta);
    if (strcmp(solsta,'MSK_SOL_STA_PRIMAL_INFEASIBLE_CER')==1)
        disp('Unfeasible primal problem');
        pause
    elseif (strcmp(solsta,'MSK_SOL_STA_DUAL_INFEASIBLE_CER')==1)
        disp('Primal infinite optimal value');
        pause
    else
        pis{1,1}=sol(1);
    end
    
    for t=2:T-1
        clear prob;
        aux=[];
        for j=1:M
            aux=[aux;Probabilities(t-1,j)*Demand_Scenarios(t-1,j);-Probabilities(t-1,j)*Demand_Scenarios(t-1,j);Probabilities(t-1,j)*Demand_Scenarios(t-1,j);0];
        end
        for j=1:M
            aux=[aux;Probabilities(t-1,j)];
        end
        prob.c=[aux;-penaltiesQ(t-1,iter);-penaltiesR(t-1,iter)];
        Lower=[];
        Upper=[];
        for i=1:M
            %Lower=[Lower;-large_bound;-Cost_b(t);-Cost_h(t);-large_bound];
            Lower=[Lower;-large_bound;-Cost_b(t);-Cost_h(t);-inf];
            Upper=[Upper;large_bound;zeros(3,1)];
        end
        Lower=[Lower;-inf*ones(M,1);0;0];
        Upper=[Upper;inf*ones(M+2,1)];
        prob.blx=Lower;
        prob.bux=Upper;
        dynamic_buc{1,t}(M+1)=-Cost_C(t)+pis{1,t-1};
        prob.blc=[dynamic_buc{1,t}(1:M+1);-inf*ones(M*iter,1)];
        prob.buc=dynamic_buc{1,t};
        prob.a = sparse(dynamic_subi_a{1,t},dynamic_subj_a{1,t},dynamic_valij_a{1,t},M+1+M*iter,5*M+2);
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
                Index=M;
            end
            pis{1,t}=sol((Index-1)*4+1);
        end
    end
    
    for t=T:-1:2
        if (t==T)
            clear prob;
            aux=[];
            for i=1:M
                aux=[aux;-Probabilities(T-1,i)*Demand_Scenarios(T-1,i);Probabilities(T-1,i)*Demand_Scenarios(T-1,i);0];
            end
            prob.c=[aux;-penaltiesQ(T-1,iter);-penaltiesR(T-1,iter)];
            Lower=[];
            for i=1:M
                Lower=[Lower;-Cost_b(T);-Cost_h(T);-inf];
            end
            Lower=[Lower;0;0];
            prob.blx=Lower;
            prob.bux=[zeros(3*M,1);inf;inf];
            prob.blc=[dynamic_buc{1,T};-Cost_C(T)+pis{1,T-1}];
            prob.buc=prob.blc;
            prob.a = sparse(dynamic_subi_a{1,T},dynamic_subj_a{1,T},dynamic_valij_a{1,T},M+1,3*M+2);
            [~,res]=mosekopt('maximize echo(0)',prob);
            sol=res.sol.bas.xx;
            
            if t == 2
                sol
                sol'*prob.c - 23.660254
            end
            
            
            solsta=strcat('MSK_SOL_STA_', res.sol.bas.solsta);
            if (strcmp(solsta,'MSK_SOL_STA_PRIMAL_INFEASIBLE_CER')==1)
                disp('Unfeasible primal problem');
                pause
            elseif (strcmp(solsta,'MSK_SOL_STA_DUAL_INFEASIBLE_CER')==1)
                disp('Primal infinite optimal value');
                pause
            else
                Errors{1,T-1}=[Errors{1,T-1};abs(sol(3*M+1))+abs(sol(3*M+2))];
                dual1=res.sol.bas.slc;
                dual2=res.sol.bas.suc;
                dual3=res.sol.bas.slx;
                dual4=res.sol.bas.sux;
                interceptD=Cost_C(T)*sum(dual1(1:M)-dual2(1:M))-Cost_C(T)*(dual1(M+1)-dual2(M+1));
                for k=1:M
                    interceptD=interceptD-Cost_b(T)*dual3(1+3*(k-1))-Cost_h(T)*dual3(2+3*(k-1));
                end
                slope=dual1(M+1)-dual2(M+1);
               
                
            end
        else
            clear prob;
            aux=[];
            for j=1:M
                aux=[aux;Probabilities(t-1,j)*Demand_Scenarios(t-1,j);-Probabilities(t-1,j)*Demand_Scenarios(t-1,j);Probabilities(t-1,j)*Demand_Scenarios(t-1,j);0];
            end
            for j=1:M
                aux=[aux;Probabilities(t-1,j)];
            end
            prob.c=[aux;-penaltiesQ(t-1,iter);-penaltiesR(t-1,iter)];
            Lower=[];
            Upper=[];
            for i=1:M
                %Lower=[Lower;-large_bound;-Cost_b(t);-Cost_h(t);-large_bound];
                Lower=[Lower;-large_bound;-Cost_b(t);-Cost_h(t);-inf];
                Upper=[Upper;large_bound;zeros(3,1)];
            end
            %Bounds on Psi
            Lower=[Lower;-inf*ones(M,1);0;0];
            Upper=[Upper;inf*ones(M+2,1)];
            prob.blx=Lower;
            prob.bux=Upper;
            dynamic_buc{1,t}(M+1)=-Cost_C(t)+pis{1,t-1};
            prob.blc=[dynamic_buc{1,t}(1:M+1);-inf*ones(M*(iter+1),1)];
            prob.buc=dynamic_buc{1,t};
            prob.a = sparse(dynamic_subi_a{1,t},dynamic_subj_a{1,t},dynamic_valij_a{1,t},M+1+M*(iter+1),5*M+2);
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
                Errors{1,t-1}=[Errors{1,t-1};abs(sol(5*M+1))+abs(sol(5*M+2))];
                dual1=res.sol.bas.slc;
                dual2=res.sol.bas.suc;
                dual3=res.sol.bas.slx;
                dual4=res.sol.bas.sux;
                interceptD=Cost_C(t)*sum(dual1(1:M)-dual2(1:M))-Cost_C(t)*(dual1(M+1)-dual2(M+1))-(prob.buc(M+2:M+1+M*(iter+1)))'*(dual2(M+2:M+1+M*(iter+1)));
                for k=1:M
                    interceptD=interceptD-Cost_b(t)*dual3(2+4*(k-1))-Cost_h(t)*dual3(3+4*(k-1))-Upper(4*(k-1)+1)*dual4(4*(k-1)+1)+Lower(4*(k-1)+1)*dual3(4*(k-1)+1);
                end
                slope=dual1(M+1)-dual2(M+1);
                %pause
                
            end
        end
        if (t==2)
            dynamic_subi_a{1,1}=[dynamic_subi_a{1,1},iter+2,iter+2];
            dynamic_subj_a{1,1}=[dynamic_subj_a{1,1},1,5];
            dynamic_valij_a{1,1}=[dynamic_valij_a{1,1},-slope,1];
            dynamic_buc{1,1}=[dynamic_buc{1,1};interceptD];
            slope
            interceptD
        else
            dynamic_buc{1,t-1}=[dynamic_buc{1,t-1};interceptD*ones(M,1)];
            for i=1:M
                dynamic_subi_a{1,t-1}=[dynamic_subi_a{1,t-1},M+1+M*iter+i,M+1+M*iter+i];
                dynamic_subj_a{1,t-1}=[dynamic_subj_a{1,t-1},4*(i-1)+1,4*M+i];
                dynamic_valij_a{1,t-1}=[dynamic_valij_a{1,t-1},-slope,1];
            end
        end
    end
    
    clear prob;
    prob.c=[Demand_First_Stage;-Demand_First_Stage;Demand_First_Stage;-Init_Stock;1];
    prob.blx=[-large_bound;-Cost_b(1);-Cost_h(1);-large_bound;-inf];
    prob.bux=[large_bound;0;0;0;inf];
    prob.blc=[Cost_C(1);-inf*ones(iter+1,1)];
    prob.buc=dynamic_buc{1,1};
    prob.a = sparse(dynamic_subi_a{1,1},dynamic_subj_a{1,1},dynamic_valij_a{1,1},2+iter,5);
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

upper_bounds=upper_bounds-Cost_C(1)*Init_Stock;
