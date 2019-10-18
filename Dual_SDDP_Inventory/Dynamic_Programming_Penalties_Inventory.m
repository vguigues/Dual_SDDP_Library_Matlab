

function [time,optimal_values,Bellman_Functions,bool_prfeas]=Dynamic_Programming_Penalties_Inventory(Cost_C,Init_Stock,Cost_b,Cost_h,Probabilities,Demand_First_Stage,Demand_Scenarios,T,M,large_bound,xmin,xmax,Nbpoints,gammas)

L=length(gammas);
Bellman_Functions=cell(T-1,L);
optimal_values=-inf*ones(1,L);

time=[];


for index_gamma=1:L
    
    dynamic_subi=cell(1,T);
    dynamic_subj=cell(1,T);
    dynamic_valij=cell(1,T);
    dynamic_buc=cell(1,T);
    
    dynamic_buc{1,T}=[Cost_C(T)*ones(M,1);0];
    dynamic_subi{1,T}=[];
    dynamic_subj{1,T}=[];
    dynamic_valij{1,T}=[];
    dynamic_c=cell(1,T);
    
    dynamic_c{1,T}=[];
    for i=1:M
        dynamic_subi{1,T}=[dynamic_subi{1,T},i,i,i];
        dynamic_subj{1,T}=[dynamic_subj{1,T},3*(i-1)+1,3*(i-1)+2,3*i];
        dynamic_valij{1,T}=[dynamic_valij{1,T},-1,1,-1];
        dynamic_c{1,T}=[dynamic_c{1,T};-Probabilities(T-1,i)*Demand_Scenarios(T-1,i);Probabilities(T-1,i)*Demand_Scenarios(T-1,i);0];
    end
    dynamic_c{1,T}=[dynamic_c{1,T};-gammas(index_gamma);-gammas(index_gamma)];
    
    for i=1:M
        dynamic_subi{1,T}=[dynamic_subi{1,T},M+1];
        dynamic_subj{1,T}=[dynamic_subj{1,T},3*i];
        dynamic_valij{1,T}=[dynamic_valij{1,T},Probabilities(T-1,i)];
    end
    
    dynamic_subi{1,T}=[dynamic_subi{1,T},M+1,M+1];
    dynamic_subj{1,T}=[dynamic_subj{1,T},3*M+1,3*M+2];
    dynamic_valij{1,T}=[dynamic_valij{1,T},1,-1];
    
    for t=2:T-1
        dynamic_buc{1,t}=[];
        dynamic_subi{1,t}=[];
        dynamic_subj{1,t}=[];
        dynamic_valij{1,t}=[];
        for i=1:M
            dynamic_subi{1,t}=[dynamic_subi{1,t},i,i,i,i];
            dynamic_subj{1,t}=[dynamic_subj{1,t},4*(i-1)+1,4*(i-1)+2,4*(i-1)+3,4*i];
            dynamic_valij{1,t}=[dynamic_valij{1,t},1,-1,1,-1];
            dynamic_c{1,t}=[dynamic_c{1,t};Probabilities(t-1,i)*Demand_Scenarios(t-1,i);-Probabilities(t-1,i)*Demand_Scenarios(t-1,i);Probabilities(t-1,i)*Demand_Scenarios(t-1,i);0];
        end
        for i=1:M
            dynamic_subi{1,t}=[dynamic_subi{1,t},M+1];
            dynamic_subj{1,t}=[dynamic_subj{1,t},4*i];
            dynamic_valij{1,t}=[dynamic_valij{1,t},Probabilities(t-1,i)];
            dynamic_c{1,t}=[dynamic_c{1,t};Probabilities(t-1,i)];
        end
        dynamic_subi{1,t}=[dynamic_subi{1,t},M+1,M+1];
        dynamic_subj{1,t}=[dynamic_subj{1,t},5*M+1,5*M+2];
        dynamic_valij{1,t}=[dynamic_valij{1,t},1,-1];
        dynamic_buc{1,t}=[Cost_C(t)*ones(M,1);0];
        dynamic_c{1,t}=[dynamic_c{1,t};-gammas(index_gamma);-gammas(index_gamma)];
    end
    
    dynamic_subi{1,1}=[1,1,1,1];
    dynamic_subj{1,1}=[1,2,3,4];
    dynamic_valij{1,1}=[1,-1,1,-1];
    dynamic_buc{1,1}=[Cost_C(1)];
    dynamic_c{1,1}=[Demand_First_Stage;-Demand_First_Stage;Demand_First_Stage;-Init_Stock;1];
    
    bool_prfeas=1;
    bool_dualfeas=1;
    
    tic
    
    pas=(xmax(T-1)-xmin(T-1))/Nbpoints;
    abscisses=[xmin(T-1):pas:xmax(T-1)];
    index=1;
    clear prob;
    prob.c=dynamic_c{1,T};
    dynamic_buc{1,T}(M+1)=-Cost_C(T)+abscisses(index);
    prob.blc=dynamic_buc{1,T};
    prob.buc=prob.blc;
    prob.blx=[repmat([-Cost_b(T);-Cost_h(T);-inf],M,1);0;0];
    prob.bux=[zeros(3*M,1);inf;inf];
    prob.a = sparse(dynamic_subi{1,T},dynamic_subj{1,T},dynamic_valij{1,T},M+1,3*M+2);
    [r,res]=mosekopt('maximize echo(0)',prob);
    sol=res.sol.bas.xx;
    solsta=strcat('MSK_SOL_STA_', res.sol.bas.solsta);
    if (strcmp(solsta,'MSK_SOL_STA_PRIMAL_INFEASIBLE_CER')==1)
        disp('Unfeasible primal problem');
        bool_prfeas=0;
    elseif (strcmp(solsta,'MSK_SOL_STA_DUAL_INFEASIBLE_CER')==1)
        disp('Primal infinite optimal value');
        bool_dualfeas=0;
        optimal_values(index_gamma)=inf;
    end
    
    if ((bool_prfeas==1)&&(bool_dualfeas==1))
        
        prev_V=sol'*prob.c;
        Bellman_Functions{T-1,index_gamma}=[prev_V];
        nb_cuts=zeros(T-1,1);
        index=2;
        
        while ((index<=Nbpoints+1)&&(bool_prfeas==1)&&(bool_dualfeas==1))
            clear prob;
            prob.c=dynamic_c{1,T};
            dynamic_buc{1,T}(M+1)=-Cost_C(T)+abscisses(index);
            prob.blc=dynamic_buc{1,T};
            prob.buc=prob.blc;
            prob.blx=[repmat([-Cost_b(T);-Cost_h(T);-inf],M,1);0;0];
            prob.bux=[zeros(3*M,1);inf;inf];
            prob.a = sparse(dynamic_subi{1,T},dynamic_subj{1,T},dynamic_valij{1,T},M+1,3*M+2);
            [r,res]=mosekopt('maximize echo(0)',prob);
            sol=res.sol.bas.xx;
            solsta=strcat('MSK_SOL_STA_', res.sol.bas.solsta);
            if (strcmp(solsta,'MSK_SOL_STA_PRIMAL_INFEASIBLE_CER')==1)
                disp('Unfeasible primal problem');
                bool_prfeas=0;
            elseif (strcmp(solsta,'MSK_SOL_STA_DUAL_INFEASIBLE_CER')==1)
                disp('Primal infinite optimal value');
                bool_dualfeas=0;
                optimal_values(index_gamma)=inf;
            else
                current_V=sol'*prob.c;
                beta=(current_V-prev_V)/(abscisses(index)-abscisses(index-1));
                theta=prev_V-beta*abscisses(index-1);
                Bellman_Functions{T-1,index_gamma}=[Bellman_Functions{T-1,index_gamma};current_V];
                if (T==2)
                    dynamic_subi{1,1}=[dynamic_subi{1,1},2+nb_cuts(1),2+nb_cuts(1)];
                    dynamic_subj{1,1}=[dynamic_subj{1,1},1,5];
                    dynamic_valij{1,1}=[dynamic_valij{1,1},-beta,1];
                    dynamic_buc{1,1}=[dynamic_buc{1,1};theta];
                else
                    for i=1:M
                        dynamic_subi{1,T-1}=[dynamic_subi{1,T-1},M+1+nb_cuts(T-1)*M+i,M+1+nb_cuts(T-1)*M+i];
                        dynamic_subj{1,T-1}=[dynamic_subj{1,T-1},4*(i-1)+1,4*M+i];
                        dynamic_valij{1,T-1}=[dynamic_valij{1,T-1},-beta,1];
                    end
                    dynamic_buc{1,T-1}=[dynamic_buc{1,T-1};theta*ones(M,1)];
                end
                prev_V=current_V;
                nb_cuts(T-1)=nb_cuts(T-1)+1;
            end
            index=index+1;
        end
        
        t=T-1;
        
        while ((t>=2)&&(bool_prfeas==1)&&(bool_dualfeas==1))
            t
            pas=(xmax(t-1)-xmin(t-1))/Nbpoints;
            abscisses=[xmin(t-1):pas:xmax(t-1)];
            index=1;
            clear prob;
            prob.c=dynamic_c{1,t};
            dynamic_buc{1,t}(M+1)=-Cost_C(t)+abscisses(index);
            prob.blc=[dynamic_buc{1,t}(1:M+1);-inf*ones(nb_cuts(t)*M,1)];
            prob.buc=dynamic_buc{1,t};
            prob.blx=[repmat([-inf;-Cost_b(t);-Cost_h(t);-inf],M,1);-inf*ones(M,1);0;0];
            prob.bux=[repmat([inf;zeros(3,1)],M,1);inf*ones(M+2,1)];
            prob.a = sparse(dynamic_subi{1,t},dynamic_subj{1,t},dynamic_valij{1,t},M+1+nb_cuts(t)*M,5*M+2);
            [r,res]=mosekopt('maximize echo(0)',prob);
            sol=res.sol.bas.xx;
            solsta=strcat('MSK_SOL_STA_', res.sol.bas.solsta);
            if (strcmp(solsta,'MSK_SOL_STA_PRIMAL_INFEASIBLE_CER')==1)
                disp('Unfeasible primal problem');
                bool_prfeas=0;
            elseif (strcmp(solsta,'MSK_SOL_STA_DUAL_INFEASIBLE_CER')==1)
                disp('Primal infinite optimal value');
                bool_dualfeas=0;
                optimal_values(index_gamma)=inf;
            else 
                prev_V=sol'*prob.c;
                Bellman_Functions{t-1,index_gamma}=[prev_V];
                index=2;
                while ((index<=Nbpoints+1)&&(bool_prfeas==1)&&(bool_dualfeas==1))
                    clear prob;
                    prob.c=dynamic_c{1,t};
                    dynamic_buc{1,t}(M+1)=-Cost_C(t)+abscisses(index);
                    prob.blc=[dynamic_buc{1,t}(1:M+1);-inf*ones(nb_cuts(t)*M,1)];
                    prob.buc=dynamic_buc{1,t};
                    prob.blx=[repmat([-inf;-Cost_b(t);-Cost_h(t);-inf],M,1);-inf*ones(M,1);0;0];
                    prob.bux=[repmat([inf;zeros(3,1)],M,1);inf*ones(M+2,1)];
                    prob.a = sparse(dynamic_subi{1,t},dynamic_subj{1,t},dynamic_valij{1,t},M+1+nb_cuts(t)*M,5*M+2);
                    [r,res]=mosekopt('maximize echo(0)',prob);
                    sol=res.sol.bas.xx;
                    solsta=strcat('MSK_SOL_STA_',res.sol.bas.solsta);
                    if (strcmp(solsta,'MSK_SOL_STA_PRIMAL_INFEASIBLE_CER')==1)
                        disp('Unfeasible primal problem');
                        bool_prfeas=0;
                    elseif (strcmp(solsta,'MSK_SOL_STA_DUAL_INFEASIBLE_CER')==1)
                        disp('Primal infinite optimal value');
                        bool_dualfeas=0;
                        optimal_values(index_gamma)=inf;
                    else
                        current_V=sol'*prob.c;
                        beta=(current_V-prev_V)/(abscisses(index)-abscisses(index-1));
                        theta=prev_V-beta*abscisses(index-1);
                        Bellman_Functions{t-1,index_gamma}=[Bellman_Functions{t-1,index_gamma};current_V];
                        if (t==2)
                            dynamic_subi{1,1}=[dynamic_subi{1,1},2+nb_cuts(1),2+nb_cuts(1)];
                            dynamic_subj{1,1}=[dynamic_subj{1,1},1,5];
                            dynamic_valij{1,1}=[dynamic_valij{1,1},-beta,1];
                            dynamic_buc{1,1}=[dynamic_buc{1,1};theta];
                        else
                            for i=1:M
                                dynamic_subi{1,t-1}=[dynamic_subi{1,t-1},M+1+nb_cuts(t-1)*M+i,M+1+nb_cuts(t-1)*M+i];
                                dynamic_subj{1,t-1}=[dynamic_subj{1,t-1},4*(i-1)+1,4*M+i];
                                dynamic_valij{1,t-1}=[dynamic_valij{1,t-1},-beta,1];
                            end
                            dynamic_buc{1,t-1}=[dynamic_buc{1,t-1};theta*ones(M,1)];
                        end
                        prev_V=current_V;
                        nb_cuts(t-1)=nb_cuts(t-1)+1;
                    end
                    index=index+1;
                end
            end
            t=t-1;
        end
        clear prob;
        prob.c=dynamic_c{1,1};
        prob.buc=dynamic_buc{1,1};
        prob.blc=[Cost_C(1);-inf*ones(nb_cuts(1),1)];
        prob.blx=[-inf;-Cost_b(1);-Cost_h(1);-inf;-inf];
        prob.bux=[inf;zeros(3,1);inf];
        prob.a = sparse(dynamic_subi{1,1},dynamic_subj{1,1},dynamic_valij{1,1},1+nb_cuts(1),5);
        [r,res]=mosekopt('maximize echo(0)',prob);
        sol=res.sol.bas.xx;
        solsta=strcat('MSK_SOL_STA_', res.sol.bas.solsta);
        if (strcmp(solsta,'MSK_SOL_STA_PRIMAL_INFEASIBLE_CER')==1)
            disp('Unfeasible primal problem');
            bool_prfeas=0;
        elseif (strcmp(solsta,'MSK_SOL_STA_DUAL_INFEASIBLE_CER')==1)
            disp('Primal infinite optimal value');
            bool_dualfeas=0;
            optimal_values(index_gamma)=inf;
        else 
            optimal_values(index_gamma)=sol'*prob.c-Cost_C(1)*Init_Stock;
        end
    end
    time=[time;toc];
end



