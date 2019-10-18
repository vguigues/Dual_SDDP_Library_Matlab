
function [time,upper_bounds,lower_pi,upper_pi]=Primal_Dual_SDDP_Inventory_Management_Feasibility_Cuts(Cost_C,Init_Stock,Cost_b,Cost_h,Probabilities,tol,Demand_First_Stage,Demand_Scenarios,T,M,talpha,large_bound,big_m,iter_max)

upper_bounds=[];
time=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialisations for feasbility cuts
%for stages 2,3,...,T.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The problem data for checking feasibility at t=2,...,T, is stored
%in the structures for t-1, for instance dynamic_feas_cuts_subi{1,t-1}+

dynamic_feas_cuts_subi=cell(1,T-1);
dynamic_feas_cuts_subj=cell(1,T-1);
dynamic_feas_cuts_valij=cell(1,T-1);
dynamic_feas_cuts_buc=cell(1,T-1);

dynamic_feas_cuts_subi{1,T-1}=[];
dynamic_feas_cuts_subj{1,T-1}=[];
dynamic_feas_cuts_valij{1,T-1}=[];

for t=2:T
    dynamic_feas_cuts_buc{1,t-1}=Cost_C(t)*ones(M,1);
end

for i=1:M
    dynamic_feas_cuts_subi{1,T-1}=[dynamic_feas_cuts_subi{1,T-1},i,i,i];
    dynamic_feas_cuts_subj{1,T-1}=[dynamic_feas_cuts_subj{1,T-1},3*(i-1)+1,3*(i-1)+2,3*i];
    dynamic_feas_cuts_valij{1,T-1}=[dynamic_feas_cuts_valij{1,T-1},1,-1,1];
end
for i=1:M
    dynamic_feas_cuts_subi{1,T-1}=[dynamic_feas_cuts_subi{1,T-1},M+1];
    dynamic_feas_cuts_subj{1,T-1}=[dynamic_feas_cuts_subj{1,T-1},3*i];
    dynamic_feas_cuts_valij{1,T-1}=[dynamic_feas_cuts_valij{1,T-1},-Probabilities(T-1,i)];
end
dynamic_feas_cuts_subi{1,T-1}=[dynamic_feas_cuts_subi{1,T-1},M+1,M+1];
dynamic_feas_cuts_subj{1,T-1}=[dynamic_feas_cuts_subj{1,T-1},3*M+1,3*M+2];
dynamic_feas_cuts_valij{1,T-1}=[dynamic_feas_cuts_valij{1,T-1},1,-1];

for t=2:T-1
    for i=1:M
        dynamic_feas_cuts_subi{1,t-1}=[dynamic_feas_cuts_subi{1,t-1},i,i,i,i];
        dynamic_feas_cuts_subj{1,t-1}=[dynamic_feas_cuts_subj{1,t-1},4*(i-1)+1,4*(i-1)+2,4*(i-1)+3,4*i];
        dynamic_feas_cuts_valij{1,t-1}=[dynamic_feas_cuts_valij{1,t-1},1,1,-1,1];
    end
    for i=1:M
        dynamic_feas_cuts_subi{1,t-1}=[dynamic_feas_cuts_subi{1,t-1},M+1];
        dynamic_feas_cuts_subj{1,t-1}=[dynamic_feas_cuts_subj{1,t-1},4*i];
        dynamic_feas_cuts_valij{1,t-1}=[dynamic_feas_cuts_valij{1,t-1},-Probabilities(t-1,i)];
    end
    dynamic_feas_cuts_subi{1,t-1}=[dynamic_feas_cuts_subi{1,t-1},M+1,M+1];
    dynamic_feas_cuts_subj{1,t-1}=[dynamic_feas_cuts_subj{1,t-1},4*M+1,4*M+2];
    dynamic_feas_cuts_valij{1,t-1}=[dynamic_feas_cuts_valij{1,t-1},1,-1];
end

Cum_Probas=zeros(T-1,M+1);
for t=1:T-1
    Cum_Probas(t,:)=[0,cumsum(Probabilities(t,:))];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initializations for forward pass
%To check

upper_pi=large_bound*ones(T-1,1);
lower_pi=-large_bound*ones(T-1,1);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cum_Probas=cell(1,T-1);
for t=1:T-1
    Cum_Probas{1,t}=[0,cumsum(Probabilities(t,:))];
end

iter=1;

while (iter<=iter_max)
       iter
       t=1;
       pis=zeros(1,T-1);
       tic
       while (t<=T)
             if (t==1)
                clear prob;
                prob.c=[Demand_First_Stage;-Demand_First_Stage;Demand_First_Stage;-Init_Stock;1];
                prob.blx=[lower_pi(1);-Cost_b(1);-Cost_h(1);-large_bound;-inf];
                prob.bux=[upper_pi(1);0;0;0;inf];
                prob.blc=[Cost_C(1);-inf*ones(iter,1)];
                prob.buc=dynamic_buc{1,1};
                prob.a = sparse(dynamic_subi_a{1,1},dynamic_subj_a{1,1},dynamic_valij_a{1,1},1+iter,5);
                [~,res]=mosekopt('maximize echo(0)',prob);
                sol=res.sol.bas.xx;
                pis(1)=sol(1);
                t=t+1;
             elseif (t==T)
                    clear prob;
                    prob.c=[zeros(3*M,1);1;1];
                    Lower=[];
                    Upper=[];
                    for i=1:M
                        Lower=[Lower;0;0;0];
                        Upper=[Upper;Cost_b(T);Cost_h(T);inf];
                    end
                    prob.blx=[Lower;0;0];
                    prob.bux=[Upper;inf;inf];
                    prob.blc=[dynamic_feas_cuts_buc{1,T-1};-Cost_C(T)+pis(T-1)];
                    prob.buc=prob.blc;
                    prob.a = sparse(dynamic_feas_cuts_subi{1,T-1},dynamic_feas_cuts_subj{1,T-1},dynamic_feas_cuts_valij{1,T-1},1+M,3*M+2);
                    [~,res]=mosekopt('minimize echo(0)',prob);
                    sol=res.sol.bas.xx;
                    aux=sol(3*M+1)+sol(3*M+2);
                    
                    if (aux>10^(-6))
                        dual1=res.sol.bas.slc;
                        dual2=res.sol.bas.suc;
                        dual3=res.sol.bas.sux;
                        constant_coeff=0;
                        for i=1:M
                            constant_coeff=constant_coeff-Cost_b(T)*dual3(3*(i-1)+1)-Cost_h(T)*dual3(3*(i-1)+2);
                        end
                        constant_coeff=constant_coeff+Cost_C(T)*sum(dual1(1:M)-dual2(1:M))-Cost_C(T)*(dual1(M+1)-dual2(M+1));
                        slopef=dual1(M+1)-dual2(M+1);
                        if (slopef>0)
                            bound=-constant_coeff/slopef;
                            upper_pi(T-1)=min(upper_pi(T-1),bound);
                        elseif (slopef<0)
                            bound=-constant_coeff/slopef;
                            lower_pi(T-1)=max(lower_pi(T-1),bound);
                        else
                            disp('Null slope');
                            pause
                        end
                        t=t-1;
                    else
                        t=t+1;
                    end     
             else
                    clear prob;  
                    
%                       clear prob;
%         prob.c=[zeros(4*M,1);1;1];
%         prob.blc=[Cost_C(t)*ones(M,1);-Cost_C(t)+abscisses(index)];
%         prob.buc=prob.blc;
%         prob.blx=[repmat([lower_pi(t);-Cost_b(t);-Cost_h(t);-inf],M,1);0;0];
%         prob.bux=[repmat([upper_pi(t);zeros(3,1)],M,1);inf;inf];
%         prob.a=sparse(dynamic_feas_subi{1,t-1},dynamic_feas_subj{1,t-1},dynamic_feas_valij{1,t-1},M+1,4*M+2);
%         prob.blx
%         prob.bux
%         prob.buc
%         full(prob.a)
%         [r,res]=mosekopt('minimize echo(0)',prob);
%         sol=res.sol.bas.xx;
%         error=abs(sol(4*M+1))+abs(sol(4*M+2));
% %                     
                    
                    prob.c=[zeros(4*M,1);1;1];  
                    Lower=[];
                    Upper=[];
                    for i=1:M
                        Lower=[Lower;lower_pi(t);0;0;0];
                        Upper=[Upper;upper_pi(t);Cost_b(t);Cost_h(t);inf];
                    end
                    prob.blx=[Lower;0;0];
                    prob.bux=[Upper;inf;inf];
                    prob.blc=[dynamic_feas_cuts_buc{1,t-1};-Cost_C(t)+pis(t-1)];
                    prob.buc=prob.blc;
                    prob.a = sparse(dynamic_feas_cuts_subi{1,t-1},dynamic_feas_cuts_subj{1,t-1},dynamic_feas_cuts_valij{1,t-1},1+M,4*M+2);
                    [~,res]=mosekopt('minimize echo(0)',prob);
                    sol=res.sol.bas.xx;
                    %t 
                    %iter
                    %[prob.blx prob.bux]
                    %prob.buc
                    %full(prob.a)
                    %pis(t-1)
                    %sol(4*M+1)
                    %sol(4*M+2)
                    aux=sol(4*M+1)+sol(4*M+2);
                    %pause
                    if (aux>10^(-6))
                        dual1=res.sol.bas.slc;
                        dual2=res.sol.bas.suc;
                        dual3=res.sol.bas.slx;
                        dual4=res.sol.bas.sux;
                        constant_coeff=Cost_C(t)*sum(dual1(1:M)-dual2(1:M))-Cost_C(t)*(dual1(M+1)-dual2(M+1));
                        for i=1:M
                            constant_coeff=constant_coeff-Cost_b(t)*dual4(4*(i-1)+2)-Cost_h(t)*dual4(4*(i-1)+3);
                            constant_coeff=constant_coeff+lower_pi(t)*dual3(4*(i-1)+1)-upper_pi(t)*dual4(4*(i-1)+1);
                        end
                        slopef=dual1(M+1)-dual2(M+1);
                        if (slopef>0)
                            bound=-constant_coeff/slopef;
                            upper_pi(t-1)=min(upper_pi(t-1),bound);
                        elseif (slopef<0)
                            bound=-constant_coeff/slopef;
                            lower_pi(t-1)=max(lower_pi(t-1),bound);
                        else
                            disp('Null slope');
                            pause
                        end
                        t=t-1;
                    else
                       clear prob;
                       aux=[];
                       for j=1:M
                           aux=[aux;Probabilities(t-1,j)*Demand_Scenarios(t-1,j);-Probabilities(t-1,j)*Demand_Scenarios(t-1,j);Probabilities(t-1,j)*Demand_Scenarios(t-1,j);0];
                       end
                       for j=1:M
                           aux=[aux;Probabilities(t-1,j)];
                       end
                       prob.c=aux;
                       Lower=[];
                       Upper=[];
                       for i=1:M
                           Lower=[Lower;lower_pi(t);-Cost_b(t);-Cost_h(t);-inf];
                           Upper=[Upper;upper_pi(t);zeros(3,1)];
                       end
                       Lower=[Lower;zeros(M,1)];
                       Upper=[Upper;inf*ones(M,1)];
                       prob.blx=Lower;
                       prob.bux=Upper;
                       dynamic_buc{1,t}(M+1)=-Cost_C(t)+pis(t-1);
                       prob.blc=[dynamic_buc{1,t}(1:M+1);-large_bound*ones(M*iter,1)];
                       prob.buc=dynamic_buc{1,t};
                       prob.a = sparse(dynamic_subi_a{1,t},dynamic_subj_a{1,t},dynamic_valij_a{1,t},M+1+M*iter,5*M);
                       [~,res]=mosekopt('maximize echo(0)',prob);
                       sol=res.sol.bas.xx;
                       Alea_Uniform=rand;
                       [~,Index] = histc(Alea_Uniform,Cum_Probas{1,t-1});
                       if (Alea_Uniform==1)
                           Index=M;
                       end
                       pis(t)=sol((Index-1)*4+1);
                       t=t+1;
                    end     
             end
       end
       
       for t=T:-1:2
        if (t==T)
            clear prob;
            aux=[];
            for i=1:M
                aux=[aux;-Probabilities(T-1,i)*Demand_Scenarios(T-1,i);Probabilities(T-1,i)*Demand_Scenarios(T-1,i);0];
            end
            prob.c=aux;
            Lower=[];
            for i=1:M
                Lower=[Lower;-Cost_b(T);-Cost_h(T);-inf];
            end
            prob.blx=Lower;
            prob.bux=zeros(3*M,1);
            prob.blc=[dynamic_buc{1,T};-Cost_C(T)+pis(T-1)];
            prob.buc=prob.blc;
            prob.a = sparse(dynamic_subi_a{1,T},dynamic_subj_a{1,T},dynamic_valij_a{1,T},M+1,3*M);
            [~,res]=mosekopt('maximize echo(0)',prob);
            sol=res.sol.bas.xx;
            dual1=res.sol.bas.slc;
            dual2=res.sol.bas.suc;
            dual3=res.sol.bas.slx;
            dual4=res.sol.bas.sux;
            interceptD=Cost_C(T)*sum(dual1(1:M)-dual2(1:M))-Cost_C(T)*(dual1(M+1)-dual2(M+1));
            for k=1:M
                interceptD=interceptD-Cost_b(T)*dual3(1+3*(k-1))-Cost_h(T)*dual3(2+3*(k-1));
            end
            slope=dual1(M+1)-dual2(M+1);  
        else
            clear prob;
            aux=[];
            for j=1:M
                aux=[aux;Probabilities(t-1,j)*Demand_Scenarios(t-1,j);-Probabilities(t-1,j)*Demand_Scenarios(t-1,j);Probabilities(t-1,j)*Demand_Scenarios(t-1,j);0];
            end
            for j=1:M
                aux=[aux;Probabilities(t-1,j)];
            end
            prob.c=aux;
            Lower=[];
            Upper=[];
            for i=1:M
                Lower=[Lower;lower_pi(t);-Cost_b(t);-Cost_h(t);-inf];
                Upper=[Upper;upper_pi(t);zeros(3,1)];
            end
            Lower=[Lower;-inf*ones(M,1)];
            Upper=[Upper;inf*ones(M,1)];
            prob.blx=Lower;
            prob.bux=Upper;
            dynamic_buc{1,t}(M+1)=-Cost_C(t)+pis(t-1);
            prob.blc=[dynamic_buc{1,t}(1:M+1);-inf*ones(M*(iter+1),1)];
            prob.buc=dynamic_buc{1,t};
            prob.a = sparse(dynamic_subi_a{1,t},dynamic_subj_a{1,t},dynamic_valij_a{1,t},M+1+M*(iter+1),5*M);
            [~,res]=mosekopt('maximize echo(0)',prob);
            sol=res.sol.bas.xx;
            dual1=res.sol.bas.slc;
            dual2=res.sol.bas.suc;
            dual3=res.sol.bas.slx;
            dual4=res.sol.bas.sux;
            interceptD=Cost_C(t)*sum(dual1(1:M)-dual2(1:M))-Cost_C(t)*(dual1(M+1)-dual2(M+1))-(prob.buc(M+2:M+1+M*(iter+1)))'*(dual2(M+2:M+1+M*(iter+1)));
            for k=1:M
                interceptD=interceptD-Cost_b(t)*dual3(2+4*(k-1))-Cost_h(t)*dual3(3+4*(k-1))-Upper(4*(k-1)+1)*dual4(4*(k-1)+1)+Lower(4*(k-1)+1)*dual3(4*(k-1)+1);
            end
            slope=dual1(M+1)-dual2(M+1);
        end
        
        if (t==2)
            dynamic_subi_a{1,1}=[dynamic_subi_a{1,1},iter+2,iter+2];
            dynamic_subj_a{1,1}=[dynamic_subj_a{1,1},1,5];
            dynamic_valij_a{1,1}=[dynamic_valij_a{1,1},-slope,1];
            dynamic_buc{1,1}=[dynamic_buc{1,1};interceptD];   
        else
            %t=3,...,T
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
    prob.blx=[lower_pi(1);-Cost_b(1);-Cost_h(1);-large_bound;-inf];
    prob.bux=[upper_pi(1);0;0;0;inf];
    prob.blc=[Cost_C(1);-inf*ones(iter+1,1)];
    prob.a = sparse(dynamic_subi_a{1,1},dynamic_subj_a{1,1},dynamic_valij_a{1,1},2+iter,5);
    %disp('buc')
    prob.buc=dynamic_buc{1,1};
    %disp('size prob.a')
    %size(prob.a)
    %[full(prob.a) prob.blc prob.buc]
    %pause
    [~,res]=mosekopt('maximize echo(0)',prob);
    sol=res.sol.bas.xx;
    ub=sol'*prob.c;
    upper_bounds=[upper_bounds;ub];
    time=[time;toc];
    iter=iter+1;
end

upper_bounds=upper_bounds-Cost_C(1)*Init_Stock;