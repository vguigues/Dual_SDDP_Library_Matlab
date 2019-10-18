

function [time,optimal_value,Bellman_Functions,List_Abs_Bell]=Dynamic_Programming_Inventory(Cost_C,Init_Stock,Cost_b,Cost_h,Probabilities,Demand_First_Stage,Demand_Scenarios,T,M,large_bound,xmin,xmax,Nbpoints)

Bellman_Functions=cell(1,T-1);
List_Abs_Bell=cell(1,T-1);

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

for i=1:M
    dynamic_subi{1,T}=[dynamic_subi{1,T},M+1];
    dynamic_subj{1,T}=[dynamic_subj{1,T},3*i];
    dynamic_valij{1,T}=[dynamic_valij{1,T},Probabilities(T-1,i)];
end

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
        dynamic_valij{1,t}=[dynamic_valij{1,t},Probabilities(t-1,  i)];
        dynamic_c{1,t}=[dynamic_c{1,t};Probabilities(t-1,i)];
    end
    dynamic_buc{1,t}=[Cost_C(t)*ones(M,1);0];
end

dynamic_subi{1,1}=[1,1,1,1];
dynamic_subj{1,1}=[1,2,3,4];
dynamic_valij{1,1}=[1,-1,1,-1];
dynamic_buc{1,1}=[Cost_C(1)];
dynamic_c{1,1}=[Demand_First_Stage;-Demand_First_Stage;Demand_First_Stage;-Init_Stock;1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dynamic_feas_subi=cell(1,T-1);
dynamic_feas_subj=cell(1,T-1);
dynamic_feas_valij=cell(1,T-1);

upper_pi=large_bound*ones(T-1,1);
lower_pi=-large_bound*ones(T-1,1);

dynamic_feas_subi{1,T-1}=[];
dynamic_feas_subj{1,T-1}=[];
dynamic_feas_valij{1,T-1}=[];

for i=1:M
    dynamic_feas_subi{1,T-1}=[dynamic_feas_subi{1,T-1},i,i,i];
    dynamic_feas_subj{1,T-1}=[dynamic_feas_subj{1,T-1},3*(i-1)+1,3*(i-1)+2,3*i];
    dynamic_feas_valij{1,T-1}=[dynamic_feas_valij{1,T-1},-1,1,-1];
end

for i=1:M
    dynamic_feas_subi{1,T-1}=[dynamic_feas_subi{1,T-1},M+1];
    dynamic_feas_subj{1,T-1}=[dynamic_feas_subj{1,T-1},3*i];
    dynamic_feas_valij{1,T-1}=[dynamic_feas_valij{1,T-1},Probabilities(T-1,i)];
end
dynamic_feas_subi{1,T-1}=[dynamic_feas_subi{1,T-1},M+1,M+1];
dynamic_feas_subj{1,T-1}=[dynamic_feas_subj{1,T-1},3*M+1,3*M+2];
dynamic_feas_valij{1,T-1}=[dynamic_feas_valij{1,T-1},1,-1];


for t=2:T-1

dynamic_feas_subi{1,t-1}=[];
dynamic_feas_subj{1,t-1}=[];
dynamic_feas_valij{1,t-1}=[];

for i=1:M
    dynamic_feas_subi{1,t-1}=[dynamic_feas_subi{1,t-1},i,i,i,i];
    dynamic_feas_subj{1,t-1}=[dynamic_feas_subj{1,t-1},4*(i-1)+1,4*(i-1)+2,4*(i-1)+3,4*i];
    dynamic_feas_valij{1,t-1}=[dynamic_feas_valij{1,t-1},1,-1,1,-1];
end

for i=1:M
    dynamic_feas_subi{1,t-1}=[dynamic_feas_subi{1,t-1},M+1];
    dynamic_feas_subj{1,t-1}=[dynamic_feas_subj{1,t-1},4*i];
    dynamic_feas_valij{1,t-1}=[dynamic_feas_valij{1,t-1},Probabilities(t-1,i)];
end

dynamic_feas_subi{1,t-1}=[dynamic_feas_subi{1,t-1},M+1,M+1];
dynamic_feas_subj{1,t-1}=[dynamic_feas_subj{1,t-1},4*M+1,4*M+2];
dynamic_feas_valij{1,t-1}=[dynamic_feas_valij{1,t-1},1,-1];

end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pas=(xmax(T-1)-xmin(T-1))/Nbpoints;
abscisses=[xmin(T-1):pas:xmax(T-1)];
found_first_feasible=0;
index=1;
tic
while ((found_first_feasible==0)&&(index<=Nbpoints+1))
   %index
   %Solve the problem for T at abscisses(index)
   clear prob;
   prob.c=[zeros(3*M,1);1;1];
   dynamic_buc{1,T}(M+1)=-Cost_C(T)+abscisses(index);
   prob.blc=dynamic_buc{1,T};
   prob.buc=prob.blc;
   prob.blx=[repmat([-Cost_b(T);-Cost_h(T);-inf],M,1);0;0];
   prob.bux=[zeros(3*M,1);inf;inf];
   prob.a=sparse(dynamic_feas_subi{1,T-1},dynamic_feas_subj{1,T-1},dynamic_feas_valij{1,T-1},M+1,3*M+2);
   [r,res]=mosekopt('minimize echo(0)',prob);
   sol=res.sol.bas.xx;
   error=sol(3*M+1)+sol(3*M+2);
   if (error==0)
       found_first_feasible=1;
       List_Abs_Bell{1,T-1}=[abscisses(index)];
       lower_pi(T-1)=abscisses(index);
       clear prob;
       prob.c=dynamic_c{1,T};
       dynamic_buc{1,T}(M+1)=-Cost_C(T)+abscisses(index);
       prob.blc=dynamic_buc{1,T};
       prob.buc=prob.blc;
       prob.blx=repmat([-Cost_b(T);-Cost_h(T);-inf],M,1);
       prob.bux=zeros(3*M,1);
       prob.a = sparse(dynamic_subi{1,T},dynamic_subj{1,T},dynamic_valij{1,T},M+1,3*M);
       [r,res]=mosekopt('maximize echo(0)',prob);
       sol=res.sol.bas.xx;
       prev_V=sol'*prob.c;
       Bellman_Functions{1,T-1}=[prev_V];
   end 
   index=index+1;
end        


found_first_unfeasible=0;
nb_cuts=zeros(T-1,1);

while ((found_first_unfeasible==0)&&(index<=Nbpoints+1))
   %index
    %Solve the problem for T at abscisses(index)
   clear prob;
   prob.c=[zeros(3*M,1);1;1];
   dynamic_buc{1,T}(M+1)=-Cost_C(T)+abscisses(index);
   prob.blc=dynamic_buc{1,T};
   prob.buc=prob.blc;
   prob.blx=[repmat([-Cost_b(T);-Cost_h(T);-inf],M,1);0;0];
   prob.bux=[zeros(3*M,1);inf;inf];
   prob.a = sparse(dynamic_feas_subi{1,T-1},dynamic_feas_subj{1,T-1},dynamic_feas_valij{1,T-1},M+1,3*M+2);
   [r,res]=mosekopt('minimize echo(0)',prob);
   sol=res.sol.bas.xx;
   error=sol(3*M+1)+sol(3*M+2);
   if (error>(10^(-10)))
       upper_pi(T-1)=abscisses(index-1);
       found_first_unfeasible=1;
   else
       List_Abs_Bell{1,T-1}=[List_Abs_Bell{1,T-1};abscisses(index)];
       clear prob;
       prob.c=dynamic_c{1,T};
       dynamic_buc{1,T}(M+1)=-Cost_C(T)+abscisses(index);
       prob.blc=dynamic_buc{1,T};
       prob.buc=prob.blc;
       prob.blx=repmat([-Cost_b(T);-Cost_h(T);-inf],M,1);
       prob.bux=zeros(3*M,1);
       prob.a = sparse(dynamic_subi{1,T},dynamic_subj{1,T},dynamic_valij{1,T},M+1,3*M);
       [r,res]=mosekopt('maximize echo(0)',prob);
       sol=res.sol.bas.xx;
       current_V=sol'*prob.c;
       beta=(current_V-prev_V)/(abscisses(index)-abscisses(index-1));
       theta=prev_V-beta*abscisses(index-1);
       Bellman_Functions{1,T-1}=[Bellman_Functions{1,T-1};current_V];
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

for t=T-1:-1:2
    t
    pas=(xmax(t-1)-xmin(t-1))/Nbpoints;
    abscisses=[xmin(t-1):pas:xmax(t-1)];
    found_first_feasible=0;
    index=1;
    while ((found_first_feasible==0)&&(index<=Nbpoints+1))
        %Solve the problem for t at abscisses(index)        
        clear prob;
        prob.c=[zeros(4*M,1);1;1];
        prob.blc=[Cost_C(t)*ones(M,1);-Cost_C(t)+abscisses(index)];
        prob.buc=prob.blc;
        prob.blx=[repmat([lower_pi(t);-Cost_b(t);-Cost_h(t);-inf],M,1);0;0];
        prob.bux=[repmat([upper_pi(t);zeros(3,1)],M,1);inf;inf];
        prob.a=sparse(dynamic_feas_subi{1,t-1},dynamic_feas_subj{1,t-1},dynamic_feas_valij{1,t-1},M+1,4*M+2);
        [r,res]=mosekopt('minimize echo(0)',prob);
        sol=res.sol.bas.xx;
        error=abs(sol(4*M+1))+abs(sol(4*M+2));
        if (error==0)
            found_first_feasible=1;
            List_Abs_Bell{1,t-1}=[abscisses(index)];
            lower_pi(t-1)=abscisses(index);
            clear prob;
            prob.c=dynamic_c{1,t};
            dynamic_buc{1,t}(M+1)=-Cost_C(t)+abscisses(index);
            prob.blc=[dynamic_buc{1,t}(1:M+1);-inf*ones(nb_cuts(t)*M,1)];
            prob.buc=dynamic_buc{1,t};
            prob.blx=[repmat([lower_pi(t);-Cost_b(t);-Cost_h(t);-inf],M,1);-inf*ones(M,1)];
            prob.bux=[repmat([upper_pi(t);zeros(3,1)],M,1);inf*ones(M,1)];
            prob.a = sparse(dynamic_subi{1,t},dynamic_subj{1,t},dynamic_valij{1,t},M+1+nb_cuts(t)*M,5*M);
            [r,res]=mosekopt('maximize echo(0)',prob);
            sol=res.sol.bas.xx;
            prev_V=sol'*prob.c;
            Bellman_Functions{1,t-1}=[prev_V];
        end
        index=index+1;
    end
    
    found_first_unfeasible=0;
     
    while ((found_first_unfeasible==0)&&(index<=Nbpoints+1))
        %Solve the problem for t at abscisses(index)
        clear prob;
        prob.c=[zeros(4*M,1);1;1];
        prob.blc=[Cost_C(t)*ones(M,1);-Cost_C(t)+abscisses(index)];
        prob.buc=prob.blc;
        prob.blx=[repmat([lower_pi(t);-Cost_b(t);-Cost_h(t);-inf],M,1);0;0];
        prob.bux=[repmat([upper_pi(t);zeros(3,1)],M,1);inf;inf];
        prob.a=sparse(dynamic_feas_subi{1,t-1},dynamic_feas_subj{1,t-1},dynamic_feas_valij{1,t-1},M+1,4*M+2);
        [r,res]=mosekopt('minimize echo(0)',prob);
        sol=res.sol.bas.xx;
        error=abs(sol(4*M+1))+abs(sol(4*M+2));
        if (error>0)
            upper_pi(t-1)=abscisses(index-1);
            found_first_unfeasible=1;
        else
            List_Abs_Bell{1,t-1}=[List_Abs_Bell{1,t-1};abscisses(index)];
            clear prob;
            prob.c=dynamic_c{1,t};
            dynamic_buc{1,t}(M+1)=-Cost_C(t)+abscisses(index);
            prob.blc=[dynamic_buc{1,t}(1:M+1);-inf*ones(nb_cuts(t)*M,1)];
            prob.buc=dynamic_buc{1,t};
            prob.blx=[repmat([lower_pi(t);-Cost_b(t);-Cost_h(t);-inf],M,1);-inf*ones(M,1)];
            prob.bux=[repmat([upper_pi(t);zeros(3,1)],M,1);inf*ones(M,1)];
            prob.a = sparse(dynamic_subi{1,t},dynamic_subj{1,t},dynamic_valij{1,t},M+1+nb_cuts(t)*M,5*M);
            [r,res]=mosekopt('maximize echo(0)',prob);
            sol=res.sol.bas.xx;
            current_V=sol'*prob.c;
            beta=(current_V-prev_V)/(abscisses(index)-abscisses(index-1));
            theta=prev_V-beta*abscisses(index-1);
            Bellman_Functions{1,t-1}=[Bellman_Functions{1,t-1};current_V];
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

clear prob;
prob.c=dynamic_c{1,1};
prob.buc=dynamic_buc{1,1};
prob.blc=[Cost_C(1);-inf*ones(nb_cuts(1),1)];
%prob.blx=[lower_pi(1);-Cost_b(1);-Cost_h(1);-large_bound;-inf];
prob.blx=[lower_pi(1);-Cost_b(1);-Cost_h(1);-inf;-inf];
prob.bux=[upper_pi(1);zeros(3,1);inf];
prob.a = sparse(dynamic_subi{1,1},dynamic_subj{1,1},dynamic_valij{1,1},1+nb_cuts(1),5);
[r,res]=mosekopt('maximize echo(0)',prob);
sol=res.sol.bas.xx;
optimal_value=sol'*prob.c-Cost_C(1)*Init_Stock;

time=toc;




  