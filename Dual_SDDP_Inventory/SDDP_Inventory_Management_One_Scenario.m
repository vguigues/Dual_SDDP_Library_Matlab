
         
function [Time,Optimal_Value,Zsups,Zinfs]=SDDP_Inventory_Management_One_Scenario(Cost_C,Init_Stock,Cost_b,Cost_h,Probabilities,tol,Demand_First_Stage,Demand_Scenarios,T,M,talpha)

Iter=1;

Zsups=[];
Zinfs=[];

Time=[]; 

zsup=inf;
zinf=-inf;

Cut_Slopes=cell(T-1,1);
Cut_Intercepts=cell(T-1,1);

DynamicSubi=cell(1,T-1);
DynamicSubj=cell(1,T-1);
DynamicValij=cell(1,T-1);

blc_Dynamic=cell(1,T-1);

End_Algo=1;

Nb_Trial_Points=0;

Trial_States=cell(1,T-1);

Cum_Probas=zeros(T-1,M+1);
for t=1:T-1
    Cum_Probas(t,:)=[0,cumsum(Probabilities(t,:))];
end

Costs=[];

%while End_Algo
while End_Algo
Iter
    CurrentCost=0;
    Nb_Trial_Points=Nb_Trial_Points+1;
    tic    
    for t=1:T
        if (t==1)
            Sampled_Demand=Demand_First_Stage;
        else
            Alea_Uniform=rand;
            [~,Index] = histc(Alea_Uniform,Cum_Probas(t-1,:));
            if (Alea_Uniform==1)
                Sampled_Demand=Demand_Scenarios(t-1,M);
            else
                Sampled_Demand=Demand_Scenarios(t-1,Index);
            end
        end
        
        if ((Iter>1)&(t<T))
            
            c=[Cost_C(t);Cost_b(t);Cost_h(t);0;1];
            
            subi=[1,1];
            subj=[1,3];
            valij=[-1,1];
            blc=[-Sampled_Demand];
            buc=[inf];
            
            subi=[subi,2,2];
            subj=[subj,1,2];
            valij=[valij,1,1];
            blc=[blc;Sampled_Demand];
            buc=[buc;inf];
            
            subi=[subi,3,3];
            subj=[subj,1,4];
            valij=[valij,1,-1];
            blc=[blc;Sampled_Demand];
            buc=[buc;Sampled_Demand];
            
            subi=[subi,DynamicSubi{1,t}];
            subj=[subj,DynamicSubj{1,t}];
            valij=[valij,DynamicValij{1,t}];
            
            
            blc=[blc;blc_Dynamic{1,t}];
            buc=[buc;inf*ones(Iter-1,1)];
            
            
            if (t==1)
                blx=[Init_Stock;0;0;-inf;0];
            else
                blx=[Trial_States{1,t-1}(Nb_Trial_Points);0;0;-inf;0];
            end
            bux=inf*ones(5,1);
            
        else
            if (t==T)
                c=[Cost_C(t);Cost_b(t);Cost_h(t)];
                
                subi=[1,1];
                subj=[1,3];
                valij=[-1,1];
                blc=[-Sampled_Demand];
                buc=[inf];
                
                subi=[subi,2,2];
                subj=[subj,1,2];
                valij=[valij,1,1];
                blc=[blc;Sampled_Demand];
                buc=[buc;inf];
                
                blx=[Trial_States{1,t-1}(Nb_Trial_Points);0;0];
                bux=inf*ones(3,1);
                
            else
                c=[Cost_C(t);Cost_b(t);Cost_h(t);0;1];
                
                subi=[1,1];
                subj=[1,3];
                valij=[-1,1];
                blc=[-Sampled_Demand];
                buc=[inf];
                
                subi=[subi,2,2];
                subj=[subj,1,2];
                valij=[valij,1,1];
                blc=[blc;Sampled_Demand];
                buc=[buc;inf];
                
                subi=[subi,3,3];
                subj=[subj,1,4];
                valij=[valij,1,-1];
                blc=[blc;Sampled_Demand];
                buc=[buc;Sampled_Demand];
                
                if (t==1)
                    blx=[Init_Stock;0;0;-inf;0];
                else
                    blx=[Trial_States{1,t-1}(Nb_Trial_Points);0;0;-inf;0];
                end
                bux=inf*ones(5,1);
            end
        end
        
        clear prob;
        prob.c=c;
        prob.blc=blc;
        prob.buc=buc;
        prob.blx=blx;
        prob.bux=bux;
        if ((Iter>1)&(t<T))
            prob.a = sparse(subi,subj,valij,2+Iter,5);
        elseif (t==T)
            prob.a = sparse(subi,subj,valij,2,3);
        else
            prob.a = sparse(subi,subj,valij,3,5);
        end
        
        [r,res]=mosekopt('minimize echo(0)',prob);
        if (res.sol.bas.solsta==5)
            disp('Infeasible')
        end
        Solution=res.sol.bas.xx;
        CurrentCost=CurrentCost+Cost_C(t)*Solution(1)+Cost_b(t)*Solution(2)+Cost_h(t)*Solution(3);
        if (t==1)
            CurrentCost=CurrentCost-Cost_C(1)*Init_Stock;
        else
            CurrentCost=CurrentCost-Cost_C(t)*Trial_States{1,t-1}(Nb_Trial_Points);
        end
        if (t<T)
            Trial_States{1,t}(Nb_Trial_Points)=Solution(4);
        end
    end
    Costs=[Costs;CurrentCost];
    
    Mean_Cost=mean(Costs);
    Sigma_Cost=sqrt(var(Costs));
    zsup=Mean_Cost+talpha*Sigma_Cost/sqrt(Iter);
    Zsups=[Zsups,zsup];
    
    for t=T:-1:2
        %t
        Coupe=Cost_C(t);
        Intercept=0;
        for j=1:M
            Sampled_Demand=Demand_Scenarios(t-1,j);
            
            subi=[1,1];
            subj=[1,3];
            valij=[-1,1];
            blc=[-Sampled_Demand];
            buc=[inf];
            
            subi=[subi,2,2];
            subj=[subj,1,2];
            valij=[valij,1,1];
            blc=[blc;Sampled_Demand];
            buc=[buc;inf];
            
            
            if (t<T)
                
                subi=[subi,3,3];
                subj=[subj,1,4];
                valij=[valij,1,-1];
                blc=[blc;Sampled_Demand];
                buc=[buc;Sampled_Demand];
                
                subi=[subi,DynamicSubi{1,t}];
                subj=[subj,DynamicSubj{1,t}];
                valij=[valij,DynamicValij{1,t}];
                
                blc=[blc;blc_Dynamic{1,t}];
                buc=[buc;inf*ones(Iter,1)];
                
            end
            if (t==T)
                c=[Cost_C(t);Cost_b(t);Cost_h(t)];
                blx=[Trial_States{1,T-1}(Iter);0;0];
                bux=inf*ones(3,1);
            else
                c=[Cost_C(t);Cost_b(t);Cost_h(t);0;1];
                blx=[Trial_States{1,t-1}(Iter);0;0;-inf;1];
                bux=inf*ones(5,1);
            end
            
            clear prob;
            prob.c=c;
            prob.blc=blc;
            prob.buc=buc;
            prob.blx=blx;
            prob.bux=bux;
            
            if (t==T)
                prob.a = sparse(subi,subj,valij,2,3);
            else
                prob.a = sparse(subi,subj,valij,3+Iter,5);
            end
            
            [r,res]=mosekopt('minimize echo(0)',prob);
            
            Solution=res.sol.bas.xx;
            dual1=res.sol.bas.slc;
            dual2=res.sol.bas.suc;
            dual3=res.sol.bas.slx;
            Coupe=Coupe-Probabilities(t-1,j)*dual3(1);
            
            Intercept=Intercept+Probabilities(t-1,j)*(-Sampled_Demand*dual1(1)+Sampled_Demand*dual1(2));
            
            if (t<T)
                Intercept=Intercept+Probabilities(t-1,j)*Sampled_Demand*(dual1(3)-dual2(3));
                for i=1:Iter
                    Intercept=Intercept+blc_Dynamic{1,t}(i)*Probabilities(t-1,j)*dual1(3+i);
                end
            end
        end
        Cut_Slopes{t-1,1}=[Cut_Slopes{t-1,1};-Coupe];
        Cut_Intercepts{t-1,1}=[Cut_Intercepts{t-1,1};Intercept];
        DynamicSubi{1,t-1}=[DynamicSubi{1,t-1},(3+Iter)*ones(1,2)];
        DynamicSubj{1,t-1}=[DynamicSubj{1,t-1},4,5];
        DynamicValij{1,t-1}=[DynamicValij{1,t-1},Coupe,1];
        blc_Dynamic{1,t-1}=[blc_Dynamic{1,t-1};Intercept];
    end
    
    Sampled_Demand=Demand_First_Stage;
    
    subi=[1,1];
    subj=[1,3];
    valij=[-1,1];
    blc=[-Sampled_Demand];
    buc=[inf];
    
    subi=[subi,2,2];
    subj=[subj,1,2];
    valij=[valij,1,1];
    blc=[blc;Sampled_Demand];
    buc=[buc;inf];
    
    subi=[subi,3,3];
    subj=[subj,1,4];
    valij=[valij,1,-1];
    blc=[blc;Sampled_Demand];
    buc=[buc;Sampled_Demand];
    
    
    subi=[subi,DynamicSubi{1,1}];
    subj=[subj,DynamicSubj{1,1}];
    valij=[valij,DynamicValij{1,1}];
    
    blc=[blc;blc_Dynamic{1,1}];
    buc=[buc;inf*ones(Iter,1)];
    
    c=[Cost_C(1);Cost_b(1);Cost_h(1);0;1];
    
    blx=[Init_Stock;0;0;-inf;0];
    bux=inf*ones(5,1);
    
    clear prob;
    prob.c=c;
    prob.blc=blc;
    
    prob.buc=buc;
    prob.blx=blx;
    prob.bux=bux;
    prob.a = sparse(subi,subj,valij,3+Iter,5);
    
    [r,res]=mosekopt('minimize echo(0)',prob);
    
    Solution=res.sol.bas.xx;
    zinf=c'*Solution-Cost_C(1)*Init_Stock;
    Zinfs=[Zinfs;zinf];
    End_Algo=abs((zsup-zinf)/zsup)>tol;
    if (Iter>1000)
        End_Algo=0;
    end    
    Iter=Iter+1;
    Time=[Time;toc];
end

Optimal_Value=zinf;



