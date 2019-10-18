
%addpath 'C:\Users\vguigues6\Dropbox\Articles_Math\Vincent_Alex\Dual_SDDP\Dual_SDDP_Library\Inventory\Dual_SDDP_Inventory'
addpath 'C:\Users\vince\Dropbox\Articles_Math\Vincent_Alex\Dual_SDDP\Dual_SDDP_Library\Inventory\Dual_SDDP_Inventory'
addpath 'C:\Program Files\Mosek\9.0'
addpath 'C:\Program Files\Mosek\9.0\toolbox\r2015a'

rng('shuffle')
T=20;
Vec=[1:1:T]';
Cost_C=1.5+cos((pi/6)*Vec);
Demands=5+0.5*Vec;
Init_Stock=10;
Cost_b=2.8*ones(T,1);
Cost_h=0.2*ones(T,1);
M=10;
Probabilities=zeros(T-1,M);
for i=1:(T-1)
    Probabilities(i,:)=(1/M)*ones(1,M);
end
tol=0.05;
Demand_First_Stage=Demands(1);
Demand_Scenarios=zeros(T-1,M);
for t=1:(T-1)
    for j=1:M
        Demand_Scenarios(t,j)=Demands(t+1)*(1.5+0.1*randn);
    end
end
talpha=1.96;
nb_iter_max=600;
penaltiesQ=ones(T-1,nb_iter_max);
penaltiesR=ones(T-1,nb_iter_max);
for t=1:T-1
    for i=1:nb_iter_max
        penaltiesQ(t,i)=1000;
        penaltiesR(t,i)=1000;
    end
end
large_bound=500;
big_m=20000;
xmin=-5*ones(T-1,1);
xmax=5*ones(T-1,1);
Nbpoints=2000;



[optimal_value,Bellman_Functions,List_Abs_Bell]=Dynamic_Programming_Inventory(Cost_C,Init_Stock,Cost_b,Cost_h,Probabilities,Demand_First_Stage,Demand_Scenarios,T,M,large_bound,xmin,xmax,Nbpoints);
subplot(3,2,1);

for t=1:T-1
    Bellman_Functionsp{t,2}=Bellman_Functionsp100{t,1};
end

Bellman_Functionsp100



plot([1:10],Zinfs(1:10),'-k')
hold on
plot([1:10],upper_bounds(1:10),'r--');
hold on
plot([1:10],upper_boundsf(1:10),'m:');
legend(['Primal SDDP                    ';'Dual SDDP with penalties       ';'Dual SDDP with feasibility cuts';])
hold on
plot([1:10],Zsups(1:10),'-k')

plot([1:600],cumsum(Time),'-k');
hold on
plot([1:600],cumsum(time),'r--');
hold on
plot([1:600],cumsum(timef),'m:');
hold on

Z1=[Zinfs(2);Zinfs(3);Zinfs(5);Zinfs(10);Zinfs(50);Zinfs(100);Zinfs(200);Zinfs(300);Zinfs(400);Zinfs(500);Zinfs(600)];
Z2=[Zsups(2);Zsups(3);Zsups(5);Zsups(10);Zsups(50);Zsups(100);Zsups(200);Zsups(300);Zsups(400);Zsups(500);Zsups(600)];
Z3=[upper_boundsf(2);upper_boundsf(3);upper_boundsf(5);upper_boundsf(10);upper_boundsf(50);upper_boundsf(100);upper_boundsf(200);upper_boundsf(300);upper_boundsf(400);upper_boundsf(500);upper_boundsf(600)]
Z4=[upper_bounds(2);upper_bounds(3);upper_bounds(5);upper_bounds(10);upper_bounds(50);upper_bounds(100);upper_bounds(200);upper_bounds(300);upper_bounds(400);upper_bounds(500);upper_bounds(600)]
[Z1 Z2 Z3 Z4]

plot([11:100],Zinfs(11:100),'-k')
hold on
plot([11:100],upper_bounds(11:100),'r--');
hold on
plot([11:100],upper_boundsf(11:100),'m:');
legend(['Primal SDDP                    ';'Dual SDDP with penalties       ';'Dual SDDP with feasibility cuts';])
hold on
plot([11:100],Zsups(11:100),'-k')

upper_boundsp
timep

upper_bounds 
time

Zinfs
Zsups
Time

err=(pupper_bounds(1)-plower_bounds(1))/pupper_bounds(1);
index=2;
while (err>0.07)
       err=(pupper_bounds(index)-plower_bounds(index))/pupper_bounds(index);
       index=index+1;
end
TimeSDDP=sum(ptime(1:index-1));

cumsumU=cumsum(time109);
%cumsumU=cumsum(timef);
cumsumL=cumsum(ptime);

pu=10^(10);
tc=0;
pl=-inf;
iu=1;
il=1;
err=inf;
while (err>0.07)
   if (cumsumU(iu)<cumsumL(il))
       tc=cumsumU(iu)
       pu=upper_bounds109(iu);
       err=(pu-pl)/pu;
       iu=iu+1;
   elseif (cumsumL(il)<cumsumU(iu))
       tc=cumsumL(il)
       pl=plower_bounds(il);
       err=(pu-pl)/pu;
       il=il+1;
   else
       tc=cumsumL(il)
       pu=upper_bounds109(iu);
       pl=plower_bounds(il);
       err=(pu-pl)/pu;
       il=il+1;
       iu=iu+1;
   end 
end

Lmins=zeros(1,T-1);
Lmaxs=zeros(1,T-1);
Lmean=zeros(1,T-1);
for t=1:T-1
    Lmins(t)=min(Errors{1,t});
    Lmaxs(t)=max(Errors{1,t});
    Lmean(t)=mean(Errors{1,t});
end

plot([2:T],Lmins);
hold on
plot([2:T],Lmaxs)
hold on
plot([2:T],Lmean)

for i=1:200
    traj=[];
    for t=1:T-1
        traj=[traj;Errors{1,t}(i)];
    end
    plot([2:T],traj,'-k');
    hold on
end


t=2;
hold on
plot(List_Abs_Bell{1,t},Bellman_Functions{1,t},'-k','Linewidth',4);
hold on
pas=10/Nbpoints;
abscisses=[-5:pas:5];
plot(abscisses,Bellman_Functionsp{t,1},'--r');
hold on
plot(abscisses,Bellman_Functionsp{t,2},'-.m');
hold on
plot(abscisses,Bellman_Functionsp{t,3},':b');
hold on

plot([10:100],Zinfs(10:100),'k-');
hold on
plot([10:100],upper_boundsf(10:100),'m:');
hold on
plot([10:100],upper_bounds(10:100),'r--');
hold on
plot([10:100],Zsups(10:100),'k-');



subplot(3,2,2);
plot(List_Abs_Bell{1,2},Bellman_Functions{1,2});
subplot(3,2,3);
plot(List_Abs_Bell{1,3},Bellman_Functions{1,3});
subplot(3,2,4);
plot(List_Abs_Bell{1,4},Bellman_Functions{1,4});
subplot(3,2,5);
plot(List_Abs_Bell{1,5},Bellman_Functions{1,5});
subplot(3,2,6);
plot(List_Abs_Bell{1,6},Bellman_Functions{1,6});



