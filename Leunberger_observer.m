%Second Part--Observability
close all
sim_time=10;

%Declaration of C matrices/vector for corresponding output vectors
Cx=[0 0 0 1 0 0]  %----observable;
Ct1t2=[0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 1 0;0 0 0 0 0 1];
Cxt2=[0 0 0 1 0 0;0 0 0 0 0 1] %----observable;
Cxt1t2=[0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1] %----observable;


%Checking the condition for observability
O1=rank([Cx; Cx*A_lin; Cx*A_lin^2; Cx*A_lin^3; Cx*A_lin^4; Cx*A_lin^5]);
O2=rank([Ct1t2; Ct1t2*A_lin; Ct1t2*A_lin^2; Ct1t2*A_lin^3; Ct1t2*A_lin^4; Ct1t2*A_lin^5]);
O3=rank([Cxt2; Cxt2*A_lin; Cxt2*A_lin^2; Cxt2*A_lin^3; Cxt2*A_lin^4; Cxt2*A_lin^5]);
O4=rank([Cxt1t2; Cxt1t2*A_lin; Cxt1t2*A_lin^2; Cxt1t2*A_lin^3; Cxt1t2*A_lin^4; Cxt1t2*A_lin^5]);

%Step input(for simulation)
x_fin2=[];
syms tau
integral=vpa(int(expm(-A_lin).*tau));
for t=0:time_step:sim_time
    x_fin_temp3 = expm((A_lin).*t)*x_init + expm((A_lin).*t)*subs(integral,tau ,t)*B_lin ;
    x_fin2= [x_fin2 x_fin_temp3];
    t
end
figure
hold on
title('')
plot([0:time_step:sim_time],x_fin2(1,:),'b');
plot([0:time_step:sim_time],x_fin2(2,:),'g');
plot([0:time_step:sim_time],x_fin2(3,:),'r');

%Leunberger Observer

%For the first controllability matrix(first vector of output)
x_error=x_init;
L=[ 30;
    12.5;
    8 ;
    11;
    0.2 ;
    0.08];

% L=[ 30;
%     0.7;
%     10 ;
%     6;
%     1 ;
%     0.01];
x_error_dot=(A_lin-L*Cx)*x_error;

% Solution for first output vector
x_fin3=[];
for t=0:time_step:sim_time
    x_fin_temp = expm((A_lin-L*Cx).*t)*x_error;
    x_fin3= [x_fin3 x_fin_temp];
    
end

figure
title('For First output vector(Linear)')
hold on
plot([0:time_step:sim_time],x_fin2(1,:)-x_fin3(1,:),'--','Color',[1,0,0]);
plot([0:time_step:sim_time],x_fin2(1,:),'r');

plot([0:time_step:sim_time],x_fin2(2,:)-x_fin3(2,:),'--','Color',[0,1,0]);
plot([0:time_step:sim_time],x_fin2(2,:),'g');

plot([0:time_step:sim_time],x_fin2(3,:)-x_fin3(3,:),'--','Color',[0,0,1]);
plot([0:time_step:sim_time],x_fin2(3,:),'b');

legend('xchat','xc','t1hat','t1','t2hat','t2')

%plot([0:time_step:sim_time],x_fin2(3,:)-x_fin3(3,:),'r');

eig(A_lin-L*Cx)


% For the third controllability matrix(third vector of output)

L_3=[ 12 12; 
     0.55 0.55;
     12 12;
     10 10;
     0.75 0.75;
     0.05  0.05];

% L_3=[ 7.5 7.5;
%      0.3 0.3;
%      8  8;
%      6  6;
%      1  1;
%      0.2  0.2];
x_error_dot3=(A_lin-L_3*Cxt2)*x_error;

% Solution for third output vector
x_fin4=[];
for t=0:time_step:sim_time
    x_fin_temp = expm((A_lin-L_3*Cxt2).*t)*x_error;
    x_fin4= [x_fin4 x_fin_temp];
    
end

figure
title('For Third output vector(Linear)')
hold on
plot([0:time_step:sim_time],x_fin2(1,:)-x_fin4(1,:),'--','Color',[1,0,0]);
plot([0:time_step:sim_time],x_fin2(1,:),'r');

plot([0:time_step:sim_time],x_fin2(2,:)-x_fin4(1,:),'--','Color',[0,1,0]);
plot([0:time_step:sim_time],x_fin2(2,:),'g');

plot([0:time_step:sim_time],x_fin2(3,:)-x_fin4(1,:),'--','Color',[0,0,1]);
plot([0:time_step:sim_time],x_fin2(3,:),'b');
legend('xchat','xc','t1hat','t1','t2hat','t2')

%plot([0:time_step:sim_time],x_fin2(3,:)-x_fin3(3,:),'r');

eig(A_lin-L_3*Cxt2)



% For the fourth controllability matrix(fourth vector of output)
L_4=[ 10 10 10;
      0.5 0.5 0.5;
      10 10 10;
      6 6 6;
      1 1 1;
      0.01 0.01 0.01];
x_error_dot4=(A_lin-L_4*Cxt1t2)*x_error;

% Solution for fourth output vector
x_fin5=[];
for t=0:time_step:sim_time
    x_fin_temp = expm((A_lin-L_4*Cxt1t2).*t)*x_error;
    x_fin5= [x_fin5 x_fin_temp];
    
end

figure
title('For Fourth output vector(Linear)')
hold on
plot([0:time_step:sim_time],x_fin2(1,:)-x_fin5(1,:),'--','Color',[1,0,0]);
plot([0:time_step:sim_time],x_fin2(1,:),'r');
plot([0:time_step:sim_time],x_fin2(2,:)-x_fin5(1,:),'--','Color',[0,1,0]);
plot([0:time_step:sim_time],x_fin2(2,:),'g');
plot([0:time_step:sim_time],x_fin2(3,:)-x_fin5(1,:),'--','Color',[0,0,1]);
plot([0:time_step:sim_time],x_fin2(3,:),'b');
legend('xchat','xc','t1hat','t1','t2hat','t2')

%plot([0:time_step:sim_time],x_fin2(3,:)-x_fin3(3,:),'r');

eig(A_lin-L_4*Cxt1t2)



%L=[0 0 0  10 0 0;
%   0 0 0 0.5 0 0;
%  0 0 0 10 0 0;
%   0 0 0 5 0 0;
%   0 0 0 1 0 0;
%   0 0 0 0.01 0 0];

%L_3=[0 0 0  10 0 10;
%   0 0 0 0.5 0 0.5;
%   0 0 0 10 0 10;
%    0 0 0 5 0 5;
%    0 0 0 1 0 1;
%   0 0 0 0.01 0 0.01];

%L_4=[0 0 0  10 10 10;
%    0 0 0 0.5 0.5 0.5;
%    0 0 0 10 10 10;
%    0 0 0 5 5 5;
%    0 0 0 1 1 1;
%    0 0 0 0.01 0.01 0.01];
