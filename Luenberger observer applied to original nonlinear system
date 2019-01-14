
%Implementing the Luenberger observer for non-linear system
control_input2 = 1;
f1_sub2 = simplify(vpa(subs(f1, [F l1 l2 M m1 m2 g],[ 1 20 10 1000 100 100 9.8])),5);
f2_sub2 = simplify(vpa(subs(f2, [F l1 l2 M m1 m2 g],[ 1 20 10 1000 100 100 9.8])),5);
f3_sub2 = simplify(vpa(subs(f3, [F l1 l2 M m1 m2 g],[ 1 20 10 1000 100 100 9.8])),5);
f1_sub2_t = subs(f1_sub2,[ x1 x2 x3 x4 x5 x6],[x1t x2t x3t x4t x5t x6t]);
f2_sub2_t = subs(f2_sub2,[ x1 x2 x3 x4 x5 x6],[x1t x2t x3t x4t x5t x6t]);
 f3_sub2_t = subs(f3_sub2,[ x1 x2 x3 x4 x5 x6],[x1t x2t x3t x4t x5t x6t]);
syms x1hat(t) x2hat(t) x3hat(t) x4hat(t) x5hat(t) x6hat(t) real
f1_sub2_hat = subs(f1_sub2,[ x1 x2 x3 x4 x5 x6],[x1hat x2hat x3hat x4hat x5hat x6hat]);
f2_sub2_hat = subs(f2_sub2,[ x1 x2 x3 x4 x5 x6],[x1hat x2hat x3hat x4hat x5hat x6hat]);
f3_sub2_hat = subs(f3_sub2,[ x1 x2 x3 x4 x5 x6],[x1hat x2hat x3hat x4hat x5hat x6hat]);
LnC = L*Cx;
vec_xt= [ x1t; x2t; x3t; x4t;x5t;x6t];
vec_xhat =[ x1hat; x2hat; x3hat; x4hat; x5hat; x6hat];
diff_eqn1 = diff(x1hat,t,1)==f1_sub2_hat+ LnC(1,:)*(vec_xt-vec_xhat);
diff_eqn2 = diff(x2hat,t,1)==f2_sub2_hat+ LnC(2,:)*(vec_xt-vec_xhat);
diff_eqn3 = diff(x3hat,t,1)==f3_sub2_hat+ LnC(3,:)*(vec_xt-vec_xhat);
diff_eqn4= diff(x4hat,t,1)==x1hat + LnC(4,:)*(vec_xt-vec_xhat);
diff_eqn5 = diff(x5hat,t,1)==x2hat + LnC(5,:)*(vec_xt-vec_xhat);
diff_eqn6 = diff(x6hat,t,1)==x3hat + LnC(6,:)*(vec_xt-vec_xhat);
diff_eqn7 = diff(x1t,t,1)==  f1_sub2_t;
diff_eqn8 = diff(x2t,t,1)== f2_sub2_t;
diff_eqn9 = diff(x3t,t,1) == f3_sub2_t;
diff_eqn10= diff(x4t,t,1) == x1t;
diff_eqn11 = diff(x5t,t,1) == x2t;
diff_eqn12 = diff(x6t,t,1) == x3t;
state_vec1= [vec_xhat; vec_xt];
eqns2 = [ diff_eqn1; diff_eqn2; diff_eqn3; diff_eqn4; diff_eqn5; diff_eqn6; diff_eqn7; diff_eqn8; diff_eqn9; diff_eqn10; diff_eqn11;diff_eqn12]; 
[M3,F3] = massMatrixForm(eqns2,state_vec1);
f = M3\F3;
ode_fun = odeFunction(f,state_vec1);
x_init2 =[zeros(6,1); x_init];
[t,x_out] = ode45(ode_fun,[0 sim_time],x_init2);
figure;
title('For First output vector(non-linear)')
hold on;
plot(t,x_out(:,10),'r');
plot(t,x_out(:,4),'--','Color',[1,0,0]);
plot(t,x_out(:,11),'g');
plot(t,x_out(:,5),'--','Color',[0,1,0]);
plot(t,x_out(:,12),'b');
plot(t,x_out(:,6),'--','Color',[0,0,1]);
legend('xc','xchat','t1','t1hat','t2','t2hat')
%plot(t,x_out(:,12),'r');


% (For third output Vector)

syms x1hat(t) x2hat(t) x3hat(t) x4hat(t) x5hat(t) x6hat(t) real

LnC_3 = L_3*Cxt2;
diff_eqn1 = diff(x1hat,t,1)==f1_sub2_hat+ LnC_3(1,:)*(vec_xt-vec_xhat);
diff_eqn2 = diff(x2hat,t,1)==f2_sub2_hat+ LnC_3(2,:)*(vec_xt-vec_xhat);
diff_eqn3 = diff(x3hat,t,1)==f3_sub2_hat+ LnC_3(3,:)*(vec_xt-vec_xhat);
diff_eqn4= diff(x4hat,t,1)==x1hat + LnC_3(4,:)*(vec_xt-vec_xhat);
diff_eqn5 = diff(x5hat,t,1)==x2hat + LnC_3(5,:)*(vec_xt-vec_xhat);
diff_eqn6 = diff(x6hat,t,1)==x3hat + LnC_3(6,:)*(vec_xt-vec_xhat);
diff_eqn7 = diff(x1t,t,1)==  f1_sub2_t;
diff_eqn8 = diff(x2t,t,1)== f2_sub2_t;
diff_eqn9 = diff(x3t,t,1) == f3_sub2_t;
diff_eqn10= diff(x4t,t,1) == x1t;
diff_eqn11 = diff(x5t,t,1) == x2t;
diff_eqn12 = diff(x6t,t,1) == x3t;

eqns2 = [ diff_eqn1; diff_eqn2; diff_eqn3; diff_eqn4; diff_eqn5; diff_eqn6; diff_eqn7; diff_eqn8; diff_eqn9; diff_eqn10; diff_eqn11;diff_eqn12]; 
[M3,F3] = massMatrixForm(eqns2,state_vec1);
f = M3\F3;
ode_fun = odeFunction(f,state_vec1);
x_init3 =[zeros(6,1); x_init];
[t,x_out] = ode45(ode_fun,[0 sim_time],x_init3);
figure;
title('For Third output vector(non-linear)')
hold on;
plot(t,x_out(:,10),'r');
plot(t,x_out(:,4),'--','Color',[1,0,0]);
plot(t,x_out(:,11),'g');
plot(t,x_out(:,5),'--','Color',[0,1,0]);
plot(t,x_out(:,12),'b');
plot(t,x_out(:,6),'--','Color',[0,0,1]);
legend('xc','xchat','t1','t1hat','t2','t2hat')
%plot(t,x_out(:,12),'r');


% (For fourth output Vector)

syms x1hat(t) x2hat(t) x3hat(t) x4hat(t) x5hat(t) x6hat(t) real

LnC_4 = L_4*Cxt1t2;
diff_eqn1 = diff(x1hat,t,1)==f1_sub2_hat+ LnC_4(1,:)*(vec_xt-vec_xhat);
diff_eqn2 = diff(x2hat,t,1)==f2_sub2_hat+ LnC_4(2,:)*(vec_xt-vec_xhat);
diff_eqn3 = diff(x3hat,t,1)==f3_sub2_hat+ LnC_4(3,:)*(vec_xt-vec_xhat);
diff_eqn4= diff(x4hat,t,1)==x1hat + LnC_4(4,:)*(vec_xt-vec_xhat);
diff_eqn5 = diff(x5hat,t,1)==x2hat + LnC_4(5,:)*(vec_xt-vec_xhat);
diff_eqn6 = diff(x6hat,t,1)==x3hat + LnC_4(6,:)*(vec_xt-vec_xhat);
diff_eqn7 = diff(x1t,t,1)==  f1_sub2_t;
diff_eqn8 = diff(x2t,t,1)== f2_sub2_t;
diff_eqn9 = diff(x3t,t,1) == f3_sub2_t;
diff_eqn10= diff(x4t,t,1) == x1t;
diff_eqn11 = diff(x5t,t,1) == x2t;
diff_eqn12 = diff(x6t,t,1) == x3t;

eqns2 = [ diff_eqn1; diff_eqn2; diff_eqn3; diff_eqn4; diff_eqn5; diff_eqn6; diff_eqn7; diff_eqn8; diff_eqn9; diff_eqn10; diff_eqn11;diff_eqn12]; 
[M3,F3] = massMatrixForm(eqns2,state_vec1);
f = M3\F3;
ode_fun = odeFunction(f,state_vec1);
x_init4 =[zeros(6,1); x_init];
[t,x_out] = ode45(ode_fun,[0 sim_time],x_init4);
figure;
title('For Fourth output vector(non-linear)')
hold on;
plot(t,x_out(:,10),'r');
plot(t,x_out(:,4),'--','Color',[1,0,0]);
plot(t,x_out(:,11),'g');
plot(t,x_out(:,5),'--','Color',[0,1,0]);
plot(t,x_out(:,12),'b');
plot(t,x_out(:,6),'--','Color',[0,0,1]);
legend('xc','xchat','t1','t1hat','t2','t2hat')

%plot(t,x_out(:,12),'r');

