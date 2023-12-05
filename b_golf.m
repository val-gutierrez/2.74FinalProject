function b = b_golf(in1,in2,in3)
%B_GOLF
%    B = B_GOLF(IN1,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    04-Dec-2023 14:38:27

c1 = in3(3,:);
c2 = in3(7,:);
dth1 = in1(3,:);
dth2 = in1(4,:);
g = in3(10,:);
k = in3(11,:);
l1 = in3(4,:);
l2 = in3(8,:);
m1 = in3(1,:);
m2 = in3(5,:);
mE = in3(9,:);
tau1 = in2(1,:);
tau2 = in2(2,:);
th1 = in1(1,:);
th2 = in1(2,:);
t2 = sin(th1);
t3 = sin(th2);
t4 = th1+th2;
t5 = dth1.^2;
t6 = l1.*t2;
t7 = sin(t4);
b = [tau1+dth2.*(c2.*dth1.*l1.*m2.*t3.*2.0+c2.*dth2.*l1.*m2.*t3+dth1.*l1.*l2.*mE.*t3.*2.0+dth2.*l1.*l2.*mE.*t3)+g.*m2.*(t6+c2.*t7)+g.*mE.*(t6+l2.*t7)+c1.*g.*m1.*t2;tau2+dth2.*(c2.*dth1.*l1.*m2.*t3+dth1.*l1.*l2.*mE.*t3)-k.*th2+c2.*g.*m2.*t7+g.*l2.*mE.*t7-c2.*l1.*m2.*t3.*t5-l1.*l2.*mE.*t3.*t5-c2.*dth1.*dth2.*l1.*m2.*t3-dth1.*dth2.*l1.*l2.*mE.*t3];
end
