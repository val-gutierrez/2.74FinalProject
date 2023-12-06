function E = energy_golf(in1,in2)
%ENERGY_GOLF
%    E = ENERGY_GOLF(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    05-Dec-2023 23:04:12

I1 = in2(2,:);
I2 = in2(6,:);
c1 = in2(3,:);
c2 = in2(7,:);
dth1 = in1(3,:);
dth2 = in1(4,:);
g = in2(10,:);
k = in2(11,:);
l1 = in2(4,:);
l2 = in2(8,:);
m1 = in2(1,:);
m2 = in2(5,:);
mE = in2(9,:);
th1 = in1(1,:);
th2 = in1(2,:);
t2 = cos(th1);
t3 = cos(th2);
t4 = th1+th2;
t5 = c2.^2;
t6 = dth1.^2;
t7 = dth2.^2;
t8 = l1.^2;
t9 = l2.^2;
t10 = l1.*t2;
t11 = cos(t4);
E = (I1.*t6)./2.0+(I2.*t6)./2.0+(I2.*t7)./2.0+(k.*th2.^2)./2.0+(c1.^2.*m1.*t6)./2.0+g.*m2.*(t10+c2.*t11)+g.*mE.*(t10+l2.*t11)+I2.*dth1.*dth2+(m2.*t5.*t6)./2.0+(m2.*t5.*t7)./2.0+(m2.*t6.*t8)./2.0+(mE.*t6.*t8)./2.0+(mE.*t6.*t9)./2.0+(mE.*t7.*t9)./2.0+dth1.*dth2.*m2.*t5+dth1.*dth2.*mE.*t9+c1.*g.*m1.*t2+c2.*l1.*m2.*t3.*t6+l1.*l2.*mE.*t3.*t6+c2.*dth1.*dth2.*l1.*m2.*t3+dth1.*dth2.*l1.*l2.*mE.*t3;
end
