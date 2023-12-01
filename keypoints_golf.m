function keypoints = keypoints_golf(in1,in2)
%KEYPOINTS_GOLF
%    KEYPOINTS = KEYPOINTS_GOLF(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    30-Nov-2023 21:36:40

c1 = in2(3,:);
c2 = in2(7,:);
l1 = in2(4,:);
l2 = in2(8,:);
th1 = in1(1,:);
th2 = in1(2,:);
t2 = cos(th1);
t3 = sin(th1);
t4 = th1+th2;
t5 = l1.*t2;
t6 = cos(t4);
t7 = l1.*t3;
t8 = sin(t4);
t9 = -t5;
keypoints = reshape([t7,t9,t7+l2.*t8,t9-l2.*t6,c1.*t3,-c1.*t2,t7+c2.*t8,t9-c2.*t6],[2,4]);
end
