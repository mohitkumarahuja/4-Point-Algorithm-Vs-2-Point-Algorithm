%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%Motion estimation using 4 points algorithm%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%                 vs                       %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%      2 points algorithm with IMU         %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%        planar scene : homography         %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name :    Mohit Ahuja
%            MSCV - 2
%      University of Burgundy
%        Le Creusot Campus
%
%             TEST-2

close all;
clear all;
clc;

%% Test 2 : Example with different datas
% Propose a test with different camera position
% R1 = I, T1 = 0
% Angles of rotation of camera 2 between 0° and 45°
% Translation of 0 to 100
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('     %%%%%%%%%%%%%%% TEST-2: %%%%%%%%%%%%%%%%');
disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  Data Generation  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lowerbound=-300;
upperbound=300;
nbpoints = 50;
zposition = 500;
P = [lowerbound + (upperbound-lowerbound).*rand(nbpoints,1),...
    lowerbound + (upperbound-lowerbound).*rand(nbpoints,1),...
    zposition*ones(nbpoints,1),...
    ones(nbpoints,1)];

% the data belong on a plane  of equation N'X+d=0 d=-zposition N = [0,0,1]'

% Camera Parameter (Camera is Calibrated)
f=1; u0 = 0; v0 = 0;
K=[f 0 u0;0 f v0;0 0 1];
K1=[f 0 u0 0;0 f v0 0;0 0 1 0];

%%%%%%%%%%%%%%%%%%%%%%%%%   Camera 1   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Camera position at time 1
tx1=10; ty1=-2; tz1=20;
% Rotation angles in degree
roll1=5; pitch1=5; yaw1=30;

% Rotation of the camera 1
Rp1 = [cos(deg2rad(pitch1)) 0 sin(deg2rad(pitch1));...
       0 1 0;...
       -sin(deg2rad(pitch1)) 0 cos(deg2rad(pitch1))];
Ry1 = [cos(deg2rad(yaw1)) -sin(deg2rad(yaw1)) 0;...
       sin(deg2rad(yaw1)) cos(deg2rad(yaw1)) 0;...
       0 0 1];
Rr1 = [1 0 0;...
       0 cos(deg2rad(roll1)) -sin(deg2rad(roll1));...
       0 sin(deg2rad(roll1)) cos(deg2rad(roll1))];
R1 = Ry1*Rp1*Rr1;

T1 = [tx1,ty1,tz1]';

% Camera Image 1 :
for i =1 : nbpoints
    P1(i,:) = K1*[R1' , -R1'*T1; 0 0 0 1]*P(i,:)';
end

%%%%%%%%%%%%%%%%%%%%%%%   Camera 2   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Camera position at time 2
tx2=18; ty2=-8; tz2=80;
% Rotation angles in degree
roll2=5; pitch2=5; yaw2=8;

% Rotation of the camera 2
Rp2=[cos(deg2rad(pitch2)) 0 sin(deg2rad(pitch2));...
     0 1 0;...
     -sin(deg2rad(pitch2)) 0 cos(deg2rad(pitch2))];
Ry2=[cos(deg2rad(yaw2)) -sin(deg2rad(yaw2)) 0;...
     sin(deg2rad(yaw2)) cos(deg2rad(yaw2)) 0;...
     0 0 1];
Rr2=[1 0 0;...
     0 cos(deg2rad(roll2)) -sin(deg2rad(roll2));...
     0 sin(deg2rad(roll2)) cos(deg2rad(roll2))];
R2 =Ry2*Rp2*Rr2;

T2 = [tx2,ty2,tz2]';

% Camera Image 2 :
for i =1 : nbpoints
    P2(i,:) = K1*[R2' , -R2'*T2; 0 0 0 1]*P(i,:)';
end

%%%%%%%%%%%%%%%%%%%%%%%    Display Data     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on
plot3(P(:,1),P(:,2),P(:,3),'*')
a = plot3(tx1,ty1,tz1,'g*');
Rx = R1*[100,0,0]';
line([tx1,tx1+Rx(1)],[ty1,ty1+Rx(2)],[tz1,tz1+Rx(3)],'Color','b')
Ry = R1*[0,100,0]';
line([tx1,tx1+Ry(1)],[ty1,ty1+Ry(2)],[tz1,tz1+Ry(3)],'Color','g')
Rz = R1*[0,0,100]';
line([tx1,tx1+Rz(1)],[ty1,ty1+Rz(2)],[tz1,tz1+Rz(3)],'Color','r')
b = plot3(tx2,ty2,tz2,'r*');
Rx = R2*[100,0,0]';
line([tx2,tx2+Rx(1)],[ty2,ty2+Rx(2)],[tz2,tz2+Rx(3)],'Color','b')
Ry = R2*[0,100,0]';
line([tx2,tx2+Ry(1)],[ty2,ty2+Ry(2)],[tz2,tz2+Ry(3)],'Color','g')
Rz = R2*[0,0,100]';
line([tx2,tx2+Rz(1)],[ty2,ty2+Rz(2)],[tz2,tz2+Rz(3)],'Color','r')
axis equal
legend([a, b], 'Camera 1', 'Camera 2')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%    4 pts algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%% 4 Point Algorithm: %%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(' ');

% Theoretical Homography H :
N = [0,0,1]';
disp('%%%%%    Displacement (D)     %%%%%');
d = (N'*T1 - zposition) 

% Explain why d is expressed like this! 
%
% Solution: 1. As "d" is the Displacement of the camera from the Norm "N" 
%              of the plane, so it is computed as the difference between 
%              the camera's height and the height of Points. 
%              Also, 
%           2. Because with 4-Point Algorithm, the solution of AE = 0 is a
%              2-dimentional space, so we will take only x and y in our
%              case.


% Homography Estimation using Translation and Rotation
Rt = R2'*R1;
Tt = R2'*(T1-T2);
H = Rt-Tt*(R1'*N)'/d;
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('Homography estimation using the Translation and Rotation');
disp('%%%%%    Homography H = Rt-Tt*(R1*N)/d  :  %%%%%');
H = H/H(3,3)

% Verification of Computed Homography
HP1 = H*P1(1,:)';   % As P2(1,:) = H*P1(1,:), So HP1 should be = P2(1,:)
HP1 = HP1/HP1(3)    % As explained above, In 4-Point Algorithm, we only take 
                    % [X, Y] Coordinates. so we are dividing x and y by z.          
disp('%% Ground Truth Should be = HP1 to validate the result %%');
disp(' Which Means, it is satisfying the equation P2(1,:) = H*P1(1,:)');
GroundTruth = (P2(1,:)/P2(1,3))'   % Ground Truth  = P2(1,:)
disp(' ');


% Homography estimation using points in two images Using SVD.
H4pt = homography2d(P1',P2');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('Homography estimation using Points in Two Images');
disp('%%%%% Using SVD,    P2(1,:) = H*P1(1,:)  %%%%%');
H4pt = H4pt/H4pt(3,3)
disp('You can see, the H and H4pt is exactly the same');
% It can be seen the Command Window that using both techniques, the
% Homography Matrix is the same.

% Homography decomposition
% solutions = invhomog(H4pt);
% solutions(1).T(:,4)/solutions(1).T(3,4);
% solutions(1).T(1:3,1:3)
% solutions(2).T(:,4)/solutions(2).T(3,4);
% solutions(2).T(1:3,1:3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%    2 pts algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' '); 
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%% 2 Point Algorithm: %%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

% HYPOTHESIS roll and pitch angles are known :
% Virtual image from camera 2 (Z-axis correpsonds to the vertical)
PV1 = (Rp1*Rr1*P1(:,1:3)')';
PV2 = (Rp2*Rr2*P2(:,1:3)')';
N =[0 0 1]';    % Norm of the Plane

% Homography estimation using Translation and Rotation.
disp(' ');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('Homography estimation using the Translation and Rotation');
disp('%%%%%    Homography H in 2-Point Algo:  %%%%%');
HV2 = Ry2'*Ry1-Ry2'*(T1-T2)/(N'*T1-zposition)*(Ry1'*N)'
disp(' Validation of the Computed Homography: ');
H_Valid = HV2*PV1(2,:)'/norm(HV2*PV1(2,:)')
disp('%% Ground Truth Should be = H_Valid to validate the result %%');
GroundTruth = (PV2(2,:)/norm(PV2(2,:)))'
disp('%% True Homography : %%');
TrueHomography = HV2/HV2(3,3)

% Homography estimation using points in two images Using SVD.
H = homography2d2Pt(PV1',PV2');
disp('%% Here, The Estimated H should be = True Homography to validate our result %%');
Estimated_H = H/H(3,3)

%%%%%%%%%%%%%%% Yaw Estimation (see exercice 1)
disp(' ');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('Yaw Angle Computation:');
Yaw = atan2(H(2,1),H(1,1))           % Computed Yaw Angle
YawGroundTruth = deg2rad(yaw1-yaw2)  % Real Yaw Angle
disp(' Both Yaw and YawGroundTruth should be Equal');

Azimuth=atan2(H(2,3),H(1,3));
SEalpha=-((sin(Yaw)/H(2,1))-1);
CEalpha=(-sin(Yaw)/H(2,1)*H(1,3))/cos(Azimuth);
Elevation=atan2(SEalpha,CEalpha);
T2E = [cos(Azimuth)*cos(Elevation);sin(Azimuth)*cos(Elevation);sin(Elevation)];
T2t = Ry2'*(T1-T2);
% Validating the Translation with the real translation
disp(' ');
disp('%%%%%%%%%%% Validation of the translation %%%%%%%%%%%');
Computed_T = T2t/T2t(3) % Computed T
Original_T = T2E/T2E(3) % Original T
disp(' '); disp(' '); 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%   END of TEST-2  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%