%%%%% Testing the 5 point algorithm
%%%%% Last modified 2010-09-10
%%%%% Nister algorithm is different from his paper

% % clc;
% % clear all;
% % close all;
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%% Structure and Motion data

Radi = 1;
CamDist = 10;
nView = 3;
nPt3D = 10;

RDcmOrig = zeros(3,3,nView);
RDcmOrig(:,:,1) = eye(3);
RDcmOrig(:,:,2) = Func_RotX(10*pi/180);
RDcmOrig(:,:,3) = Func_RotY(15*pi/180);
RQuatOrig = Func_DcmToQuat(RDcmOrig)';

CamOrig = [[0 0 0]' [0.5 0.5 0]' [-0.5 -0.5 0]'];
TOrig = zeros(3,nView);
for iV = 1:nView
    TOrig(:,iV) = -RDcmOrig(:,:,iV)*CamOrig(:,iV);
end

zPt = -Radi + 2*Radi*rand(1,nPt3D);
tPt = 2*pi*rand(1,nPt3D);
rPt = sqrt(Radi^2-zPt.^2);
xPt = rPt.*cos(tPt);
yPt = rPt.*sin(tPt);

Pt3D = [xPt;yPt;zPt+CamDist];

PtSph = zeros(3,nPt3D,nView);
PtSphRAP = zeros(3,nPt3D,nView);
PtImg = [zeros(2,nPt3D,nView);ones(1,nPt3D,nView)];

for iV = 1:nView
    PtSph(:,:,iV) = normc(RDcmOrig(:,:,iV)*Pt3D+TOrig(:,iV)*ones(1,nPt3D));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vl = 1;
vr = 2;

Pt1 = PtSph(:,:,vl);
Pt2 = PtSph(:,:,vr);

solutionE = calibrated_fivepoint(Pt1(:,1:end),Pt2(:,1:end));
%%% each column of matE forms 3 columns (not 3 rows) of an essential matrix

tempVect = zeros(1,size(solutionE,2));

for iC = 1:size(solutionE,2)
    tempVect(iC) = mean(abs((diag(Pt1'*reshape(solutionE(:,iC),3,3)*Pt2))));
end

[temp,EInd] = min(tempVect);

E = reshape(solutionE(:,EInd),3,3);

[uE,dE,vE] = svd(E);

if det(uE)<0
    uE = -uE;
end
if det(vE)<0
    vE = -vE;
end

%%% Verify the cheirality

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Hartley

W = [0 -1 0;1 0 0;0 0 1];

P = zeros(3,4,4);

P(:,:,1) = [(uE*W*vE')' -(uE*W*vE')'*uE(:,end)];
P(:,:,2) = [(uE*W*vE')' (uE*W*vE')'*uE(:,end)];
P(:,:,3) = [(uE*W'*vE')' -(uE*W'*vE')'*uE(:,end)];
P(:,:,4) = [(uE*W'*vE')' (uE*W'*vE')'*uE(:,end)];

flagRT = 0;  %%%%% check if [R T] can be recovered

for i = 1:4
    matD = [eye(3) [0 0 0]' -Pt1(:,end) [0 0 0]';P(:,:,i) [0 0 0]' -Pt2(:,end)];
    [UD,DD,VD] = svd(matD);

    if VD(5,end)/VD(6,end)>0 && VD(3,end)/VD(4,end)>0
        %%%%% VD(5,end) and VD(6,end) should be positive as they are the
        %%%%% scale of the 3D point to the spherical points
        %%%%% use VD(5,end)/VD(6,end) due to up-to-scale solution
        %%%%% VD(3,end)/VD(4,end) give Z of the 3D point
        matOut = P(:,:,i);
        flagRT = 1;
        break;
    end
end

if flagRT == 0
    fprintf('Can not recover [R T] from 5 point algorithm \n');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Nister
% % % D = [0 1 0;-1 0 0;0 0 1];
% % % 
% % % PA = [(uE*D*vE')' -(uE*D*vE')'*normc(uE(:,end))];
% % % 
% % % Ht = [eye(3) [0 0 0]';-2*vE(:,end)' -1];
% % % Hr = diag([1 1 1 -1]);
% % % 
% % % matD = [eye(3) [0 0 0]' -Pt1(:,end) [0 0 0]';...
% % %     PA [0 0 0]' -Pt2(:,end)];
% % % 
% % % [UD,DD,VD] = svd(matD);
% % % 
% % % 
% % % c1 = VD(3,end)*VD(4,end);
% % % c2 = [0 0 1]*(PA*VD(1:4,end))*VD(4,end);
% % % c3 = VD(3,end)*[0 0 0 1]*[eye(3) [0 0 0]';-2*vE(:,end)' -1]*VD(1:4,end);
% % % 
% % % if c1>0 && c2>0
% % %     matOut = PA;
% % % elseif c1<0 && c2<0
% % %     matOut = PA*Hr;
% % % elseif c1*c2<0 && c3>0
% % %     matOut = [(uE*D'*vE')' (uE*D'*vE')'*normc(uE(:,end))];
% % % elseif c1*c2<0 && c3<0
% % %     matOut = [(uE*D'*vE')' -(uE*D'*vE')'*normc(uE(:,end))];
% % % end
% % % 
% % % fprintf('*********** \n iC %d c1 %f c2 %f c3 %f \n\n',iC,c1,c2,c3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Rotation angle error');
disp(rad2deg(abs(2*(acos(Func_DcmToQuat(matOut(1:3,1:3))*[1 0 0 0]')-acos(RQuatOrig(1,vr))))));
disp('Rotation axis error');
disp(rad2deg(acos(dot(normc([[0 0 0]' eye(3)]*Func_DcmToQuat(matOut(1:3,1:3))'),normc(RQuatOrig(2:4,vr))))));
disp('Translation direction error');
disp(rad2deg(acos(dot(normc(matOut(1:3,4)),normc(TOrig(:,vr))))));
