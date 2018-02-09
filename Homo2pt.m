close all
clear all

borneinf=-300;
bornesup=300;
altitude=-5;
nbpoints=500;

index=1;
bruitI=0;
bruit=0;
pourcentbruit=0;
bruitnormale=0;
%for bruitI=0:0.01:0.1,
%for pourcentbruit=0:1:25,
%for bruitnormale=0.05:0.05:1,    
   % bruit=abs(altitude*pourcentbruit/100);
    
 %   for k=1:1:1000,
     k=1   
        %Points 3D
        P=[borneinf + (bornesup-borneinf).*rand(nbpoints,1)+bruit.*rand(nbpoints,1),borneinf + (bornesup-borneinf).*rand(nbpoints,1)+bruit.*rand(nbpoints,1),altitude*ones(nbpoints,1)+bruit.*rand(nbpoints,1),ones(nbpoints,1)];
        %for k=1:1:1000,
        %param?tres cam?ra 1
        tx1=10;
        ty1=-5;
        tz1=23;
        roulis1=-90;
        tangage1=10;
        lacet1=15;
        
        %focale
        f=1;
        %f=820;
        K=[f 0 0;0 f 0;0 0 1];
        %param?tres cam?ra 2
        %tx2=24;
        tx2=30*randn(1);
        %ty2=15;
        ty2=30*randn(1);
        %tz2=43;
        tz2=43;%30*randn(1);
        roulis2=10+(5*randn(1));
        tangage2=20+(5*randn(1));
        lacet2=-40+(5*randn(1));
        
        
        %Rotation camŽra 1
        Rt=[cos(deg2rad(tangage1)) 0 sin(deg2rad(tangage1));0 1 0;-sin(deg2rad(tangage1)) 0 cos(deg2rad(tangage1))];
        Rl=[cos(deg2rad(lacet1)) -sin(deg2rad(lacet1)) 0;sin(deg2rad(lacet1)) cos(deg2rad(lacet1)) 0;0 0 1];
        Rr=[1 0 0;0 cos(deg2rad(roulis1)) -sin(deg2rad(roulis1));0 sin(deg2rad(roulis1)) cos(deg2rad(roulis1))];
        R1=Rl*Rt*Rr;
        d1=R1'*[0 0 1]';
        axe1=R1'*[1 0 0]';
        %Points Sphere 1
        M1=K*[R1',-R1'*[tx1 ty1 tz1]'];
        %p1=(K*((R1*P')-(R1*[tx1*ones(nbpoints,1) ty1*ones(nbpoints,1) tz1*ones(nbpoints,1)]')))';
        p1=(M1*P')';
        
        %Rotation camŽra 2
        Rt2=[cos(deg2rad(tangage2)) 0 sin(deg2rad(tangage2));0 1 0;-sin(deg2rad(tangage2)) 0 cos(deg2rad(tangage2))];
        Rl2=[cos(deg2rad(lacet2)) -sin(deg2rad(lacet2)) 0;sin(deg2rad(lacet2)) cos(deg2rad(lacet2)) 0;0 0 1];
        Rr2=[1 0 0;0 cos(deg2rad(roulis2)) -sin(deg2rad(roulis2));0 sin(deg2rad(roulis2)) cos(deg2rad(roulis2))];
        R2=Rl2*Rt2*Rr2;
        d2=R2'*[0 0 1]';
        axe2=R2'*[1 0 0]';
        %Points Sphere 2
        M2=K*[R2',-R2'*[tx2 ty2 tz2]'];
        p2=(M2*P')';
        
        %Projection des points 3D sur l'image normalisée ou la sphère
        for i=1:1:nbpoints,
            %version plan normalisé
            %q1(i,1:3)=p1(i,1:3)/p1(i,3);
            %q2(i,1:3)=p2(i,1:3)/p2(i,3);
            %version sphère
            q1(i,1:3)=p1(i,1:3)/norm(p1(i,:));
            q2(i,1:3)=p2(i,1:3)/norm(p2(i,:));
            
        end
        
        %Ajout du bruit
        Bruit(:,1)=bruitI*randn(nbpoints,1);
        Bruit(:,2)=bruitI*randn(nbpoints,1);
        Bruit(:,3)=bruitI*randn(nbpoints,1);
        Bruit(:,4)=bruitI*randn(nbpoints,1);
        
        q1(:,1:2)=q1(:,1:2)+Bruit(:,1:2);
        q2(:,1:2)=q2(:,1:2)+Bruit(:,3:4);
        
        %q1=(inv(K)*q1')';
        %q2=(inv(K)*q2')';
        
        %for i=1:1:nbpoints,
            %version plan normalisé
            %q1(i,1:3)=q1(i,1:3)/q1(i,3);
            %q2(i,1:3)=q2(i,1:3)/q2(i,3);
        %end
        
        q1sauv=q1;
        q2sauv=q2;

        %Ajout de bruit aux normales d1 et d2
        [Azd1,Eld1,Rd1]=cart2sph(d1(1),d1(2),d1(3));
        Azd1=rad2deg(Azd1);
        Eld1=rad2deg(Eld1);
        bruitAzd1=bruitnormale*randn(1000,1);
        bruitEld1=bruitnormale*randn(1000,1);
        [Azd2,Eld2,Rd2]=cart2sph(d2(1),d2(2),d2(3));
        Azd2=rad2deg(Azd2);
        Eld2=rad2deg(Eld2);
        bruitAzd2=bruitnormale*randn(1000,1);
        bruitEld2=bruitnormale*randn(1000,1);
        
        
            [d1(1),d1(2),d1(3)]=sph2cart(deg2rad(Azd1+bruitAzd1(k)),deg2rad(Eld1+bruitEld1(k)),1);
            [d2(1),d2(2),d2(3)]=sph2cart(deg2rad(Azd2+bruitAzd2(k)),deg2rad(Eld2+bruitEld2(k)),1);
            
        %Dérotation des points sur la sphère
        figure
plot3(q1(:,1)',q1(:,2)',q1(:,3)','*b');
hold on;
line([0 0],[0 0],[0 2]);
%line([0 2*d1(1)],[0 2*d1(2)],[0 2*d1(3)]);
plot3(2,0,0,'r*');
plot3(0,2,0,'g*');
%plot3(P(:,1),P(:,2),P(:,3),'r*');        
        
        [q1,q2,Run,Rdeux]=derotate(q1sauv',q2sauv',d1,d2);
        q1=q1';
        q2=q2';
        
        
        figure
plot3(q1(:,1)',q1(:,2)',q1(:,3)','*b');
hold on;
line([0 0],[0 0],[0 2]);
plot3(2,0,0,'r*');
plot3(0,2,0,'g*');
        %Ajout de points bruités
        % Ptsbruit1=rand(100,3);
        % Ptsbruit2=rand(100,3);
        % for i=1:1:size(Ptsbruit1,1)
        % q1(2*i,:)=Ptsbruit1(i,:)/norm(Ptsbruit1(i,:));
        % q2(2*i,:)=Ptsbruit2(i,:)/norm(Ptsbruit2(i,:));
        % end
        
        
        %Version avec RANSAC
        [Hpinv, inliers] = ransacfithomography2pt(q1', q2', 0.1);
        
        
        %Estimation de l'homographie à partir de deux points
        
%         for i=1:1:2,
%             W(2*i-1,:)=[-q1(i,2)*q2(i,3) -q1(i,1)*q2(i,3) 0 -q1(i,3)*q2(i,3)];
%             W(2*i,:)=[q1(i,1)*q2(i,3) -q1(i,2)*q2(i,3) q1(i,3)*q2(i,3) 0];
%             Q(2*i-1,1)= -q1(i,3)*q2(i,2);
%             Q(2*i,1)= q1(i,3)*q2(i,1);
%         end
%         h=pinv(W)*Q;
%         Hpinv=[h(1) -h(2) h(3);h(2) h(1) h(4);0 0 1];
%         
%         Hpinv=inv(Hpinv);
        
        %Estimation du lacet
        
        Lacet=atan2(Hpinv(2,1),Hpinv(1,1));
        
        
        %Estimation de l'élévation de la translation
        Azimuth=atan2(Hpinv(2,3),Hpinv(1,3));
        
        
        SEalpha=-((sin(Lacet)/Hpinv(2,1))-1);
        CEalpha=(-sin(atan(Hpinv(2,1)/Hpinv(1,1)))/Hpinv(2,1)*Hpinv(1,3))/cos(Azimuth);
        Elevation=atan2(SEalpha,CEalpha)
        SEalpha=Hpinv(2,1)-sin(Lacet);
        CEalpha=-Hpinv(1,3)*(sin(Lacet)/cos(Azimuth));
        Elevation=atan2(SEalpha,CEalpha)
        
        
        [A,E,R]=cart2sph(tx2,ty2,tz2);
        lacet2-lacet1;
        %resul(bruit,k,1)=rad2deg(Lacet);
        %Azimuth
        %Elevation
        
        [TX,TY,TZ]=sph2cart(Azimuth,Elevation,1);
        TCORR=R1*Run'*[TX;TY;TZ];
        [Acorr,Ecorr,Rcorr]=cart2sph(TCORR(1),TCORR(2),TCORR(3));
        
        A;
        %resul(bruit,k,2)=Acorr;
        E;
        %resul(bruit,k,3)=Ecorr;
        %k
        %end
        %bruit
        %end
        
        TCORR=TCORR/TCORR(3);
        %TCORR=TCORR*tz2;
        TCORR=TCORR/norm(TCORR);
        TREF=[tx2-tx1;ty2-ty1;tz2-tz1];
        TREF=TREF/norm(TREF);
        erreur(k,index) = min(real(acos(dot(TREF,TCORR))),abs(pi-abs(real(acos(dot(TREF,TCORR))))));
        erreur_angle(k,index)=real(acos(dot(axe1,Run'*[cos(Lacet) -sin(Lacet) 0;sin(Lacet) cos(Lacet) 0;0 0 1]*Rdeux*axe2)));
        %[real(acos(dot(TREF,TCORR))),abs(pi-real(acos(dot(TREF,TCORR))))]
    %end
    index=index+1;
%end

mean(erreur);
std(erreur);

%plot(mean(erreur),'b');
%hold on;
%plot(mean(erreur)+std(erreur),'*b');
