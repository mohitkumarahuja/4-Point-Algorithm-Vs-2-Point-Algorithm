function [Q,Qp,R,Rp] = derotate(K,Kp,g,gp)

% derotate directional constraints:
c = [0,0,1]';
c = c/norm(c);

a = cross(g,c);
if norm(a) > 0.000001
    a = a/norm(a);
    %theta = acos(dot(g,c))+pi/2
    theta = acos(dot(g,c));
    R = expm(skew3(a)*theta);
    Q = R*K;
    %Q = Q./repmat(Q(3,:), 3,1);
else
    R = eye(3);
    Q = K;
end

ap = cross(gp,c);
if norm(ap) > 0.000001
    ap = ap/norm(ap);
    %thetap = acos(dot(gp,c))+pi/2
    thetap = acos(dot(gp,c));
    Rp = expm(skew3(ap)*thetap);
    Qp = Rp*Kp;
    %Qp = Qp./repmat(Qp(3,:), 3,1);
else
    Rp = eye(3);
    Qp = Kp;
end
