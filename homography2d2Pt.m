

function H = homography2d2Pt(varargin)

[x1, x2] = checkargs(varargin(:));

% Attempt to normalise each set of points so that the origin 
% is at centroid and mean distance from origin is sqrt(2).
    [x1, T1] = normalise2dpts(x1);
    [x2, T2] = normalise2dpts(x2);
    

Npts = size(x1,2);

for i = 1:Npts
    W(2*i-1,:)=[-x1(2,i)*x2(3,i) -x1(1,i)*x2(3,i) 0 -x1(3,i)*x2(3,i)];
    W(2*i,:)=[x1(1,i)*x2(3,i) -x1(2,i)*x2(3,i) x1(3,i)*x2(3,i) 0];
    Q(2*i-1,1)= -x1(3,i)*x2(2,i);
    Q(2*i,1)= x1(3,i)*x2(1,i);
end
h=pinv(W)*Q;
Hpinv=[h(1) -h(2) h(3);h(2) h(1) h(4);0 0 1];
H=Hpinv;
% Denormalise
H = T2\H*T1;


function [x1, x2] = checkargs(arg)

if length(arg) == 2
    x1 = arg{1};
    x2 = arg{2};
    if ~all(size(x1)==size(x2))
        error('x1 and x2 must have the same size');
    elseif size(x1,1) ~= 3
        error('x1 and x2 must be 3xN');
    end
elseif length(arg) == 1
    if size(arg{1},1) ~= 6
        error('Single argument x must be 6xN');
    else
        x1 = arg{1}(1:3,:);
        x2 = arg{1}(4:6,:);
    end
else
    error('Wrong number of arguments supplied');
end
