% Get the data of the solution of 2D Incompressible NS Equations with SIMPLE 
% and show as a figure. the data: "solution.dat" = Data file containing for  
% each grid point, coordinate, velocity and pressure, in the following format:
%
%              x(1)  y(1)     u(1)  v(1)     p(1)
%              x(2)  y(2)     u(2)  v(2)     p(2)
%               .     .         .    .        .
%               .     .         .    .        .
%              x(n)  y(n)     u(n)  v(n)     p(n)
%
% Copyrhigt by Zhou Lvwen: zhou.lv.wen[at]gmail.com 

dx = 0.02; dy = 0.02;
solution = load('solution.dat');
x = solution(:,1);
y = solution(:,2);
u = solution(:,3);
v = solution(:,4);
p = solution(:,5);
I = round(x/dx)+1;
J = round(y/dy)+1;

imin=min(I); imax = max(I);
jmin=min(J); jmax = max(J);

U =  zeros(jmax-jmin+1, imax-imin+1);
V =  zeros(jmax-jmin+1, imax-imin+1);
P =  inf*ones(jmax-jmin+1, imax-imin+1);

for k = 1:length(x)
    i = I(k); j = J(k);
    U(j,i) = u(k); V(j,i) = v(k); P(j,i) = p(k);
end

U =  flipud(U); V = -flipud(V); P =  flipud(P);
[X,Y] = meshgrid([[imin:imax]-1]*dx, [[jmin:jmax]-1]*dy);
imagesc(X(:),Y(:),P)
hold on
streamslice(X,Y,U,V,10,'b');
xlabel('x');ylabel('y')
axis image

