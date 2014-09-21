% Get the data of the 1D shock tube solutions and show as a 
% figure. the data: "solution.dat" = Data file containing for  
% each grid point, coordinate, density, velocity, pressure, 
% in the following format:
%
%              x(1)        rho(1)        u(1)        p(1)
%              x(2)        rho(2)        u(2)        p(2)
%               .            .            .           .
%               .            .            .           .
%              x(n)        rho(n)        u(n)        p(n)


M = load('solution.dat');
x   = M(:,1);
rho = M(:,2);
u   = M(:,3);
p   = M(:,4);

plot(x, rho, '-', x, u, '--', x, p, '-.', 'linewidth',1.5)
xlabel('$x$', 'interpreter','latex','fontsize', 15)
ylabel('$\rho$, $u$, $p$', 'interpreter','latex','fontsize', 15)
legend('density', 'volecity', 'pressure',0)
