%%
clear all, close all, clc, format long, format compact
%%
% *Exercise 2*
%
% (a)

N=40;
h=1./N;
k=h;
T=0.2;
M=T/k;
I=eye(N-1);
xvals=linspace(0, 1, N+1);

u=@(x,t)exp((-pi^2)*t)*sin(pi*x)+(x/2).*(1-x);  % Exact solution

A=(2/h^2)*diag(ones(N-1,1)) - (1/h^2)*diag(ones(N-2,1),1) - (1/h^2)*diag(ones(N-2,1),-1);
b=ones((N-1),1);

ux0=zeros(N-1, M+1);
for i=1:N-1
    ux0(i,1)=sin(pi*xvals(i+1))+0.5*xvals(i+1)*(1-xvals(i+1));
end

% Forward Euler
uFE=ux0;
for n=1:M
    uFE(:,n+1)=uFE(:,n)+h*((-A)*uFE(:,n)+b);
end
for i=1:M+1
    plot(xvals,[0;uFE(:,i);0])
    hold on
end
xlabel('x')
ylabel('u')
ylim([0 1.2])
legend('t=0', 't=0.025', 't=0.05', 't=0.075', 't=0.1', 't=0.125', 't=0.15', 't=0.175', 't=0.2')
title('Forward Euler')

% Backward Euler
uBE=ux0;
for n=1:M
    uBE(:,n+1)=(I+h.*A)\(uBE(:,n)+h.*b);
end

figure
for i=1:M+1
    plot(xvals,[0;uBE(:,i);0])
    hold on
end
xlabel('x')
ylabel('u')
legend('t=0', 't=0.025', 't=0.05', 't=0.075', 't=0.1', 't=0.125', 't=0.15', 't=0.175', 't=0.2')
title('Backward Euler')

% Crank Nicolson
uCN=ux0;
for n=1:M
    uCN(:,n+1)=(I+(h/2).*A)\(uCN(:,n)-(h/2).*A*uCN(:,n)+h.*b);
end

figure
for i=1:M+1
    plot(xvals,[0;uCN(:,i);0])
    hold on
end
xlabel('x')
ylabel('u')
legend('t=0', 't=0.025', 't=0.05', 't=0.075', 't=0.1', 't=0.125', 't=0.15', 't=0.175', 't=0.2')
title('Crank-Nicolson')

% Exact soln
tmesh=linspace(0,0.2,M+1);
figure
for i=1:M+1
    plot(xvals,u(xvals,tmesh(i)))
    hold on
end
xlabel('x')
ylabel('u')
legend('t=0', 't=0.025', 't=0.05', 't=0.075', 't=0.1', 't=0.125', 't=0.15', 't=0.175', 't=0.2')
title('Exact solution')

% Crank-Nicolson seems to have the most accuracy, followed by Backward Euler and lastly, Forward Euler. The Forward Euler method performs very poorly when t=0.75 and t=0.2.
%%
% (b)

Nvals=[10,20,40,80,160];
hvals=1./Nvals;

E_FE=[];
E_BE=[];
E_CN=[];
for N=Nvals
    h=1./N;
    M=N*T;
    I=eye(N-1);
    xvals=linspace(0, 1, N+1);

    A=(2/h^2)*diag(ones(N-1,1)) - (1/h^2)*diag(ones(N-2,1),1) - (1/h^2)*diag(ones(N-2,1),-1);
    b=ones((N-1),1);

    ux0=zeros(N-1, M+1);
    for i=1:N-1
        ux0(i,1)=sin(pi*xvals(i+1))+0.5*xvals(i+1)*(1-xvals(i+1));
    end
    
    uFE=ux0;
    uBE=ux0;
    uCN=ux0;
    for n=1:M
        uFE(:,n+1)=uFE(:,n)+h*((-A)*uFE(:,n)+b);
        uBE(:,n+1)=(I+h.*A)\(uBE(:,n)+h.*b);
        uCN(:,n+1)=(I+(h/2).*A)\(uCN(:,n)-(h/2).*A*uCN(:,n)+h.*b);
    end
    
    err=sqrt(h)*norm(uFE(:,M+1)'-u(linspace(h, 1-h, N-1),T));
    E_FE=[E_FE err];
    err=sqrt(h)*norm(uBE(:,M+1)'-u(linspace(h, 1-h, N-1),T));
    E_BE=[E_BE err];
    err=sqrt(h)*norm(uCN(:,M+1)'-u(linspace(h, 1-h, N-1),T));
    E_CN=[E_CN err];
end

figure
loglog(hvals, E_FE)
xlabel('h')
ylabel('error')
title('Forward Euler error')

figure
loglog(hvals, E_BE)
xlabel('h')
ylabel('error')
title('Backward Euler error')

figure
loglog(hvals, E_CN)
xlabel('h')
ylabel('error')
title('Crank-Nicolson error')

% Only Backward Euler and Crank-Nicolson converge so we will study the order of convergence for these two.

BE_order = (log(E_BE(1))-log(E_BE(5)))/(log(hvals(1))-log(hvals(5)))
CN_order = (log(E_CN(1))-log(E_CN(5)))/(log(hvals(1))-log(hvals(5)))
%%
% (c)

% fix k small
Nvals=[10,20,40];
hvals=1./Nvals;

E_BE=[];
E_CN=[];
for N=Nvals
    h=1./N;
    M=6400;
    I=eye(N-1);
    xvals=linspace(0, 1, N+1);

    A=(2/h^2)*diag(ones(N-1,1)) - (1/h^2)*diag(ones(N-2,1),1) - (1/h^2)*diag(ones(N-2,1),-1);
    b=ones((N-1),1);

    ux0=zeros(N-1, M+1);
    for i=1:N-1
        ux0(i,1)=sin(pi*xvals(i+1))+0.5*xvals(i+1)*(1-xvals(i+1));
    end
    uBE=ux0;
    uCN=ux0;
    for n=1:M
        uBE(:,n+1)=(I+h.*A)\(uBE(:,n)+h.*b);
        uCN(:,n+1)=(I+(h/2).*A)\(uCN(:,n)-(h/2).*A*uCN(:,n)+h.*b);
    end
    err=sqrt(h)*norm(uBE(:,M+1)'-u(linspace(h, 1-h, N-1),T));
    E_BE=[E_BE err];
    err=sqrt(h)*norm(uCN(:,M+1)'-u(linspace(h, 1-h, N-1),T));
    E_CN=[E_CN err];
end

figure
loglog(hvals, E_BE)
xlabel('h')
ylabel('error')
title('k small Backward Euler')

figure
loglog(hvals, E_CN)
xlabel('h')
ylabel('error')
title('k small Crank-Nicolson')

BE_order_p = (log(E_BE(1))-log(E_BE(3)))/(log(hvals(1))-log(hvals(3)))
CN_order_p = (log(E_CN(1))-log(E_CN(3)))/(log(hvals(1))-log(hvals(3)))

% fix h small
Mvals=[8,16,32];
kvals=1./Mvals;

E_BE=[];
E_CN=[];
for M=Mvals
    N=640;
    h=1./N;
    I=eye(N-1);

    A=(2/h^2)*diag(ones(N-1,1)) - (1/h^2)*diag(ones(N-2,1),1) - (1/h^2)*diag(ones(N-2,1),-1);
    b=ones((N-1),1);
    xvals=linspace(0, 1, N+1);

    ux0=zeros(N-1, M+1);
    for i=1:N-1
        ux0(i,1)=sin(pi*xvals(i+1))+0.5*xvals(i+1)*(1-xvals(i+1));
    end
    uBE=ux0;
    uCN=ux0;
    for n=1:M
        uBE(:,n+1)=(I+h.*A)\(uBE(:,n)+h.*b);
        uCN(:,n+1)=(I+(h/2).*A)\(uCN(:,n)-(h/2).*A*uCN(:,n)+h.*b);
    end
    err=sqrt(h)*norm(uBE(:,M+1)'-u(linspace(h, 1-h, N-1),T));
    E_BE=[E_BE err];
    err=sqrt(h)*norm(uCN(:,M+1)'-u(linspace(h, 1-h, N-1),T));
    E_CN=[E_CN err];
end

figure
loglog(kvals, E_BE)
xlabel('k')
ylabel('error')
title('h small Backward Euler')

figure
loglog(kvals, E_CN)
xlabel('k')
ylabel('error')
title('h small Crank-Nicolson')

BE_order_q = (log(E_BE(1))-log(E_BE(3)))/(log(kvals(1))-log(kvals(3)))
CN_order_q = (log(E_CN(1))-log(E_CN(3)))/(log(kvals(1))-log(kvals(3)))
%%
% (d)

Nvals=[10,20,40,80,160,320];
eig_min=[];
eig_max=[];
for N=Nvals
    h=1./N;
    A=(2/h^2)*diag(ones(N-1,1)) - (1/h^2)*diag(ones(N-2,1),1) - (1/h^2)*diag(ones(N-2,1),-1);
    eig_min=[eig_min min(eig(A))];
    eig_max=[eig_max max(eig(A))];
end
eig_min=eig_min'
eig_max=eig_max'

% As h tends to 0, the minimum eigenvalue seems to be converging to ~10 and the maximum eigenvalue seems to be tending to infinity.
