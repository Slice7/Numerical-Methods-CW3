%%
clear all, close all, clc, format long, format compact
%%  
% *Exercise 1*
% 
% (a)

y=@(t)tanh(t);
f=@(t,y)1-y^2;
tmax=20;
y0=0;
Nvec=[19, 21, 40, 80, 160];
xvals=linspace(0, 20, 200);
e_fe=[];
e_fe_end=[];
hold on
for N=Nvec
    [t_fe, u_fe] = feuler(f,[0,tmax],y0,N);
    err=[];
    for i=1:N+1
        err=[err u_fe(i)-y((i-1)*20/N)];
    end
    e_fe=[e_fe max(abs(err))];
    e_fe_end=[e_fe_end abs(u_fe(N+1)-y(20))];
    plot(t_fe, u_fe)
end

e_fe=e_fe'
e_fe_end=e_fe_end'

plot(xvals, y(xvals))
xlabel('t')
ylabel('y')
legend('19', '21', '40', '80', '160', 'tanh(x)', 'Location', 'SouthEast')
title('Forward Euler')
%%
% (b)

e_he=[];
e_he_end=[];
figure
hold on
for N=Nvec
    [t_he, u_he] = heun(f,[0,tmax],y0,N);
    err=[];
    for i=1:N+1
        err=[err u_he(i)-y((i-1)*20/N)];
    end
    e_he=[e_he max(abs(err))];
    e_he_end=[e_he_end abs(u_he(N+1)-y(20))];
    plot(t_he, u_he)
end

e_he=e_he'
e_he_end=e_he_end'

plot(xvals, y(xvals))
xlabel('t')
ylabel('y')
legend('19', '21', '40', '80', '160', 'tanh(x)', 'Location', 'SouthEast')
title('Heun')
%%
% (c)

hvec=1./Nvec;
figure
loglog(hvec, e_fe, hvec, e_he)
xlabel('h')
ylabel('Error')
legend('Forward Euler error', 'Heun error')
title('Error')

% Finding approx. conv. orders
fe_order = (log(e_fe(1))-log(e_fe(5)))/(log(hvec(1))-log(hvec(5)))
he_order = (log(e_he(1))-log(e_he(5)))/(log(hvec(1))-log(hvec(5)))

% From the lectures, we know that Forward Euler has linear convergence and Heun has quadratic convergence. Therefore, we expect to get order 1 for Forward Euler and 2 for Heun which is what we are observing from the graphs.
%%
% (d)

% Every approximation apart from N=19 for Forward Euler and N=19 for Heun seems to reproduce the same asymptotic behaviour.

% y(20) is approximated well for those not mentioned above. When N=21 in Heun, y(20) isn't as good as the rest but seems to be close enough.

% The approximation does not converge to the exact solution when N=19 (i.e. h>1) for both methods, but they do converge for N=21 (i.e h<1) for both methods.

% For Forward Euler, the oscillations in the approximation get larger as tn gets larger. However for Heun, the approximation simply converges to a different solution.
