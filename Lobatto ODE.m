%MAIN------------------------------------------------------------------

close all;
clear; %%%clear all previous data

%TASK_1
%2y'' + 9y' + 4y = 0 for t belongs to [0,10], y (0) = 5 and y'(0) = -6
%y''(0) = 17
%implicit Lobatto IIIA order 4 method

yInitial = [5; -6];
A = [0, 1; -2, -9/2];
h = 0.01;
T = [0:0.01:10];

%Creating left matrix
LeftMatrix = zeros(6,6);
LeftMatrix(1:2,1:2) = eye(2);
LeftMatrix(3:4,1:2) = -5/24*h*A;
LeftMatrix(3:4,3:4) = eye(2)-1/3*h*A;
LeftMatrix(3:4,5:6) = 1/24*h*A;
LeftMatrix(5:6,1:2) = -1/6*h*A;
LeftMatrix(5:6,3:4) = -2/3*h*A;
LeftMatrix(5:6,5:6) = eye(2)-1/6*h*A;

Un = generateF123(LeftMatrix, A, yInitial, h);
figure('Name','Lobatto IIIA');
plot(T,Un(:,1));
hold on;
title('Lobatto IIIA');
xlabel('t');
ylabel('f(t)');
legend('f(t)');

%ode113
figure('Name','ODE113');
[t,y] = ode113(@myODE, [0,10], transpose(yInitial));
plot(t,y(:,1));
hold on;
title('ODE113');
xlabel('t');
ylabel('f(t)');
legend('f(t)');

%ode113 scatter 
figure('Name','ODE113');
options1 = odeset('RelTol', eps, 'AbsTol', eps);
[t1,y1] = ode113(@myODE, [0:0.01:10], transpose(yInitial), options1);
scatter(t1,y1(:,1));
hold on;
title('ODE113 scattered with options(RelTol, AbsTol)');
xlabel('t');
ylabel('f(t)');
legend('f(t)');

%TASK_2

iterations = 10000;

stepv2 = logspace(-4, -1, iterations); %xInitial, xFinal, numberOfSteps

twoNorm = zeros(iterations , 1);
infNorm = zeros(iterations , 1);
EtwoNorm = zeros(iterations , 1);
EinfNorm = zeros(iterations , 1);

%Generating norms
for i=1:iterations 
    steph = stepv2(i);
    %steph = h*i;
    columns = floor((10/steph)+1);
    [t,yn] = ode113(@myODE, linspace(0, (columns-1)*steph, columns), transpose(yInitial), options1);
    Unx = generateF123(LeftMatrix, A, yInitial, steph);
    twoNorm(i) = norm(Unx(:,1) - yn(:,1))/norm(yn(:,1));
    infNorm(i) = norm(Unx(:,1) - yn(:,1), inf)/norm(yn(:,1), inf);
   
    %%TASK_3
    %Euler implicit Method
    LeftMatrix2 = eye(2) - steph*A;
    UnE = zeros(columns,2);
    UnE(1,:) = yInitial;
    previous2 = yInitial;
    for j=2:columns  
        RightMatrix2 = A*previous2;
        f1 = LeftMatrix2\RightMatrix2;
        UnE(j,:) = UnE(j-1,:)+steph*transpose(f1);
        previous2 = transpose(UnE(j,:));
    end 
    EtwoNorm(i) = norm(UnE(:,1) - yn(:,1))/norm(yn(:,1));
    EinfNorm(i) = norm(UnE(:,1) - yn(:,1), inf)/norm(yn(:,1), inf);
%     if i==1
%     figure('Name','Euler explicit method');
%     plot(T, UnE(:,1));
%     hold on;
%     xlabel('t');
%     ylabel('Value');
%     title('Euler explicit method');
%     legend('f(t)', 'df(t)/dt');
%     end

    %%Euler explicit method
    UnEE = zeros(columns,2);
    UnEE(1,:) = yInitial;
    for k=2:columns
        UnEE(k,:) = UnEE(k-1,:)+steph*transpose(A*transpose(UnEE(k-1,:)));
    end
    EEtwoNorm(i) = norm(UnEE(:,1) - yn(:,1))/norm(yn(:,1));
    EEinfNorm(i) = norm(UnEE(:,1) - yn(:,1), inf)/norm(yn(:,1), inf);
    
%     if i==1
%     figure('Name','Euler explicit method');
%     z = [0:0.001:10];
%     plot(z, UnEE(:,1));
%     hold on;
%     xlabel('t');
%     ylabel('Value');
%     title('Euler explicit method');
%     legend('f(t)');
%     end
    
    
end

figure('Name','Two Norm');
loglog(stepv2, twoNorm);
hold on;
loglog(stepv2, EtwoNorm);
loglog(stepv2, EEtwoNorm);
xlabel('Step h');
ylabel('Value of norm');
title('Two norm vs step[h]');
legend('Two norm Lobatto IIIA', 'Two norm Euler implicit', 'Two norm Euler explicit');

figure('Name','Inf Norm');
loglog(stepv2, infNorm);
hold on;
loglog(stepv2, EinfNorm);
loglog(stepv2, EEinfNorm);
xlabel('Step h');
ylabel('Value of norm');
title('Infinity norm vs step[h]');
legend('Infinity norm Lobatto IIIA', 'Infinity norm Euler implicit', 'Infinity norm Euler explicit');

'Program terminated'
%END_OF_MAIN-----------------------------------------------------------

function dydt = myODE(t,y)
dydt = zeros(2,1);
dydt(1) = y(2);   
dydt(2) = (-9*y(2)-4*y(1))/2;
end

function Un = generateF123(LeftMatrix, A, yInitial, step)
    RightMatrix = zeros(6,1);
    uNminus1 = yInitial;
    columns = floor((10/step)+1);
    Un = zeros(columns,2);
    Un(1,:) = yInitial;
    
    for j=2:columns
        RightMatrix(1:2) = A*uNminus1;
        RightMatrix(3:4) = A*uNminus1;
        RightMatrix(5:6) = A*uNminus1;
        F123 = LeftMatrix\RightMatrix;
        uNminus1 = uNminus1+step*(1/6*F123(1:2)+2/3*F123(3:4)+1/6*F123(5:6));
        Un(j,:) = uNminus1;    
    end  
    
end
