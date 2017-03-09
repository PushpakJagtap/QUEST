clc;
clear all;
sDIM=10;
a=5e-2;
ae2=5e-3;ae5=5e-3; 
ae1=3.3e-3;ae3=3.3e-3;ae4=3.3e-3;ae6=3.3e-3;ae7=3.3e-3; ae8=3.3e-3; ae9=3.3e-3;ae10=3.3e-3;   
ah=3.6e-3;
te=12;
th=100;
Dt=5; %time step for ode solving
X1=zeros(sDIM,1);
%control inputs
filename = 'controller.txt';
input=importdata(filename);

for i=1:length(input)
   u(1,i)=bitand(floor(input(i)/2),1);
   u(2,i)=bitand(input(i),1);
end
% initial conditions
X_init=[18.9 19 19.1 19.5 20.8 19.7 19 19.8 19 19.8];

%initialization
X=X_init;
time = 0;
k=0;
for j = 1:length(input)*5
    k=ceil(j/5);
    X1(1)=X(1)+Dt*((-a-ae1)*X(1)+a*X(2)+ae1*te);
    X1(2)=X(2)+Dt*((-4*a-ae2-ah*u(1,k))*X(2)+a*X(1)+a*X(7)+a*X(9)+a*X(3)+ae2*te+ah*th*u(1,k));
    X1(3)=X(3)+Dt*((-2*a-ae3)*X(3)+a*X(2)+a*X(4)+ae3*te);
    X1(4)=X(4)+Dt*((-2*a-ae4)*X(4)+a*X(3)+a*X(5)+ae4*te);
    X1(5)=X(5)+Dt*((-4*a-ae5-ah*u(2,k))*X(5)+a*X(4)+a*X(8)+a*X(6)+a*X(10)+ae5*te+ah*th*u(2,k));
    X1(6)=X(6)+Dt*((-a-ae6)*X(6)+a*X(5)+ae6*te);
    X1(7)=X(7)+Dt*((-a-ae7)*X(7)+a*X(2)+ae7*te);
    X1(8)=X(8)+Dt*((-a-ae8)*X(8)+a*X(5)+ae8*te);
    X1(9)=X(9)+Dt*((-a-ae9)*X(9)+a*X(2)+ae9*te);
    X1(10)=X(10)+Dt*((-a-ae10)*X(10)+a*X(5)+ae10*te);
       
    X=X1;
    Y(1,j)=X(1);
    Y(2,j)=X(2);
    Y(3,j)=X(3);
    Y(4,j)=X(4);
    Y(5,j)=X(5);
    Y(6,j)=X(6);
    Y(7,j)=X(7);
    Y(8,j)=X(8);
    Y(9,j)=X(9);
    Y(10,j)=X(10);
    Time(1,j)=time;
    time=time+Dt;     
end
time1=0;
for l = 1:length(input)*500
    m=ceil(l/500);
    Uem(1,l)=u(1,m);
    Uem(2,l)=u(2,m);
    time1=time1+Dt/100;
    Time1(1,l)=time1;
end
figure(1);
plot(Time,Y);
title('Output Response');
figure(2);
plot(Time1,Uem(1,:));
title('control input to heater 1');
figure(3);
plot(Time1,Uem(2,:));
title('control input to heater 2');

