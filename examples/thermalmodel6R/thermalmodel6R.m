clc;
clear all;
sDIM=6;
a=5e-2;
ae1=5e-3;ae4=5e-3; 
ae2=3.3e-3;ae3=3.3e-3;ae5=3.3e-3;ae6=3.3e-3; 
ah=3.6e-3;
te=12;
th=100;
Dt=5;
lb=[0,0]
ub=[1,1]
eta=[0.5,1]
%variable Declaration
Y=zeros(sDIM,6);
X1=zeros(sDIM,1);
%control inputs
filename = 'controller.txt';
input=importdata(filename);
	
for i=1:length(lb)
	nInput(i)=floor(ub(i)/eta(i))-ceil(lb(i)/eta(i))+1;
end
for i=1:length(lb)-1
    temp1=1;
    for j=i:length(lb)-1
        temp1=temp1*nInput(j+1);
    end
    xx(i)=temp1;
end
xx;
for j=1:length(input)
    for i=1:length(lb)
        if i<length(lb)
            temp2=floor(input(j)/xx(i));
            temp3=mod(temp2,xx(i));
            u(i,j)=lb(i)+temp3*eta(i);
        elseif i==length(lb)
            temp4=mod(input(j),xx(i-1));
            u(i,j)=lb(i)+temp4*eta(i);
        end
    end
end               

% initial conditions
X_init=[16 16 16 16 16 16];
%X_init=19.5*ones(sDIM,D_tau);

%initial conditions
X=X_init;
time = 0;
k=0;
for j = 1:length(input)*Dt
    k=ceil(j/5);
    X1(1)=X(1)+Dt*((-3*a-ae1-ah*u(1,k))*X(1)+a*X(2)+a*X(3)+a*X(5)+ae1*te+ah*th*u(1,k));
    X1(2)=X(2)+Dt*((-2*a-ae2)*X(2)+a*X(1)+a*X(4)+ae2*te);
    X1(3)=X(3)+Dt*((-2*a-ae3)*X(3)+a*X(1)+a*X(4)+ae3*te);
    X1(4)=X(4)+Dt*((-3*a-ae4-ah*u(2,k))*X(4)+a*X(2)+a*X(3)+a*X(6)+ae4*te+ah*th*u(2,k));
    X1(5)=X(5)+Dt*((-a-ae5)*X(5)+a*X(1)+ae5*te);
    X1(6)=X(6)+Dt*((-a-ae6)*X(6)+a*X(4)+ae6*te);
    X=X1;  
    Y(:,j)=X1;
    Time(1,j)=time;
    time=time+Dt;   
end
time1=0;
for l = 1:length(input)*Dt*100
    m=ceil(l/500);
    Uem(1,l)=u(1,m);
    Uem(2,l)=u(2,m);
    time1=time1+Dt/100;
    Time1(1,l)=time1;
end
figure(1);
plot(Time,Y);
figure(2);
plot(Y(1,:),Y(2,:));
hold on;
plot([17.5, 17.5 ,21.5,21.5,17.5],[17.5,21.5,21.5,17.5,17.5]);
figure(3);
plot(Y(3,:),Y(4,:));
hold on;
plot([17.5, 17.5 ,21.5,21.5,17.5],[17.5,21.5,21.5,17.5,17.5]);
figure(4);
plot(Y(5,:),Y(6,:));
hold on;
plot([17.5, 17.5 ,21.5,21.5,17.5],[17.5,21.5,21.5,17.5,17.5]);
figure(5);
plot(Time1,Uem(1,:));
figure(6);
plot(Time1,Uem(2,:));

