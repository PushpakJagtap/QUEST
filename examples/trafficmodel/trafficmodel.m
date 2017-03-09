clc;
clear all;
v = 70; 
l = 0.25;
T = 0.002777; %T is in hrs = 10sec
q = 0.25;
filename = 'controller.txt';
p=importdata(filename)
x_init=[1.41783,4.92788,10.6711,9.53216,14.5308];
z=x_init;
time = 0;
Y(:,1)=z;
Time(1,1)=time;
for i=1:length(p)
    %Mapping to system inputs
    if(p(i)==0)
      b = [6,0,8,0,0];
    end
    if(p(i)==1)
      b = [6,0,0,0,0];
    end
    if(p(i)==2)
      b = [0,0,8,0,0];
    end
    % traffic model dynamics 
    x(1) = (1-T*v/l)*z(1) + b(1);
    x(2) = T*v/l*z(1) + (1-T*v/l-q) *z(2) + b(2);
    x(3) = T*v/l*z(2) + (1-T*v/l) *z(3) + b(3);
    x(4) = T*v/l*z(3) + (1-T*v/l) *z(4) + b(4);
    x(5) = T*v/l*z(4) + (1-T*v/l) *z(5) + b(5);
    z=x;
    time=time+10;
    Y(:,i+1)=z;
    Time(1,i+1)=time;
end
plot(Time,Y);


