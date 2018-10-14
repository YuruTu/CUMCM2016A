clear all;
close all;
clc
format long
syms h S Fw Ff Ff1 a b c d l L F depth n pl m x1 y1 y t distance n a1 b1;
F=[];
theta=[];?
v=24; %风速
l=105*10^(-3); %锚链每节链环的长度
L=22.05; %锚链的总长度
num=0; %通过更改不在海床上的链节数得到一个最优解
num1=round(L/l);
num2=0;
lin=0/180*pi; %第一个链节与水平方向的夹角
lin1=90/180*pi;
lin2=0;
m2=1200; %重物球质量
pg=7.7*10^3; %重物球的密度（单位：kg/m^3)
depth=20; %水深
pl=7; %锚链单位长度的质量
vh=0; %海水流速
g=9.8; %可通过改变此语句来修改重力加速度，单位为m/s^2
p=1.025*10^3; %海水密度
M=1000; %浮标质量
m=10; %钢管质量
m1=100; %设备和钢桶总质量
y=0;
d=1;
j1=0;
j2=0;
while(abs(y-d)>0.005)%在这里选择所需要的精度，
if (y>d)&&(num<round(L/l))
num1=num;
num=round((num1+num2)/2);
elseif (y<d)&&(num<round(L/l));
num2=num;
num=round((num1+num2)/2);
elseif (y<d)&&(num==round(L/l))
lin2=lin;
lin=(lin1+lin2)/2;
elseif(y>d)&&(num==round(L/l))
lin1=lin;
lin=(lin1+lin2)/2;
end
%钢桶受到的浮力
Ff1=p*g*pi*(0.3/2)^2;
%钢管收到的浮力
Ff2=p*g*pi*(0.05/2)^2;
%重物球所受浮力
Ffg=p*g*m2/pg;
%重物球所受海水水流力
Fhg=374*pi*((m2/pg/3/4)^(1/3))^2*vh^2;
%风对浮标受力面的投影面积
S=2*(2-h);
%风对浮标产生的力
Fw=0.625*S*v^2;
%浮标在水中的体积
V=pi*(2/2)^2*h;
%浮标所受到的浮力
Ff=p*g*V;
%浮标受到海水的近似水流力
Fb=374*2*h*vh^2;
%钢桶受到海水的近似水流力
Fs1=374*0.3*vh^2;
%钢管受到海水的水流力的近似值
Fs=374*0.05*vh^2;
%浮标浸没水中的高度
if num==round(L/l)
h=(m2*g+M*g+4*m*g+m1*g-Ff1-4*Ff2-Ffg+pl*L*g+(Fhg+4*Fs+Fs1)*tan(lin))/(p*g*pi-(1.25*v^2+374*vh^2)*tan(lin));
else?
h=(m2*g+M*g+4*m*g+m1*g-Ff1-4*Ff2-Ffg+num*pl*l*g)/(p*g*pi);
end
a=Fw+Fb;
b=-M*g+Ff+(Fw+Fb)*tan(lin);
if j1==0
a=eval(a);
b=eval(b);
else
end
F(1)=sqrt(a^2+b^2);
theta(1)=atan(b/a);
n=0;
for i=1:4
%钢管受到海水的水流力
Fh(i)=374*0.05*sin(theta(i));
n=n+Fh(i);
a=Fw+Fb+n;
if j1==0
a=eval(a);
else
end
b=F(i)*sin(theta(i))+p*g*pi*(50*10^(-3)/2)^2-m*g;
F(i+1)=sqrt(a^2+b^2);
theta(i+1)=atan(b/a);
end
c=0;
for i=1:5
c=c+sin(theta(i));
end
d=depth-c-h;
y1=lin;
distance=0;
if num==round(L/l)
y=l*sin(y1);
x1=Fw/sqrt(1-(sin(y1))^2);
for i=1:num-1
m=(x1*sin(y1)+i*pl*l*g)/sqrt((x1*sin(y1)+i*pl*l*g)^2+Fw^2);
m=m*l;
y=y+m;
n=Fw/sqrt((x1*y1+i*pl*l*g)^2+Fw^2)*l;
if j1==0
n=eval(n);
else
end
distance=distance+n;
if j1==0
y=eval(y);
else
end
end
else
y=y1*l;
distance=(round(L/l)-num)*l;
for i=1:num
x1=Fw/sqrt(1-(sin(y1))^2);
m=(x1*y1+i*pl*l*g)/sqrt((x1*y1+i*pl*l*g)^2+Fw^2)*l;
y=y+m;?
n=Fw/sqrt((x1*y1+i*pl*l*g)^2+Fw^2)*l;
if j1==0
y=eval(y);
n=eval(n);
else
end
distance=distance+n;
end
end
m=0;
j1=1;
j2=j2+1;
end
%钢桶受到的浮力
Ff1=p*g*pi*(0.3/2)^2;
%钢管收到的浮力
Ff2=p*g*pi*(0.05/2)^2;
%重物球所受浮力
Ffg=p*g*m2/pg;
%重物球所受海水水流力
Fhg=374*pi*((m2/pg/3/4)^(1/3))^2*vh^2;
%风对浮标受力面的投影面积
S=2*(2-h);
%风对浮标产生的力
Fw=0.625*S*v^2;
%浮标在水中的体积
V=pi*(2/2)^2*h;
%浮标所受到的浮力
Ff=p*g*V;
%浮标受到海水的近似水流力
Fb=374*2*h*vh^2;
%钢桶受到海水的近似水流力
Fs1=374*0.3*vh^2;
%钢管受到海水的水流力的近似值
Fs=374*0.05*vh^2;
%浮标浸没水中的高度
if num==round(L/l)
h=(m2*g+M*g+4*m*g+m1*g-Ff1-4*Ff2-Ffg+pl*L*g+(Fhg+4*Fs+Fs1)*tan(lin))/(p*g*pi-(1.25*v^2+374*vh^2)*tan(lin));
else?
h=(m2*g+M*g+4*m*g+m1*g-Ff1-4*Ff2-Ffg+num*pl*l*g)/(p*g*pi);
end
a=Fw+Fb;
b=-M*g+Ff+(Fw+Fb)*tan(lin);
F(1)=sqrt(a^2+b^2);
theta(1)=atan(b/a);
n=0;
for i=1:4
%钢管受到海水的水流力
Fh(i)=374*0.05*sin(theta(i));
n=n+Fh(i);
a=Fw+Fb+n;
b=F(i)*sin(theta(i))+p*g*pi*(50*10^(-3)/2)^2-m*g;
F(i+1)=sqrt(a^2+b^2);
theta(i+1)=atan(b/a);
end
disp('输出钢管和钢桶的倾斜角度（角度制）')
th=90-theta*180/pi
m=85*pi/180;
if theta(5)>m
disp('钢桶的倾斜角足够小，测量准确')
else?
disp('钢桶的倾斜角过大')
end
c=0;
for i=1:5
c=c+sin(theta(i));
end
d=depth-c-h;
y1=lin;
distance=0;
if num==round(L/l)
y=l*sin(y1);
x1=Fw/sqrt(1-(sin(y1))^2);
for i=1:num-1
m=(x1*sin(y1)+i*pl*l*g)/sqrt((x1*sin(y1)+i*pl*l*g)^2+Fw^2);
m=m*l;
y=y+m;
n=Fw/sqrt((x1*y1+i*pl*l*g)^2+Fw^2)*l;
distance=distance+n;
plot(distance,y,'o')
hold on
end
else
y=y1*l;
for i=1:round(L/l)-num
distance=i*l;
y=0;
plot(distance,y,'o')
hold on
grid on
end
for i=1:num
x1=Fw/sqrt(1-(sin(y1))^2);
m=(x1*y1+i*pl*l*g)/sqrt((x1*y1+i*pl*l*g)^2+Fw^2)*l;
y=y+m;?
n=Fw/sqrt((x1*y1+i*pl*l*g)^2+Fw^2)*l;
if j1==0
y=eval(y);
n=eval(n);
else
end
distance=distance+n;
plot(distance,y,'o')
hold on
end
end
m=0;
for i=1:5
m=m+cos(theta(i));
end
%浮标的运动半径
disp('输出浮标的运动半径')
ans=distance+m
