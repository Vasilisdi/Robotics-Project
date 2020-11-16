clear all;clc;close all;
%%
%I Meros
%prwto erwthma
%dedomena megthh kai shmeia A,B,C,D
l3=0.2;l5=0.3;l2=0.1; 
%PA=[0.3 0.4 0.6 0];PB=[0.5 0.6 0.9 degtorad(30)];PC=[0.1 0.5 0.4 degtorad(60)];PD=[0.7 0.2 0.3 degtorad(90)];  %LATHOS POREIA, DEN EINAI AYKSHTIKH
PA=[0.3 0.4 0.6 0];PB=[0.5 0.6 0.9 degtorad(30)];PC=[0.7 0.8 1 degtorad(60)];PD=[1.2 1.2 1.3 degtorad(90)]; %ontws aukshtikh poreia
%vriskw ta qi
d2tA=PA(1,3);d3tA=sqrt(PA(1,1)^2+PA(1,2)^2-l3^2)-l5;theta1A=radtodeg(atan((l5*PA(1,2)+d3tA*PA(1,2)-l3*PA(1,1))/(l5*PA(1,1)+d3tA*PA(1,1)+l3*PA(1,2))));qA=[theta1A d2tA d3tA PA(1,4)];
d2tB=PB(1,3);d3tB=sqrt(PB(1,1)^2+PB(1,2)^2-l3^2)-l5;theta1B=radtodeg(atan((l5*PB(1,2)+d3tB*PB(1,2)-l3*PB(1,1))/(l5*PB(1,1)+d3tB*PB(1,1)+l3*PB(1,2))));qB=[theta1B d2tB d3tB PB(1,4)];
d2tC=PC(1,3);d3tC=sqrt(PC(1,1)^2+PC(1,2)^2-l3^2)-l5;theta1C=radtodeg(atan((l5*PC(1,2)+d3tC*PC(1,2)-l3*PC(1,1))/(l5*PC(1,1)+d3tC*PC(1,1)+l3*PC(1,2))));qC=[theta1C d2tC d3tC PC(1,4)];
d2tD=PD(1,3);d3tD=sqrt(PD(1,1)^2+PD(1,2)^2-l3^2)-l5;theta1D=radtodeg(atan((l5*PD(1,2)+d3tD*PD(1,2)-l3*PD(1,1))/(l5*PD(1,1)+d3tD*PD(1,1)+l3*PD(1,2))));qD=[theta1D d2tD d3tD PD(1,4)];
%%
%deytero erwthma
%oi xronoi tmax
wmax=1;vmax=0.5;
t1AB=(theta1B-theta1A)/wmax;t2AB=(d2tB-d2tA)/vmax;t3AB=(d3tB-d3tA)/vmax;t4AB=(PB(1,4)-PA(1,4))/wmax;tABmax=max([t1AB,t2AB,t3AB,t4AB]);
t1BC=(theta1C-theta1B)/wmax;t2BC=(d2tC-d2tB)/vmax;t3BC=(d3tC-d3tB)/vmax;t4BC=(PC(1,4)-PB(1,4))/wmax;tBCmax=max([t1BC,t2BC,t3BC,t4BC]);
t1CD=(theta1D-theta1C)/wmax;t2CD=(d2tD-d2tC)/vmax;t3CD=(d3tD-d3tC)/vmax;t4CD=(PD(1,4)-PC(1,4))/wmax;tCDmax=max([t1CD,t2CD,t3CD,t4CD]);T=tABmax+tBCmax+tCDmax;
%%
%trito erwthma


%i=1 theta 1
syms A B C D E F
%AB
eqns=[];eqns = [E == 0, D==0 , F==qA(1,1) , 5*A*(tABmax)^4+4*B*(tABmax)^3+3*C*(tABmax)^2+2*D*(tABmax)+E==(qC(1,1)-qB(1,1))/(tBCmax) ...
    ,20*A*(tABmax)^3+12*B*(tABmax)^2+6*C*(tABmax)+2*D==0 , A*(tABmax)^5+B*(tABmax)^4+C*(tABmax)^3+D*(tABmax)^2+E*(tABmax)+F==qB(1,1)];
vars=[];vars = [A B C D E F];[solA,solB,solC,solD,solE,solF] = solve(eqns, vars);poly=[];poly=[double(solA) double(solB) double(solC) double(solD) double(solE) double(solF)];
polyd=polyder(poly);polydd=polyder(polyd);tAB=0:0.01:tABmax;sqAB1=polyval(poly,tAB);
dsqAB1=polyval(polyd,tAB);ddsqAB1=polyval(polydd,tAB);
%BC
eqns=[];eqns=[A*tABmax+B==qB(1,1),A*(tABmax+tBCmax)+B==qC(1,1)];vars=[];vars = [A B];
[solA,solB] = solve(eqns, vars);poly=[];poly=[double(solA),double(solB)];
polyd=polyder(poly);polydd=polyder(polyd);tBC=tABmax:0.01:tABmax+tBCmax;sqBC1=polyval(poly,tBC);
dsqBC1=polyval(polyd,tBC);ddsqBC1=polyval(polydd,tBC);
%CD
eqns=[];eqns = [A*(tABmax+tBCmax+tCDmax)^5+B*(tABmax+tBCmax+tCDmax)^4+C*(tABmax+tBCmax+tCDmax)^3+D*(tABmax+tBCmax+tCDmax)^2+E*(tABmax+tBCmax+tCDmax)+F==qD(1,1)...
    5*A*(tABmax+tBCmax+tCDmax)^4+4*B*(tABmax+tBCmax+tCDmax)^3+3*C*(tABmax+tBCmax+tCDmax)^2+2*D*(tABmax+tBCmax+tCDmax)+E==0 ...
    20*A*(tABmax+tBCmax+tCDmax)^3+12*B*(tABmax+tBCmax+tCDmax)^2+6*C*(tABmax+tBCmax+tCDmax)+2*D==0 ...
    5*A*(tABmax+tBCmax)^4+4*B*(tABmax+tBCmax)^3+3*C*(tABmax+tBCmax)^2+2*D*(tABmax+tBCmax)+E==(qC(1,1)-qB(1,1))/(tBCmax) ...
    20*A*(tABmax+tBCmax)^3+12*B*(tABmax+tBCmax)^2+6*C*(tABmax+tBCmax)+2*D==0 ...
    ,A*(tABmax+tBCmax)^5+B*(tABmax+tBCmax)^4+C*(tABmax+tBCmax)^3+D*(tABmax+tBCmax)^2+E*(tABmax+tBCmax)+F==qC(1,1)];
vars=[];vars = [A B C D E F];[solA,solB,solC,solD,solE,solF] = solve(eqns, vars);poly=[];poly=[double(solA) double(solB) double(solC) double(solD) double(solE) double(solF)];
polyd=polyder(poly);polydd=polyder(polyd);tCD=tABmax+tBCmax:0.01:tABmax+tBCmax+tCDmax;sqCD1=polyval(poly,tCD);
dsqCD1=polyval(polyd,tCD);ddsqCD1=polyval(polydd,tCD);
%seq gia i=1
ti1=[tAB tBC tCD];sqi1=[sqAB1 sqBC1 sqCD1];dsqi1=[dsqAB1 dsqBC1 dsqCD1];
ddsqi1=[ddsqAB1 ddsqBC1 ddsqCD1];
figure;subplot(3,1,1);plot(ti1,sqi1);title('istoria \theta 1');ylabel('\theta 1 [degrees]');grid on;
subplot(3,1,2);plot(ti1,dsqi1);title('taxuthta \theta 1');ylabel('d{\theta 1}/dt [rad/sec]');grid on;
subplot(3,1,3);plot(ti1,ddsqi1);title('epitaxunsh \theta 1');ylabel('dd{\theta 1}/ddt [rad/sec2]');xlabel('t [sec]');grid on;



%i=2 d2
syms A B C D E F
%AB
eqns=[];eqns = [E == 0, D==0 , F==qA(1,2) , 5*A*(tABmax)^4+4*B*(tABmax)^3+3*C*(tABmax)^2+2*D*(tABmax)+E==(qC(1,2)-qB(1,2))/(tBCmax) ...
    ,20*A*(tABmax)^3+12*B*(tABmax)^2+6*C*(tABmax)+2*D==0 , A*(tABmax)^5+B*(tABmax)^4+C*(tABmax)^3+D*(tABmax)^2+E*(tABmax)+F==qB(1,2)];
vars=[];vars = [A B C D E F];[solA,solB,solC,solD,solE,solF] = solve(eqns, vars);poly=[];poly=[double(solA) double(solB) double(solC) double(solD) double(solE) double(solF)];
polyd=polyder(poly);polydd=polyder(polyd);tAB=0:0.01:tABmax;sqAB2=polyval(poly,tAB);
dsqAB2=polyval(polyd,tAB);ddsqAB2=polyval(polydd,tAB);
%BC
eqns=[];eqns=[A*tABmax+B==qB(1,2),A*(tABmax+tBCmax)+B==qC(1,2)];vars=[];vars = [A B];
[solA,solB] = solve(eqns, vars);poly=[];poly=[double(solA),double(solB)];
polyd=polyder(poly);polydd=polyder(polyd);tBC=tABmax:0.01:tABmax+tBCmax;sqBC2=polyval(poly,tBC);
dsqBC2=polyval(polyd,tBC);ddsqBC2=polyval(polydd,tBC);
%CD
eqns=[];eqns = [A*(tABmax+tBCmax+tCDmax)^5+B*(tABmax+tBCmax+tCDmax)^4+C*(tABmax+tBCmax+tCDmax)^3+D*(tABmax+tBCmax+tCDmax)^2+E*(tABmax+tBCmax+tCDmax)+F==qD(1,2)...
    5*A*(tABmax+tBCmax+tCDmax)^4+4*B*(tABmax+tBCmax+tCDmax)^3+3*C*(tABmax+tBCmax+tCDmax)^2+2*D*(tABmax+tBCmax+tCDmax)+E==0 ...
    20*A*(tABmax+tBCmax+tCDmax)^3+12*B*(tABmax+tBCmax+tCDmax)^2+6*C*(tABmax+tBCmax+tCDmax)+2*D==0 ...
    5*A*(tABmax+tBCmax)^4+4*B*(tABmax+tBCmax)^3+3*C*(tABmax+tBCmax)^2+2*D*(tABmax+tBCmax)+E==(qC(1,2)-qB(1,2))/(tBCmax) ...
    20*A*(tABmax+tBCmax)^3+12*B*(tABmax+tBCmax)^2+6*C*(tABmax+tBCmax)+2*D==0 ...
    ,A*(tABmax+tBCmax)^5+B*(tABmax+tBCmax)^4+C*(tABmax+tBCmax)^3+D*(tABmax+tBCmax)^2+E*(tABmax+tBCmax)+F==qC(1,2)];
vars=[];vars = [A B C D E F];[solA,solB,solC,solD,solE,solF] = solve(eqns, vars);poly=[];poly=[double(solA) double(solB) double(solC) double(solD) double(solE) double(solF)];
polyd=polyder(poly);polydd=polyder(polyd);tCD=tABmax+tBCmax:0.01:tABmax+tBCmax+tCDmax;sqCD2=polyval(poly,tCD);
dsqCD2=polyval(polyd,tCD);ddsqCD2=polyval(polydd,tCD);
%seq gia i=1
ti2=[tAB tBC tCD];sqi2=[sqAB2 sqBC2 sqCD2];dsqi2=[dsqAB2 dsqBC2 dsqCD2];
ddsqi2=[ddsqAB2 ddsqBC2 ddsqCD2];
figure;subplot(3,1,1);plot(ti2,sqi2);title('istoria d2');ylabel('d2 [degrees]');grid on;
subplot(3,1,2);plot(ti2,dsqi2);title('taxuthta d2');ylabel('d{d2}/dt [rad/sec]');grid on;
subplot(3,1,3);plot(ti2,ddsqi2);title('epitaxunsh d2');ylabel('dd{d2}/ddt [rad/sec2]');xlabel('t [sec]');grid on;


%i=3 d3
syms A B C D E F
%AB
eqns=[];eqns = [E == 0, D==0 , F==qA(1,3) , 5*A*(tABmax)^4+4*B*(tABmax)^3+3*C*(tABmax)^2+2*D*(tABmax)+E==(qC(1,3)-qB(1,3))/(tBCmax) ...
    ,20*A*(tABmax)^3+12*B*(tABmax)^2+6*C*(tABmax)+2*D==0 , A*(tABmax)^5+B*(tABmax)^4+C*(tABmax)^3+D*(tABmax)^2+E*(tABmax)+F==qB(1,3)];
vars=[];vars = [A B C D E F];[solA,solB,solC,solD,solE,solF] = solve(eqns, vars);poly=[];poly=[double(solA) double(solB) double(solC) double(solD) double(solE) double(solF)];
polyd=polyder(poly);polydd=polyder(polyd);tAB=0:0.01:tABmax;sqAB3=polyval(poly,tAB);
dsqAB3=polyval(polyd,tAB);ddsqAB3=polyval(polydd,tAB);
%BC
eqns=[];eqns=[A*tABmax+B==qB(1,3),A*(tABmax+tBCmax)+B==qC(1,3)];vars=[];vars = [A B];
[solA,solB] = solve(eqns, vars);poly=[];poly=[double(solA),double(solB)];
polyd=polyder(poly);polydd=polyder(polyd);tBC=tABmax:0.01:tABmax+tBCmax;sqBC3=polyval(poly,tBC);
dsqBC3=polyval(polyd,tBC);ddsqBC3=polyval(polydd,tBC);
%CD
eqns=[];eqns = [A*(tABmax+tBCmax+tCDmax)^5+B*(tABmax+tBCmax+tCDmax)^4+C*(tABmax+tBCmax+tCDmax)^3+D*(tABmax+tBCmax+tCDmax)^2+E*(tABmax+tBCmax+tCDmax)+F==qD(1,3)...
    5*A*(tABmax+tBCmax+tCDmax)^4+4*B*(tABmax+tBCmax+tCDmax)^3+3*C*(tABmax+tBCmax+tCDmax)^2+2*D*(tABmax+tBCmax+tCDmax)+E==0 ...
    20*A*(tABmax+tBCmax+tCDmax)^3+12*B*(tABmax+tBCmax+tCDmax)^2+6*C*(tABmax+tBCmax+tCDmax)+2*D==0 ...
    5*A*(tABmax+tBCmax)^4+4*B*(tABmax+tBCmax)^3+3*C*(tABmax+tBCmax)^2+2*D*(tABmax+tBCmax)+E==(qC(1,3)-qB(1,3))/(tBCmax) ...
    20*A*(tABmax+tBCmax)^3+12*B*(tABmax+tBCmax)^2+6*C*(tABmax+tBCmax)+2*D==0 ...
    ,A*(tABmax+tBCmax)^5+B*(tABmax+tBCmax)^4+C*(tABmax+tBCmax)^3+D*(tABmax+tBCmax)^2+E*(tABmax+tBCmax)+F==qC(1,3)];
vars=[];vars = [A B C D E F];[solA,solB,solC,solD,solE,solF] = solve(eqns, vars);poly=[];poly=[double(solA) double(solB) double(solC) double(solD) double(solE) double(solF)];
polyd=polyder(poly);polydd=polyder(polyd);tCD=tABmax+tBCmax:0.01:tABmax+tBCmax+tCDmax;sqCD3=polyval(poly,tCD);
dsqCD3=polyval(polyd,tCD);ddsqCD3=polyval(polydd,tCD);
%seq gia i=1
ti3=[tAB tBC tCD];sqi3=[sqAB3 sqBC3 sqCD3];dsqi3=[dsqAB3 dsqBC3 dsqCD3];
ddsqi3=[ddsqAB3 ddsqBC3 ddsqCD3];
figure;subplot(3,1,1);plot(ti3,sqi3);title('istoria d3');ylabel('d3 [degrees]');grid on;
subplot(3,1,2);plot(ti3,dsqi3);title('taxuthta d3');ylabel('d{d3}/dt [rad/sec]');grid on;
subplot(3,1,3);plot(ti3,ddsqi3);title('epitaxunsh d3');ylabel('dd{d3}/ddt [rad/sec2]');xlabel('t [sec]');grid on;


%i=4 theta 4
syms A B C D E F
%AB
eqns=[];eqns = [E == 0, D==0 , F==qA(1,4) , 5*A*(tABmax)^4+4*B*(tABmax)^3+3*C*(tABmax)^2+2*D*(tABmax)+E==(qC(1,4)-qB(1,4))/(tBCmax) ...
    ,20*A*(tABmax)^3+12*B*(tABmax)^2+6*C*(tABmax)+2*D==0 , A*(tABmax)^5+B*(tABmax)^4+C*(tABmax)^3+D*(tABmax)^2+E*(tABmax)+F==qB(1,4)];
vars=[];vars = [A B C D E F];[solA,solB,solC,solD,solE,solF] = solve(eqns, vars);poly=[];poly=[double(solA) double(solB) double(solC) double(solD) double(solE) double(solF)];
polyd=polyder(poly);polydd=polyder(polyd);tAB=0:0.01:tABmax;sqAB4=polyval(poly,tAB);
dsqAB4=polyval(polyd,tAB);ddsqAB4=polyval(polydd,tAB);
%BC
eqns=[];eqns=[A*tABmax+B==qB(1,4),A*(tABmax+tBCmax)+B==qC(1,4)];vars=[];vars = [A B];
[solA,solB] = solve(eqns, vars);poly=[];poly=[double(solA),double(solB)];
polyd=polyder(poly);polydd=polyder(polyd);tBC=tABmax:0.01:tABmax+tBCmax;sqBC4=polyval(poly,tBC);
dsqBC4=polyval(polyd,tBC);ddsqBC4=polyval(polydd,tBC);
%CD
eqns=[];eqns = [A*(tABmax+tBCmax+tCDmax)^5+B*(tABmax+tBCmax+tCDmax)^4+C*(tABmax+tBCmax+tCDmax)^3+D*(tABmax+tBCmax+tCDmax)^2+E*(tABmax+tBCmax+tCDmax)+F==qD(1,4)...
    5*A*(tABmax+tBCmax+tCDmax)^4+4*B*(tABmax+tBCmax+tCDmax)^3+3*C*(tABmax+tBCmax+tCDmax)^2+2*D*(tABmax+tBCmax+tCDmax)+E==0 ...
    20*A*(tABmax+tBCmax+tCDmax)^3+12*B*(tABmax+tBCmax+tCDmax)^2+6*C*(tABmax+tBCmax+tCDmax)+2*D==0 ...
    5*A*(tABmax+tBCmax)^4+4*B*(tABmax+tBCmax)^3+3*C*(tABmax+tBCmax)^2+2*D*(tABmax+tBCmax)+E==(qC(1,4)-qB(1,4))/(tBCmax) ...
    20*A*(tABmax+tBCmax)^3+12*B*(tABmax+tBCmax)^2+6*C*(tABmax+tBCmax)+2*D==0 ...
    ,A*(tABmax+tBCmax)^5+B*(tABmax+tBCmax)^4+C*(tABmax+tBCmax)^3+D*(tABmax+tBCmax)^2+E*(tABmax+tBCmax)+F==qC(1,4)];
vars=[];vars = [A B C D E F];[solA,solB,solC,solD,solE,solF] = solve(eqns, vars);poly=[];poly=[double(solA) double(solB) double(solC) double(solD) double(solE) double(solF)];
polyd=polyder(poly);polydd=polyder(polyd);tCD=tABmax+tBCmax:0.01:tABmax+tBCmax+tCDmax;sqCD4=polyval(poly,tCD);
dsqCD4=polyval(polyd,tCD);ddsqCD4=polyval(polydd,tCD);
%seq gia i=1
ti4=[tAB tBC tCD];sqi4=[sqAB4 sqBC4 sqCD4];dsqi4=[dsqAB4 dsqBC4 dsqCD4];
ddsqi4=[ddsqAB4 ddsqBC4 ddsqCD4];
figure;subplot(3,1,1);plot(ti4,radtodeg(sqi4));title('istoria \theta 4');ylabel('\theta 4 [degrees]');grid on;
subplot(3,1,2);plot(ti4,dsqi4);title('taxuthta \theta 4');ylabel('d{\theta 4}/dt [rad/sec]');grid on;
subplot(3,1,3);plot(ti4,ddsqi4);title('epitaxunsh \theta 4');ylabel('dd{\theta 4}/ddt [rad/sec2]');xlabel('t [sec]');grid on;

%%
% tetarto erwthma
d2max=max(sqi2);d3max=max(sqi3);l1=d2max+0.3;l4=d3max+0.4;i=1;tt=ti4;
thita1=sqi1;d2=sqi2;d3=sqi3;thita4=sqi4;
vthita1=dsqi1;vd2=dsqi2;vd3=dsqi3;vthita4=dsqi4;
athita1=ddsqi1;ad2=ddsqi2;ad3=ddsqi3;athita4=ddsqi4;

for t=tt
    A04=[cos(thita4(i))*sin(thita1(i)) sin(thita1(i))*sin(thita4(i)) cos(thita1(i)) l5*cos(thita1(i))+d3(i)*cos(thita1(i))-l3*sin(thita1(i)) ;
    cos(thita4(i))*cos(thita1(i))  -cos(thita1(i))*sin(thita4(i)) sin(thita1(i)) l5*sin(thita1(i))+d3(i)*sin(thita1(i))+l3*cos(thita1(i)) ;
    sin(thita4(i)) cos(thita4(i)) 0 d2(i) ; 0 0 0 1];
    rd4=[0 0 0 1]';
    rd04= A04*rd4;
    rd04x(i)= rd04(1,1); rd04y(i)=rd04(2,1); rd04z(i)=rd04(3,1);
    i=i+1;
end

T1=tt;figure
subplot(3,1,1);plot(T1,rd04x);title('diagramma theshs D , x-dir ');ylabel('rdx [m]');grid on;
subplot(3,1,2);plot(T1,rd04y);title('diagramma theshs D , y-dir ');ylabel('rdy [m]');grid on;
subplot(3,1,3);plot(T1,rd04z);title('diagramma theshs D , z-dir ');xlabel('time [sec]');ylabel('rdz [m]');grid on;


i=1;
for t=tt
    
    A01=[cos(thita1(i)) -sin(thita1(i)) 0 0 ; sin(thita1(i)) cos(thita1(i)) 0 0 ; 0 0 1 0 ; 0 0 0 1 ];
    A12=[0 0 1 0 ; 1 0 0 l3 ; 0 1 0 d2(i) ; 0 0 0 1];
    A23=[1 0 0 0 ; 0 1 0 0 ; 0 0 1 d3(i) ; 0 0 0 1];
    A34=[cos(thita4(i)) -sin(thita4(i)) 0 0 ; sin(thita4(i)) cos(thita4(i)) 0 0 ; 0 0 1 l5 ; 0 0 0 1]; 
    A04=[cos(thita4(i))*sin(thita1(i)) sin(thita1(i))*sin(thita4(i)) cos(thita1(i)) l5*cos(thita1(i))+d3(i)*cos(thita1(i))-l3*sin(thita1(i)) ;
    cos(thita4(i))*cos(thita1(i))  -cos(thita1(i))*sin(thita4(i)) sin(thita1(i)) l5*sin(thita1(i))+d3(i)*sin(thita1(i))+l3*cos(thita1(i)) ;
    sin(thita4(i)) cos(thita4(i)) 0 d2(i) ; 0 0 0 1];

    rd4=[0 0 0 1]';
    Q1=[0 -1 0 0 ; 1 0 0 0 ; 0 0 0 0 ; 0 0 0 0];
    Q2=[0 0 0 0 ; 0 0 0 0 ; 0 0 0 1 ; 0 0 0 0];
    
    U41=Q1*A04;
    U42=A01*Q2*A12*A23*A34;
    U43=A01*A12*Q2*A23*A34;
    U44=A01*A12*A23*Q1*A34;
    
    vd0=( U41*vthita1(i) + U42*vd2(i) + U43*vd3(i) + U44*vthita4(i) )*rd4;
    vdx(i)=vd0(1,1);vdy(i)=vd0(2,1);vdz(i)=vd0(3,1);
    
    i=i+1;
end

T2=tt;figure
subplot(3,1,1);plot(T2,rd04x);title('diagramma taxuthtas D , x-dir ');ylabel('vdx [m]');grid on;
subplot(3,1,2);plot(T2,rd04y);title('diagramma taxuthtas D , y-dir ');ylabel('vdy [m]');grid on;
subplot(3,1,3);plot(T2,rd04z);title('diagramma taxuthtas D , z-dir ');xlabel('time [sec]');ylabel('vdz [m]');grid on;



i=1;
for t=tt
    
    A01=[cos(thita1(i)) -sin(thita1(i)) 0 0 ; sin(thita1(i)) cos(thita1(i)) 0 0 ; 0 0 1 0 ; 0 0 0 1 ];
    A12=[0 0 1 0 ; 1 0 0 l3 ; 0 1 0 d2(i) ; 0 0 0 1];
    A23=[1 0 0 0 ; 0 1 0 0 ; 0 0 1 d3(i) ; 0 0 0 1];
    A34=[cos(thita4(i)) -sin(thita4(i)) 0 0 ; sin(thita4(i)) cos(thita4(i)) 0 0 ; 0 0 1 l5 ; 0 0 0 1]; 
    A04=[cos(thita4(i))*sin(thita1(i)) sin(thita1(i))*sin(thita4(i)) cos(thita1(i)) l5*cos(thita1(i))+d3(i)*cos(thita1(i))-l3*sin(thita1(i)) ;
    cos(thita4(i))*cos(thita1(i))  -cos(thita1(i))*sin(thita4(i)) sin(thita1(i)) l5*sin(thita1(i))+d3(i)*sin(thita1(i))+l3*cos(thita1(i)) ;
    sin(thita4(i)) cos(thita4(i)) 0 d2(i) ; 0 0 0 1];

    rd4=[0 0 0 1]';
    Q1=[0 -1 0 0 ; 1 0 0 0 ; 0 0 0 0 ; 0 0 0 0];
    Q2=[0 0 0 0 ; 0 0 0 0 ; 0 0 0 1 ; 0 0 0 0];
    Q3=[0 0 0 0 ; 0 0 0 0 ; 0 0 0 1 ; 0 0 0 0];
    Q4=[0 -1 0 0 ; 1 0 0 0 ; 0 0 0 0 ; 0 0 0 0];
    
    U411=Q1*Q1*A04;
    U412=Q1*A01*Q2*A12*A23*A34;
    U413=Q1*A01*A12*Q3*A23*A34;
    U414=Q1*A01*A12*A23*Q4*A34;
    U421=Q1*A01*Q2*A12*A23*A34;
    U422=A01*Q2*Q2*A12*A23*A34;
    U423=A01*Q2*A12*Q3*A23*A34;
    U424=A01*Q2*A12*A23*Q4*A34;
    U431=Q1*A01*A12*Q3*A23*A34;
    U432=A01*Q2*A12*Q3*A23*A34;
    U433=A01*A12*Q3*Q3*A23*A34;
    U434=A01*A12*Q3*A23*Q4*A34;
    U441=Q1*A01*A12*A23*Q4*A34;
    U442=A01*Q2*A12*A23*Q4*A34;
    U443=A01*A12*Q3*A23*Q4*A34;
    U444=A01*A12*A34*Q4*Q4*A34;
    
    U41=Q1*A04;
    U42=A01*Q2*A12*A23*A34;
    U43=A01*A12*Q2*A23*A34;
    U44=A01*A12*A23*Q1*A34;
    
    
    ad0=( U411*vthita1(i)*vthita1(i)+ U412*vthita1(i)*vd2(i)+ U413*vthita1(i)*vthita1(i) + U414*vthita1(i)*vthita4(i)+ ...
    U421*vd2(i)*vthita1(i)+ U422*vd2(i)*vd2(i) + U423*vd2(i)*vd3(i) + U424*vd2(i)*vthita4(i) + ...
    U431*vd3(i)*vthita1(i) + U432*vd3(i)*vd2(i) + U433*vd3(i)*vd3(i) + U434*vd3(i)*vthita4(i) + ...
    U441*vthita4(i)*vthita1(i) +  U442*vthita4(i)*vd2(i) +  U443*vthita4(i)*vd3(i) +  U444*vthita4(i)*vthita4(i)+ ...
    U41*athita1(i) + U42*ad2(i) + U43*ad3(i) + U44*athita4(i) )*rd4;
    
    adx(i)=ad0(1,1);ady(i)=ad0(2,1);adz(i)=ad0(3,1);
    
    i=i+1;
end

T3=tt;figure
subplot(3,1,1);plot(T3,rd04x);title('diagramma epitaxunshs D , x-dir ');ylabel('adx [m]');grid on;
subplot(3,1,2);plot(T3,rd04y);title('diagramma epitaxunshs D , y-dir ');ylabel('ady [m]');grid on;
subplot(3,1,3);plot(T3,rd04z);title('diagramma epitaxunshs D , z-dir ');xlabel('time [sec]');ylabel('adz [m]');grid on;


%%
%pempto erwthma
i=1;
for t=tt
    
w40=[abs(vthita4(i))*cos(thita1(i));abs(vthita4(i))*sin(thita1(i));abs(vthita1(i))];

w40x(i)=w40(1,1);
w40y(i)=w40(2,1);
w40z(i)=w40(3,1);

w(i)=sqrt(w40(1,1)^2+w40(2,1)^2+w40(3,1)^2);

e40=[abs(athita4(i))*cos(thita1(i))-abs(vthita1(i))*abs(vthita4(i))*sin(thita1(i)) 
   ;abs(athita4(i))*sin(thita1(i))+abs(vthita1(i))*abs(vthita4(i))*cos(thita1(i)) 
   ;abs(athita1(i))];

e40x(i)=e40(1,1);e40y(i)=e40(2,1);e40z(i)=e40(3,1);

e(i)=sqrt(e40(1,1)^2+e40(2,1)^2+e40(3,1)^2);

i=i+1;
end


T4=tt;figure
subplot(3,1,1);plot(T4,w40x);title('diagramma apoluths gwniakhs taxythtas , x-dir ');ylabel('{\omega}_{x} [rad/sec]');grid on;
subplot(3,1,2);plot(T4,w40y);title('diagramma apoluths gwniakhs taxythtas , y-dir ');ylabel('{\omega}_{y} [rad/sec]');grid on;
subplot(3,1,3);plot(T4,w40z);title('diagramma apoluths gwniakhs taxythtas , z-dir ');xlabel('time [sec]');ylabel('{\omega}_{z} [rad/sec]');grid on;

T4=tt;figure
subplot(3,1,1);plot(T4,e40x);title('diagramma apoluths gwniakhs epitaxunshs , x-dir ');ylabel('{\alpha}_{x} [read/sec^2]');grid on;
subplot(3,1,2);plot(T4,e40y);title('diagramma apoluths gwniakhs epitaxunshs , y-dir ');ylabel('{\alpha}_{y} [read/sec^2]');grid on;
subplot(3,1,3);plot(T4,e40z);title('diagramma apoluths gwniakhs epitaxunshs , z-dir ');xlabel('time [sec]');ylabel('{\alpha}_{z} [read/sec^2]');grid on;

T4=tt;figure
subplot(2,1,1);plot(T4,w);title('diagramma apoluths gwniakhs taxuthtas ');ylabel('{\omega} [rad/sec]');grid on;
subplot(2,1,2);plot(T4,e);title('diagramma apoluths gwniakhs epitaxunshs  ');xlabel('time [sec]');ylabel('{\alpha} [read/sec^2]');grid on;

%%
%II Meros
%prwto erwthma
m1=30;m2=25;m3=20;m4=4;m=1;
r1=l2/6;
h1=l2/2;
m1b=2.7;
m1a=m1-m1b;
l1=1.18-h1;
m3a=m3-m1b;
m4a=1.5;
m4pl=0.83;
t1=1.02;
t2=0.35;
t3=2.83;
xs2=-0.055;
ys2=-0.17;
zs2=0.0144;
zs4=(-3/4*l5*m4a-l5*m4pl)/m4;

jzz1b=m1b*r1^2/2
jyy1b=m1b*(3*r1^2+h1^2)
jxx1b=m1b*(3*r1^2+h1^2)

jxx1a=(1/3)*m1a*l1^2
jyy1a=(1/3)*m1a*l1^2
jzz1a=0;

j1=[jxx1a+jxx1b 0 0 ; 0 jyy1a+jyy1b 0 ; 0 0 jzz1a+jzz1b]

jxxa=0.071;
jyya=0.282;
jzza=0.212;
jyyb=0.47;
jxxb=0.73;
jzzb=1.2;
xsa=-0.095;
ysa=-0.28;
zsa=-0.094;
xsb=0.055;
ysb=0.17;
zsb=0.056;

ma=9.4;
mb=15.6;
m2=25;

jxx2=jxxa+ma*(ysa^2+zsa^2)+jxxb+mb*(ysb^2+zsb^2)
jyy2=jyya+ma*(zsa^2+xsa^2)+jyyb+mb*(zsb^2+xsb^2)
jzz2=jzza+ma*(xsa^2+ysa^2)+jzzb+mb*(xsb^2+ysb^2)
jxy2=ma*xsa*ysa+mb*xsb*ysb
jxz2=ma*xsa*zsa+mb*xsb*zsb
jyz2=ma*ysa*zsa+mb*ysb*zsb
j2=[jxx2 -jxy2 -jxz2 ; -jxy2 jyy2 -jyz2 ; -jxz2 -jyz2 jzz2]



jzz3b=m1b*r1^2/2
jyy3b=m1b*(3*r1^2+h1^2)
jxx3b=m1b*(3*r1^2+h1^2)

jxx3a=(1/3)*m3a*l4^2
jyy3a=(1/3)*m3a*l4^2
jzz3a=0;

j3=[jxx3a+jxx3b 0 0 ; 0 jyy3a+jyy3b 0 ; 0 0 jzz3a+jzz3b]

jzz4a=m4a*(l2/10)^2/2
jyy4a=m4a*(3*(l2/10)^2+(l5/2)^2)/12
jxx4a=m4a*(3*(l2/10)^2+(l5/2)^2)/12

ixx1=(m4pl/12)*((l5/8)^2+(l5/16)^2)+m4pl*(l5/4)^2
iyy1=(m4pl/12)*((l5/2)^2+(l5/16)^2)+m4pl*(l5/4)^2
izz1=(m4pl/12)*((l5/2)^2+(l5/8)^2)+m4pl*(l5/4)^2

ixx2=(m4pl/12)*((l5/8)^2+(l5/16)^2)+m4pl*(l5/4)^2
iyy2=(m4pl/12)*((l5/2)^2+(l5/16)^2)+m4pl*(l5/4)^2
izz2=(m4pl/12)*((l5/2)^2+(l5/8)^2)+m4pl*(l5/4)^2

izz3=(m4pl/12)*((l5/2)^2+(l5/8)^2)
ixx3=(m4pl/12)*((l5/8)^2+(l5/16)^2)+m4pl*(l5/2)^2
iyy3=(m4pl/12)*((l5/2)^2+(l5/16)^2)+m4pl*(l5/2)^2

Ixx=izz1+izz2+ixx3
Iyy=iyy1+iyy2+iyy3
Izz=ixx1+ixx2+izz3

j4=[Ixx+jxx4a 0 0 ; 0 Iyy+jyy4a 0 ; 0 0 Izz+jzz4a]


%%
%deutero erwthma

i=1;

for t=tt
    
    vg10=0;
    vg20=[-vthita1(i)*ys2 vthita1(i)*xs2 abs(vd2(i))]';
    vg30=[-vthita1(i)*ys2+abs(vd3(i))*cos(thita1(i)) vthita1(i)*xs2+abs(vd3(i))*sin(thita1(i)) abs(vd2(i))]';
    vg40=[-vthita1(i)*ys2+abs(vd3(i))*cos(thita1(i))+vthita4(i)*sin(thita1(i))*zs4;...
        vthita1(i)*xs2+abs(vd3(i))*sin(thita1(i))-vthita4(i)*cos(thita1(i))*zs4;abs(vd2(i))];
    ag10=[0;0;0];
    ag20=[-athita1(i)*ys2-vthita1(i)^2*xs2;athita1(i)*xs2-vthita1(i)^2*ys2;abs(ad2(i))];
    ag30=ag20+[cos(thita1(i))*abs(ad3(i))-2*vthita1(i)*abs(vd3(i))*sin(thita1(i));...
       sin(thita1(i))*abs(ad3(i))+2*vthita1(i)*abs(vd3(i))*cos(thita1(i));0];
    ag40=ag30+[abs(athita4(i))*sin(thita1(i))*zs4+vthita1(i)*abs(vthita4(i))*cos(thita1(i))*zs4;...
        -abs(athita4(i))*cos(thita1(i))*zs4+vthita1(i)*abs(vthita4(i))*sin(thita1(i))*zs4;-zs4*(vthita4(i))^2];
    
        F1=-m1*ag10;
        F1X(i)=F1(1,1);
        F1Y(i)=F1(2,1);
        F1Z(i)=F1(3,1);
        
        F2=-m2*ag20;
        F2X(i)=F2(1,1);
        F2Y(i)=F2(2,1);
        F2Z(i)=F2(3,1);
        
        F3=-m3*ag30;
        F3X(i)=F3(1,1);
        F3Y(i)=F3(2,1);
        F3Z(i)=F3(3,1);
        
        F4=-m4*ag40;
        F4X(i)=F4(1,1);
        F4Y(i)=F4(2,1);
        F4Z(i)=F4(3,1);
        
        w10=[0;0;vthita1(i)];w20=w10;w30=w10;w40=[abs(vthita4(i))*cos(thita1(i));abs(vthita4(i))*sin(thita1(i));vthita1(i)];
        e10=[0;0;athita1(i)];e20=e10;e30=e10;e40=[abs(athita4(i))*cos(thita1(i));abs(athita4(i))*sin(thita1(i));athita1(i)];
        
        M1=-j1*e10-cross(w10,j1*w10);
        M1X(i)=M1(1,1);
        M1Y(i)=M1(2,1);
        M1Z(i)=M1(3,1);
        
        M2=-j2*e20-cross(w20,j2*w20);
        M2X(i)=M2(1,1);
        M2Y(i)=M2(2,1);
        M2Z(i)=M2(3,1);
        
        M3=-j3*e30-cross(w30,j3*w30);
        M3X(i)=M3(1,1);
        M3Y(i)=M3(2,1);
        M3Z(i)=M3(3,1);
        
        M4=-j4*e40-cross(w40,j4*w40);
        M4X(i)=M4(1,1);
        M4Y(i)=M4(2,1);
        M4Z(i)=M4(3,1);
    i=i+1;
end

T10=tt;figure
subplot(3,1,1);plot(T10,F1X);title('diagramma adraneiakhs dunamhs melous 1 - xdir');xlabel('time [sec]');ylabel('F1-x [Nt]');grid on
subplot(3,1,2);plot(T10,F1Y);title('diagramma adraneiakhs dunamhs melous 1 - ydir');xlabel('time [sec]');ylabel('F1-y [Nt]');grid on
subplot(3,1,3);plot(T10,F1Z);title('diagramma adraneiakhs dunamhs melous 1 - zdir');xlabel('time [sec]');ylabel('F1-z [Nt]');grid on

figure
subplot(3,1,1);plot(T10,F2X);title('diagramma adraneiakhs dunamhs melous 2 - xdir');xlabel('time [sec]');ylabel('F2-x [Nt]');grid on
subplot(3,1,2);plot(T10,F2Y);title('diagramma adraneiakhs dunamhs melous 2 - ydir');xlabel('time [sec]');ylabel('F2-y [Nt]');grid on
subplot(3,1,3);plot(T10,F2Z);title('diagramma adraneiakhs dunamhs melous 2 - zdir');xlabel('time [sec]');ylabel('F2-z [Nt]');grid on

figure
subplot(3,1,1);plot(T10,F3X);title('diagramma adraneiakhs dunamhs melous 3 - xdir');xlabel('time [sec]');ylabel('F3-x [Nt]');grid on
subplot(3,1,2);plot(T10,F3Y);title('diagramma adraneiakhs dunamhs melous 3 - ydir');xlabel('time [sec]');ylabel('F3-y [Nt]');grid on
subplot(3,1,3);plot(T10,F3Z);title('diagramma adraneiakhs dunamhs melous 3 - zdir');xlabel('time [sec]');ylabel('F3-z [Nt]');grid on

figure
subplot(3,1,1);plot(T10,F4X);title('diagramma adraneiakhs dunamhs melous 4 - xdir');xlabel('time [sec]');ylabel('F4-x [Nt]');grid on
subplot(3,1,2);plot(T10,F4Y);title('diagramma adraneiakhs dunamhs melous 4 - ydir');xlabel('time [sec]');ylabel('F4-y [Nt]');grid on
subplot(3,1,3);plot(T10,F4Z);title('diagramma adraneiakhs dunamhs melous 4 - zdir');xlabel('time [sec]');ylabel('F4-z [Nt]');grid on


figure
subplot(3,1,1);plot(T10,M1X);title('diagramma adraneiakhs rophs melous 1 - xdir');xlabel('time [sec]');ylabel('M1-x [Nt]');grid on
subplot(3,1,2);plot(T10,M1Y);title('diagramma adraneiakhs rophs melous 1 - ydir');xlabel('time [sec]');ylabel('M1-y [Nt]');grid on
subplot(3,1,3);plot(T10,M1Z);title('diagramma adraneiakhs rophs melous 1 - zdir');xlabel('time [sec]');ylabel('M1-z [Nt]');grid on

figure
subplot(3,1,1);plot(T10,M2X);title('diagramma adraneiakhs rophs melous 2 - xdir');xlabel('time [sec]');ylabel('M2-x [Nt]');grid on
subplot(3,1,2);plot(T10,M2Y);title('diagramma adraneiakhs rophs melous 2 - ydir');xlabel('time [sec]');ylabel('M2-y [Nt]');grid on
subplot(3,1,3);plot(T10,M2Z);title('diagramma adraneiakhs rophs melous 2 - zdir');xlabel('time [sec]');ylabel('M2-z [Nt]');grid on

figure
subplot(3,1,1);plot(T10,M3X);title('diagramma adraneiakhs rophs melous 3 - xdir');xlabel('time [sec]');ylabel('M3-x [Nt]');grid on
subplot(3,1,2);plot(T10,M3Y);title('diagramma adraneiakhs rophs melous 3 - ydir');xlabel('time [sec]');ylabel('M3-y [Nt]');grid on
subplot(3,1,3);plot(T10,M3Z);title('diagramma adraneiakhs rophs melous 3 - zdir');xlabel('time [sec]');ylabel('M3-z [Nt]');grid on

figure
subplot(3,1,1);plot(T10,M4X);title('diagramma adraneiakhs rophs melous 4 - xdir');xlabel('time [sec]');ylabel('M4-x [Nt]');grid on
subplot(3,1,2);plot(T10,M4Y);title('diagramma adraneiakhs rophs melous 4 - ydir');xlabel('time [sec]');ylabel('M4-y [Nt]');grid on
subplot(3,1,3);plot(T10,M4Z);title('diagramma adraneiakhs rophs melous 4 - zdir');xlabel('time [sec]');ylabel('M4-z [Nt]');grid on



%%
%trito erwthma


syms f01x f01y  f01z f12x f12y f12z f23x f23y f23z f34x f34y f34z M01x...
            M01y M01z M12x M12y M12z M23x M23y M23z M34x M34y M34z

g=9.81;
i=1;
for t=tt
    R01=[cos(thita1(i)) -sin(thita1(i)) 0  ; sin(thita1(i)) cos(thita1(i)) 0  ; 0 0 1 ];
    R12=[0 0 1  ; 1 0 0  ; 0 1 0 ];
    R23=[1 0 0  ; 0 1 0  ; 0 0 1 ];
    R34=[cos(thita4(i)) -sin(thita4(i)) 0 ;sin(thita4(i)) cos(thita4(i)) 0;0 0 1];
    r1=[0 0 0]';
    r1p=[0 0 l1/2]';
    r2=[l3 d2(i)+1.5*l3 0]';
    r2p=[xs2 ys2 zs2]';
    r3=[0 0 d3(i)]';
    r3p=[0 0 -l4/2]';
    r4=[0 0 l5]';
    r4p=[0 0 zs4]';
    g1=[0 0 -g]';
    g2=[0 -g 0]';
    g3=[0 -g 0]';
    g4=[g 0 0]';
    
   
    vg10=0;
    vg20=[-vthita1(i)*ys2 vthita1(i)*xs2 abs(vd2(i))]';
    vg30=[-vthita1(i)*ys2+abs(vd3(i))*cos(thita1(i)) vthita1(i)*xs2+abs(vd3(i))*sin(thita1(i)) abs(vd2(i))]';
    vg40=[-vthita1(i)*ys2+abs(vd3(i))*cos(thita1(i))+vthita4(i)*sin(thita1(i))*zs4;...
        vthita1(i)*xs2+abs(vd3(i))*sin(thita1(i))-vthita4(i)*cos(thita1(i))*zs4;abs(vd2(i))];
    ag10=[0;0;0];
    ag20=[-athita1(i)*ys2-vthita1(i)^2*xs2;athita1(i)*xs2-vthita1(i)^2*ys2;abs(ad2(i))];
    ag30=ag20+[cos(thita1(i))*abs(ad3(i))-2*vthita1(i)*abs(vd3(i))*sin(thita1(i));...
       sin(thita1(i))*abs(ad3(i))+2*vthita1(i)*abs(vd3(i))*cos(thita1(i));0];
    ag40=ag30+[abs(athita4(i))*sin(thita1(i))*zs4+vthita1(i)*abs(vthita4(i))*cos(thita1(i))*zs4;...
        -abs(athita4(i))*cos(thita1(i))*zs4+vthita1(i)*abs(vthita4(i))*sin(thita1(i))*zs4;-zs4*(vthita4(i))^2];
 
    
    cr1=cross((r1+r1p),ag10);
    cr2=cross((r2+r2p),ag20);
    cr3=cross((r3+r3p),ag30);
    cr4=cross((r4+r4p),ag40);
    cm1=cross((r1+r1p),g1);
    cm2=cross((r2+r2p),g2);
    cm3=cross((r3+r3p),g3);
    cm4=cross((r4+r4p),g4);
   
    
    eqns=[-f01x+m1*ag10(1,1)+R12(1,1)*f12x+R12(1,2)*f12y+R12(1,3)*f12z==0,...
    -f01y+m1*ag10(2,1)+R12(2,1)*f12x+R12(2,2)*f12y+R12(2,3)*f12z==0,...
    -f01z+m1*ag10(3,1)+R12(3,1)*f12x+R12(3,2)*f12y+R12(3,3)*f12z-m1*g==0,...
    -M01x+M1X(i)+m1*cr1(1,1)+R12(1,1)*M12x+R12(1,2)*M12y+R12(1,3)*M12z-m1*cm1(1,1)==0,...
    -M01y+M1Y(i)+m1*cr1(2,1)+R12(2,1)*M12x+R12(2,2)*M12y+R12(2,3)*M12z-m1*cm1(2,1)==0,...
    -M01z+M1Z(i)+m1*cr1(3,1)+R12(3,1)*M12x+R12(3,2)*M12y+R12(3,3)*M12z-m1*cm1(3,1)==0,...
    -f12x+m2*ag20(1,1)+R23(1,1)*f23x+R23(1,2)*f23y+R23(1,3)*f23z==0,...
    -f12y+m2*ag20(2,1)+R23(2,1)*f23x+R23(2,2)*f23y+R23(2,3)*f23z-m2*g==0,...
    -f12z+m2*ag20(3,1)+R23(3,1)*f23x+R23(3,2)*f23y+R23(3,3)*f23z==0,...
    -M12x+M2X(i)+m2*cr2(1,1)+R23(1,1)*M23x+R23(1,2)*M23y+R23(1,3)*M23z-m1*cm2(1,1)+r2(2,1)*f23z==0,...
    -M12y+M2Y(i)+m2*cr2(2,1)+R23(2,1)*M23x+R23(2,2)*M23y+R23(2,3)*M23z-m1*cm2(2,1)-r2(1,1)*f23z==0,...
    -M12z+M2Z(i)+m2*cr2(3,1)+R23(3,1)*M23x+R23(3,2)*M23y+R23(3,3)*M23z-m1*cm2(3,1)+r2(1,1)*f23y-f23x*r2(2,1)==0,...  
    -f23x+m3*ag30(1,1)+R34(1,1)*f34x+R34(1,2)*f34y+R34(1,3)*f34z==0,...
    -f23y+m3*ag30(2,1)+R34(2,1)*f34x+R34(2,2)*f34y+R34(2,3)*f34z-m3*g==0,...
    -f23z+m3*ag30(3,1)+R34(3,1)*f34x+R34(3,2)*f34y+R34(3,3)*f34z==0,...
    -M23x+M3X(i)+m3*cr3(1,1)+R34(1,1)*M34x+R34(1,2)*M34y+R34(1,3)*M34z-m3*cm3(1,1)+d3(i)*(f34y*cos(thita4(i))-f34x*sin(thita4(i)))==0,...
    -M23y+M3Y(i)+m3*cr3(2,1)+R34(2,1)*M34x+R34(2,2)*M34y+R34(2,3)*M34z-m3*cm3(2,1)+d3(i)*(f34x*cos(thita4(i))-f34y*sin(thita4(i)))==0,...
    -M23z+M3Z(i)+m3*cr3(3,1)+R34(3,1)*M34x+R34(3,2)*M34y+R34(3,3)*M34z-m3*cm3(3,1)==0,...  
    -f34x+(m4+2)*ag40(1,1)-(m4+2)*g==0,...
    -f34y+(m4+2)*ag40(2,1)==0,...
    -f34z+(m4+2)*ag40(3,1)==0,...
    -M34x+M4X(i)+m4*cr4(1,1)-m4*cm4(1,1)==0,...
    -M34y+M4Y(i)+m4*cr4(2,1)-m4*cm4(2,1)==0,...
    -M34z+M4Z(i)+m4*cr4(3,1)-m4*cm4(3,1)==0];
  
     S=solve(eqns);
     sol=[S.f01x; S.f01y ;S.f01z ;S.f12x; S.f12y; S.f12z; S.f23x ;S.f23y; S.f23z;S.f34x; S.f34y ;S.f34z;...
            S.M01x; S.M01y ;S.M01z ;S.M12x; S.M12y; S.M12z; S.M23x ;S.M23y; S.M23z;S.M34x; S.M34y ;S.M34z];
        [A,b] = equationsToMatrix(eqns,[f01x,f01y, f01z, f12x, f12y, f12z, f23x, f23y, f23z, f34x, f34y, f34z, M01x,...
            M01y, M01z, M12x, M12y, M12z, M23x ,M23y ,M23z ,M34x, M34y, M34z]);
     z = A\b;
     force(i,:)=double(z);
    
     i=i+1;
end
    
T11=tt;figure
subplot(3,1,1);plot(T11,force(:,1));title('diagramma apaitoumenhs kinhthrias dunamhs arhtrwshs 1 - xdir');xlabel('time [sec]');ylabel('F1-x [Nt]');grid on
subplot(3,1,2);plot(T11,force(:,2));title('diagramma apaitoumenhs kinhthrias dunamhs arhtrwshs 1 - ydir');xlabel('time [sec]');ylabel('F1-y [Nt]');grid on
subplot(3,1,3);plot(T11,force(:,3));title('diagramma apaitoumenhs kinhthrias dunamhs arhtrwshs 1 - zdir');xlabel('time [sec]');ylabel('F1-z [Nt]');grid on

figure
subplot(3,1,1);plot(T11,force(:,4));title('diagramma apaitoumenhs kinhthrias dunamhs arhtrwshs 2 - xdir');xlabel('time [sec]');ylabel('F2-x [Nt]');grid on
subplot(3,1,2);plot(T11,force(:,5));title('diagramma apaitoumenhs kinhthrias dunamhs arhtrwshs 2 - ydir');xlabel('time [sec]');ylabel('F2-y [Nt]');grid on
subplot(3,1,3);plot(T11,force(:,6));title('diagramma apaitoumenhs kinhthrias dunamhs arhtrwshs 2 - zdir');xlabel('time [sec]');ylabel('F2-z [Nt]');grid on

figure
subplot(3,1,1);plot(T11,force(:,7));title('diagramma apaitoumenhs kinhthrias dunamhs arhtrwshs 3 - xdir');xlabel('time [sec]');ylabel('F3-x [Nt]');grid on
subplot(3,1,2);plot(T11,force(:,8));title('diagramma apaitoumenhs kinhthrias dunamhs arhtrwshs 3 - ydir');xlabel('time [sec]');ylabel('F3-y [Nt]');grid on
subplot(3,1,3);plot(T11,force(:,9));title('diagramma apaitoumenhs kinhthrias dunamhs arhtrwshs 3 - zdir');xlabel('time [sec]');ylabel('F3-z [Nt]');grid on

figure
subplot(3,1,1);plot(T11,force(:,10));title('diagramma apaitoumenhs kinhthrias dunamhs arhtrwshs 4 - xdir');xlabel('time [sec]');ylabel('F4-x [Nt]');grid on
subplot(3,1,2);plot(T11,force(:,11));title('diagramma apaitoumenhs kinhthrias dunamhs arhtrwshs 4 - ydir');xlabel('time [sec]');ylabel('F4-y [Nt]');grid on
subplot(3,1,3);plot(T11,force(:,12));title('diagramma apaitoumenhs kinhthrias dunamhs arhtrwshs 4 - zdir');xlabel('time [sec]');ylabel('F4-z [Nt]');grid on


figure
subplot(3,1,1);plot(T11,force(:,13));title('diagramma apaitoumenhs kinhthrias rophs arhtrwshs 1 - xdir');xlabel('time [sec]');ylabel('M1-x [Nt]');grid on
subplot(3,1,2);plot(T11,force(:,14));title('diagramma apaitoumenhs kinhthrias rophs arhtrwshs 1 - ydir');xlabel('time [sec]');ylabel('M1-y [Nt]');grid on
subplot(3,1,3);plot(T11,force(:,15));title('diagramma apaitoumenhs kinhthrias rophs arhtrwshs 1 - zdir');xlabel('time [sec]');ylabel('M1-z [Nt]');grid on

figure
subplot(3,1,1);plot(T11,force(:,16));title('diagramma apaitoumenhs kinhthrias rophs arhtrwshs 2 - xdir');xlabel('time [sec]');ylabel('M2-x [Nt]');grid on
subplot(3,1,2);plot(T11,force(:,17));title('diagramma apaitoumenhs kinhthrias rophs arhtrwshs 2 - ydir');xlabel('time [sec]');ylabel('M2-y [Nt]');grid on
subplot(3,1,3);plot(T11,force(:,18));title('diagramma apaitoumenhs kinhthrias rophs arhtrwshs 2 - zdir');xlabel('time [sec]');ylabel('M2-z [Nt]');grid on

figure
subplot(3,1,1);plot(T11,force(:,19));title('diagramma apaitoumenhs kinhthrias rophs arhtrwshs 3 - xdir');xlabel('time [sec]');ylabel('M3-x [Nt]');grid on
subplot(3,1,2);plot(T11,force(:,20));title('diagramma apaitoumenhs kinhthrias rophs arhtrwshs 3 - ydir');xlabel('time [sec]');ylabel('M3-y [Nt]');grid on
subplot(3,1,3);plot(T11,force(:,21));title('diagramma apaitoumenhs kinhthrias rophs arhtrwshs 3 - zdir');xlabel('time [sec]');ylabel('M3-z [Nt]');grid on

figure
subplot(3,1,1);plot(T11,force(:,22));title('diagramma apaitoumenhs kinhthrias rophs arhtrwshs 4 - xdir');xlabel('time [sec]');ylabel('M4-x [Nt]');grid on
subplot(3,1,2);plot(T11,force(:,23));title('diagramma apaitoumenhs kinhthrias rophs arhtrwshs 4 - ydir');xlabel('time [sec]');ylabel('M4-y [Nt]');grid on
subplot(3,1,3);plot(T11,force(:,24));title('diagramma apaitoumenhs kinhthrias rophs arhtrwshs 4 - zdir');xlabel('time [sec]');ylabel('M4-z [Nt]');grid on

