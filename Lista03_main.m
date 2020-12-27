pkg load control
pkg load signal

% Lista Ejercicios 03 - Mamani Huanca Jhoel Rene
% SCA GrupoA

clear;clc;
sys = tf([0.00527],[1.211 1])
G = [sys 0     0;
       0  sys  0;
       0  0     sys];

[A,B,C,D] = tf2ss(G); % Del sistema

%$ %%%%%%%%%%%%%%%%%%%%%%%%% Ejercicio 01 %%%%%%%%%%%%%%%%%%%%%%%%%
%barreras de desempeño y estabilidad
%1ra barrera Baja frecuencia = 20dB, Wd = 0.5 rad/s
%2da barrera Alta frecuencia = -20dB, Wn = 100 rad/s
ad=0.1; wd=0.5;
ar=0.3; wr=0.8; 
an=0.1; wn=100;

w=logspace(-1,3,1000); %rango
g1=20*log10(1/ad)*ones(1,177);% 177 = 0.5 rad/s 
g2=zeros(1,573);
g3=20*log10(an)*ones(1,250);
gt=[g1 g2 g3];
gt2=[20*log10(1/ar)*ones(1,227) zeros(1,773)];
% barreras de estabilidad
semilogx(w, gt,'m'); 
print -djpg image1.jpg %para guardar la imagen (formato.jpg)

ylim([-40,40])
grid; hold on;
semilogx(w,gt2,'b'); title('Barreras del Sistema');xlabel('w(rad/s)');ylabel('dB');hold off;
print -djpg image2.jpg %para guardar la imagen (formato.jpg)


%%%%%%%%%%%%%%%%%%%%%%%%%%% EJERCICIO 02 %%%%%%%%%%%%%%%%%%%%%%%%%
%Valores singulares maximos y minimos en w =  1rad/s
w=logspace(-2,2,100);
sv=sigma(G,w);
Valor=sv(:,51) % 1 rad/s
sv=20*log10(sv);
Valor_db=sv(:,51) % 1 rad/s
figure; semilogx(w,sv); grid;
print -djpg image3.jpg %para guardar la imagen (formato.jpg)

[a,b,c,d] = tf2ss(G)

%%%%%%%%%%%%%%%%%%%%%%%%%%% EJERCICIO 03 %%%%%%%%%%%%%%%%%%%%%%%%%
%Calculo de la matriz de ganancia G con LQR.
Q = C'*C;
R = 0.1*eye(3);
GG = lqr(G,Q,R)

%%%%%%%%%%%%%%%%%%%%%%%%%%% EJERCICIO 04 %%%%%%%%%%%%%%%%%%%%%%%%%
%Calculo de la matriz de ganancia H con filtro de Kalman.
Q1=eye(3,3);
R1 = 0.01*eye(3);
HH = lqe(G,Q1,R1)

%%%%%%%%%%%%%%%%%%%%%%%%%% EJERCICIO 5 %%%%%%%%%%%%%%%%%%%%%%%%% 
%Escriba la función de transferencia del compensador K(s) 
%incluye a G y H.
K11=[A-B*GG]; 
K12=[-B*GG]; 
K21=[HH*C]; 
K22=[A-HH*C-B*GG];  

AK=[K11 K12;K21 K22];  
BK=[zeros(3,3);-HH];     
CK=[C zeros(3,3)];     
DK=[0 0 0;0 0 0;0 0 0];
% % % % % % % % % % % % % % % 
sysK=ss(AK,BK,CK,DK)

%%%%%%%%%%%%%%%%%%%%%%%%%% EJERCICIO 6 %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% sensibilidad
al=[zeros(3,3) zeros(3,3); B A];
bl=[eye(3,3) ;zeros(3,3)];
cl=[zeros(3,3) C];
dl=eye(3,3);

w=logspace(-2,3,100);
%%% sensibilidad
sv = sigma(ss(al-bl*cl, bl, -cl,dl),w);
sv = 20*log10(sv);


figure
semilogx(w, sv);title('Valores singulares de sensibilidad')
grid; xlabel('Frecuencia (rad/seg)'); ylabel('Valores Singulares (dB)')
print -djpg image4.jpg %para guardar la imagen (formato.jpg)

sv = sigma(ss(al-bl*cl, bl, -cl, 0*dl),w);
sv = 20*log10(sv);


figure
semilogx(w, sv)
title('Valores singulares de sensibilidad complementaria')
grid; xlabel('Frecuencia (rad/seg)'); ylabel('Valores Singulares (dB)')
print -djpg image5.jpg %para guardar la imagen (formato.jpg)
