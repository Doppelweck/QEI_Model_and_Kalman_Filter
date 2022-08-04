%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Digitale Signalverarbeitung                                %
%     Reducing of Encoder Noise by using Kalman Filter                    %
%       Sebastian Friedrich, Christina Liesenfeld                         %
%                       29.8.2018                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% delete memory ----------------------------------------------------------%
clc;
close all;
clear all;

% System times -----------------------------------------------------------%
T_QEI=0.001;                                 %Sample Time QEI
T = 0.001;                                   %Maximum Step Size Simulink Solver

% Motor Specifications----------------------------------------------------%
k_M=1/(466);                                 % Vs
J= 33.3 *10^-7+ ((1/121) *0.5* 28.65*10^-6); % kg m^2
RM=0.21;                                     % Ohm

% Discrete State Space System --------------------------------------------%
A= [1 T_QEI (T_QEI^2)/2; 0 1 T_QEI; 0 0 1];
B= [0 0 0]';                                 % No input signal available 
C= [1 0 0; 0 1 0];
I= eye(3);

% Error -----------------------------------------------------------------%
deltaP=35.714285*10^-6;                      % m per flank
varP=((deltaP+1)^2-1)/12;                    % Variance position
deltaV=deltaP/T_QEI;                         % m per s
varV=((deltaV+1)^2-1)/12;                    % Variance velocity

deltaAC=19.387*10^(-1);                      % m per s^2
Q=(deltaAC/3)^2;                             % Variance of system noise
R=[varP 0; 0 varV];                          % Measurement noise variance
G=[(T_QEI^2)/2 ; T_QEI ; 1];                 % Stochastic control matrix


%% transfer to simulink
mysim=simset ('MaxStep',T,'InitialState',[0 0 0]);
[tb,pos,y]=sim('QEI_Model_2016a');

%%
mesP=squeeze(QEI_Pos.signals.values);
mesV=squeeze(QEI_Vel.signals.values);
mesT=QEI_Pos.time;
mesIsValid=squeeze(QEI_NewValue.signals.values);
PosReal=squeeze(PosReal.signals.values);

%% estimated initial values
close all
xprp(:,1)=[0; 0; 0];
Pp = [1 0 0; 0 0.5 0; 0 0 0.1];
xT(1)=0;

UseInvalidMeasurments = 0;                      % True=1 false=0

% Kalman filter algorithm-------------------------------------------------%
for i = 1:((length(mesP)-1))
    
    if (mesIsValid(i)==0 && UseInvalidMeasurments)
        C= [0 0 0; 0 0 0];
    else
        C= [1 0 0; 0 1 0];
    end
    
    xmess(:,i+1)=[mesP((i)+1) ;mesV((i)+1); 0]'; % Measurment Value
    
    %Time-Update:----------------------------------------------------------
    xprm(:,i+1)=A*xprp(:,i);
    Pm=A*Pp*A'+G*Q*G';
    
    %Measurement-Update:---------------------------------------------------
    K=Pm*(C')*inv((C*Pm*(C'))+R);
    xprp(:,i+1)=xprm(:,i+1)+K*((C*xmess(:,i+1))-C*xprm(:,i+1));
    Pp=(I-K*C)*Pm;
    xT(i+1)=1*10^-2;
    
    K=Pm*(C')*inv((C*Pm*(C'))+R);
    xprp(:,i+1)=xprm(:,i+1)+K*((C*xmess(:,i+1))-C*xprm(:,i+1));
    Pp=(I-K*C)*Pm;
    
    Pm1(i,:)=Pm(1,:);
    Pm2(i,:)=Pm(2,:);
    Pm3(i,:)=Pm(3,:);
    
    K1(i,:)=K(1,:);
    K2(i,:)=K(2,:);
    K3(i,:)=K(3,:);
    
    Pp1(i,:)=Pp(1,:);
    Pp2(i,:)=Pp(2,:);
    Pp3(i,:)=Pp(3,:);
    
    if (squeeze(QEI_NewValue.signals.values(i+1)))==0
        xT(i+1)=0;
    else
        xT(i+1)=1*-0.01;
    end;
end;

% Mean square error ------------------------------------------------------%

schiff.orig = squeeze(U.signals.values);
schiff.mess = mesV;
schiff.x = xprp(2,:)';

schiff.std_orig = mean((schiff.mess-schiff.orig).^2);
schiff.std_karl = mean((schiff.x-schiff.orig).^2);
schiff.std_orig_dB = 10*log10(schiff.std_orig)*2
schiff.std_karl_dB = 10*log10(schiff.std_karl)*2

% Plots ------------------------------------------------------------------%
% measured position
figure(1)                                   
stairs(mesT,mesP)
title('measured position','Fontsize',12)
ylabel('p_{m} [m]')
xlabel('time [s]')
grid minor;

% figure(2)
% plot(mesT,mesV)
% %axis([ 0 tsim 0 xfp])
% title('measured velocity','Fontsize',12)
% ylabel('v [m/s]')
% xlabel('time [s]')
% grid minor
%%
% measured, calculated and real velocity
f3=figure(3)

f3p2=plot(mesT,mesV,'-','LineWidth',1,'Color',[0, 0.4470, 0.7410]);
hold on
f3p3=plot(mesT,xprp(2,:),'-','LineWidth',1,'Color',[0.8500, 0.3250, 0.0980]);
hold on
f3p1=plot(squeeze(U.time),squeeze(U.signals.values),':k','LineWidth',1);
hold on
%f3p4=stem(mesT(xT==0),xT(xT==0),'-r','LineWidth',0.05); 
%hold on
%f3p4=stem(mesT(xT==1),0*xT(xT==1),'-g','LineWidth',0.05); 
%hold on
f3p1.Color(4)=0.99;
f3p2.Color(4)=0.7;
f3p3.Color(4)=0.8;
hold on

hold off
x1=squeeze(U.signals.values);
x2=xprp(2,:);
l1=xT(xT==0);
tl1=mesT(xT==0);
l2=xT(xT==1);
tl2=mesT(xT==1);
grid minor
ylabel('v [m/s]')
xlabel('time [s]')
lgd = legend({'measured','calculated','real'},'Location','northwest','FontSize',8); %northwest
title(lgd,'velocity:')
% create a new pair of axes inside current figure
% axes('position',[.45 .37 .45 .55])
% box on % put box around new pair of axes
%
% indexOfInterest = (mesT < 2.2) & (mesT > 1.75); % range of t near perturbation
%
% f3xp2=plot(mesT(indexOfInterest),mesV(indexOfInterest),'-','LineWidth',1,'Color',[0, 0.4470, 0.7410]); % plot on new axes
% hold on
% f3xp3=plot(mesT(indexOfInterest),x2(indexOfInterest),'-','LineWidth',2.5,'Color',[0.8500, 0.3250, 0.0980]);
% hold on
% f3xp1=plot(mesT(indexOfInterest),x1(indexOfInterest),':k','LineWidth',2);
% hold on
% % f3xp4=stem(mesT(indexOfInterest),xT(indexOfInterest),'-','LineWidth',0.05,'Color',	[0, 0.5, 0]);%hold on
% % hold on
% % f3xp1.Color(4)=0.99;
% % f3xp2.Color(4)=0.7;
% % f3xp3.Color(4)=0.6;
% axis tight
%%

% measured, calculated and real position
figure(4)
stairs(mesT,mesP);
hold on
f4p1=plot(mesT,xprp(1,:));
hold on
plot(tb,PosReal);
%title('position','Fontsize',12)
grid minor
ylabel('p [m]')
%ylabel({'p [m]'},'Interpreter','latex')
xlabel('Zeit [s]')
lgd=legend('measured','calculated','real','Location','southeast');
title(lgd,'position:')

%%
for i = 1:((length(mesT)-1))
    tS2(i,1)=mesT(i,1);
end;
%%
% Kalman Gains
figure(5)
subplot(2,3,1);
plot(tS2,K1(:,1)')
set(gca,'xtick',1:6); xlim([0 6]); 
title('K11','Fontsize',8)
grid minor
hold on
subplot(2,3,4);
plot(tS2,K1(:,2)')
set(gca,'xtick',1:6); xlim([0 6]);
title('K12','Fontsize',8)
xlabel('Time [s]')
grid minor
subplot(2,3,2);
plot(tS2,K2(:,1)')
set(gca,'xtick',1:6); xlim([0 6]);
title('K21','Fontsize',8)
grid minor
hold on
subplot(2,3,5);
plot(tS2,K2(:,2)')
set(gca,'xtick',1:6); xlim([0 6]);
title('K22','Fontsize',8)
xlabel('Time [s]')
grid minor
hold on
subplot(2,3,3);
plot(tS2,K3(:,1)')
set(gca,'xtick',1:6); xlim([0 6]);
title('K31','Fontsize',8)
grid minor
subplot(2,3,6);
hold on
plot(tS2,K3(:,2)')
set(gca,'xtick',1:6); xlim([0 6]);
title('K32','Fontsize',8)
xlabel('Time [s]')
grid minor
set(gcf,'Position',[649 333 500 250])

%%
% error covariance of estimation
figure(6)
subplot(3,3,1);
plot(tS2,Pp1(:,1)')
title('Pp11','Fontsize',12)
grid minor
subplot(3,3,4);
plot(tS2,Pp2(:,1)')
title('Pp21','Fontsize',12)
grid minor
subplot(3,3,7);
plot(tS2,Pp3(:,1)')
title('Pp31','Fontsize',12)
xlabel('time [s]')
grid minor
subplot(3,3,2);
plot(tS2,Pp1(:,2)')
title('Pp12','Fontsize',12)
grid minor
subplot(3,3,5);
plot(tS2,Pp2(:,2)')
title('Pp22','Fontsize',12)
grid minor
subplot(3,3,8);
plot(tS2,Pp3(:,2)')
xlabel('time [s]')
title('Pp32','Fontsize',12)
grid minor
subplot(3,3,3);
plot(tS2,Pp1(:,3)')
title('Pp13','Fontsize',12)

grid minor
subplot(3,3,6);
plot(tS2,Pp2(:,3)')
title('Pp23','Fontsize',12)
%xlabel('time [s]')
grid minor
subplot(3,3,9);
plot(tS2,Pp3(:,3)')
title('Pp33','Fontsize',12)
xlabel('time [s]')
grid minor