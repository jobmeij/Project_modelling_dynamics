clear all;
close all hidden;
clc;

Condition = 1:3;
M0 = [0;3;3];
H0 = [1;4;8];
R0 = [1;4;4]; 
for c = Condition
    switch c
        case 1
            a = [4.5 1 0.5 0.4 0.7 0.2 0.4];
        case 2
            a = [3.5 1 5 0.4 0.7 0.1 0.1];
        case 3
            a = [3 1 4.8 0.4 3.7 1.9 0.1];
        otherwise
            disp('Invalid condition, choose 1, 2 or 3.');
            return;
    end

    for ic = 1:size(M0,1)
        
        dt = 0.01;
        t = 0:dt:100;
        
        M = zeros(1,length(t));
        H = zeros(1,length(t));
        R = zeros(1,length(t));
        
        M(1) = M0(ic); %density Tumor cells
        H(1) = H0(ic); %density Active cells
        R(1) = R0(ic); %density Resting cells

        for i = 1:length(t)-1;
            M_cur = M(i);
            H_cur = H(i);
            R_cur = R(i);

            M_dot = 1+a(1)*M_cur*(1-M_cur)-a(2)*M_cur*H_cur;
            H_dot = a(3)*H_cur*R_cur-a(4)*H_cur;
            R_dot = a(5)*R_cur*(1-R_cur)-a(6)*H_cur*R_cur-a(7)*R_cur;

            M(i+1) = M_cur+M_dot*dt;
            H(i+1) = H_cur+H_dot*dt;
            R(i+1) = R_cur+R_dot*dt;
        end

        style = '-';
        switch ic
            case 1
                style = '-';
            case 2
                style = '-.';
            case 3
                style = '--';
            otherwise
                disp('Invalid condition.');
                return;
        end

        figure(c);
        hold on;
        plot(t,M,strcat('r',style));
        hold on;
        plot(t,H,strcat('g',style));
        plot(t,R,strcat('b',style));
        grid on;
        xlabel('Time');
        ylabel('Density');
        legend('Tumor cells', 'Active cells', 'Resting cells');
        title(strcat('Time plot condition ',num2str(c)));
    end
end

%% TIME DELAY

clear all;
close all hidden;
clc;

a = [3 1 4.8 0.4 3.7 1.9 0.1]; % Case 3

R1 = 0;
H1 = 0;
M1 = 0.5+0.5*sqrt((a(1)+4)/a(1));

R2 = (a(5)-a(7))/a(5);
H2 = 0;
M2 = 0.5+0.5*sqrt((a(1)+4)/a(1));

R3 = a(4)/a(3); %0.0833
H3 = (a(5)*(1-R3)-a(7))/(a(6)); %1.7325
M3 = 0.5-((a(2))/(2*a(1)))*H3+sqrt((a(1)^2+a(2)^2*H3^2-2*a(1)*a(2)*H3+4*a(1))/(4*a(1)^2)); %0.8260

FP_vect = [M3 H3 R3;0 H3 R3;0 1 0.4;M3 H3+0.5 R3;2 4 0.6];
%FP_vect = [M1 H1 R1;M2 H2 R2;M3 H3 R3];
       
delays = [0.7];

dt = 0.002;
t = 0:dt:500;

Ntrajectories = size(FP_vect,1)*length(delays);

M = zeros(Ntrajectories,length(t));
H = zeros(Ntrajectories,length(t));
R = zeros(Ntrajectories,length(t));

j = 0;

for FP = FP_vect'
    for d = 1:size(delays,1)
        j = j + 1;
        delay = delays(d);

        M(j,1) = FP(1); %density Tumor cells
        H(j,1) = FP(2); %density Active cells
        R(j,1) = FP(3); %density Resting cells

        for i = 1:length(t)-1;
            M_cur = M(j,i);
            H_cur = H(j,i);
            R_cur = R(j,i);

            if(i > round(delay/dt,0))
            H_cur_delay = H(j,i-round(delay/dt,0));
            R_cur_delay = R(j,i-round(delay/dt,0));
            else
            H_cur_delay = H(j,1);
            R_cur_delay = R(j,1);
            end

            M_dot = 1+a(1)*M_cur*(1-M_cur)-a(2)*M_cur*H_cur;
            H_dot = a(3)*H_cur_delay*R_cur_delay-a(4)*H_cur;
            R_dot = a(5)*R_cur*(1-R_cur)-a(6)*H_cur*R_cur-a(7)*R_cur;

            M(j,i+1) = M_cur+M_dot*dt;
            H(j,i+1) = H_cur+H_dot*dt;
            R(j,i+1) = R_cur+R_dot*dt;
        end

        style = '-';
        switch d
            case 1
                style = '-';
                col = 'k';
            case 2
                style = '-.';
                col = 'b';
            case 3
                style = '--';
                col = 'm';
            otherwise
                disp('Invalid condition.');
                return;
        end
    end
end 

style = '-';
col = 'k';

figure(1);
hold on;
plot(t,M,strcat('r',style));
hold on;
plot(t,H,strcat('g',style));
plot(t,R,strcat('b',style));
grid on;
xlabel('Time');
ylabel('Density');
legend('Tumor cells', 'Active cells', 'Resting cells');
title('Time plot condition 3');

figure(2);
Tr = plot3(R',H',M',strcat(col,'.'),'MarkerSize',0.1);
hold on;
Ep = plot3(R(:,end),H(:,end),M(:,end),'r.','MarkerSize',30);
hold on;
Sp = plot3(R(:,1),H(:,1),M(:,1),'g.','MarkerSize',20);
grid on;
xlabel('R');
ylabel('H');
zlabel('M');
title('State plot condition 3');
legend([Ep Sp],{'End point','Start point'},'Location','best');

% figure(3);
% plot3(t',R',H',strcat(col,'.'),'MarkerSize',0.3);
% hold on;
% plot3(t(1),R(1),H(1),'g.','MarkerSize',30);
% plot3(t(end),R(end),H(end),'r.','MarkerSize',30);
% grid on;
% xlabel('Time');
% ylabel('R');
% zlabel('H');
% title('State/Time plot condition 3');
% 
% figure(4);
% plot3(t',R',M',strcat(col,'.'),'MarkerSize',0.3);
% hold on;
% plot3(t(1),R(1),M(1),'g.','MarkerSize',30);
% plot3(t(end),R(end),M(end),'r.','MarkerSize',30);
% grid on;
% xlabel('Time');
% ylabel('R');
% zlabel('M');
% title('State/Time plot condition 3');
% 
% figure(5);
% plot3(t',H',M',strcat(col,'.'),'MarkerSize',0.3);
% hold on;
% plot3(t(1),H(1),M(1),'g.','MarkerSize',30);
% plot3(t(end),H(end),M(end),'r.','MarkerSize',30);
% grid on;
% xlabel('Time');
% ylabel('H');
% zlabel('M');
% title('State/Time plot condition 3');
