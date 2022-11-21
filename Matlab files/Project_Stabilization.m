clear all;
close all hidden;
clc;

a = [3 1 4.8 0.4 3.7 1.9 0.1]; %Case 3
us = 1.73646;
syms M H R;

f1 = 1+a(1)*M*(1-M)-a(2)*M*H;
f2 = a(3)*H*R - a(4)*H;
f3 = a(5)*R*(1-R)-a(6)*H*R-a(7)*R+us;
assume([M>=0 H>=0 R>=0]);
S = vpasolve([f1==0 f2==0 f3==0],[M H R]);
S.M;
S.H;
S.R;

Fixedpoints = [double(S.M) double(S.H) double(S.R)];

stability = [];
for s = 1:size(Fixedpoints,1)
M0 = Fixedpoints(s,1);
H0 = Fixedpoints(s,2);
R0 = Fixedpoints(s,3);

A = [a(1)-2*a(1)*M0-a(2)*H0 -a(2)*M0     0;
     0                      a(3)*R0-a(4) a(3)*H0;
     0                      -a(6)*R0     a(5)-2*a(5)*R0-a(6)*H0-a(7)];
 
E = eig(A);
if(max(real(E)) < 0)
    stability = [stability; 1];
    Mf = M0;
    Hf = H0;
    Rf = R0;
else
    stability = [stability; 0];
end
end

Fixedpoints = [Fixedpoints stability];

disp('Fixedpoints');
disp('|    m    |    h    |    r    |stability|');
disp(Fixedpoints);

dt = 0.001;
t = 0:dt:100;

%      |M - M*|
%u = F*|H - H*| + us
%      |R - R*|
F_vect = [1 0 -1];

for F = F_vect'
    A2 = [a(1)-2*a(1)*Mf-a(2)*Hf -a(2)*Mf       0;
          0                      a(3)*Rf-a(4)   a(3)*Hf;
          F(1)                   -a(6)*Rf+F(2)  a(5)-2*a(5)*Rf-a(6)*Hf-a(7)+F(3)];
    E2 = eig(A2);
    if(max(real(E2)) < 0)
        disp(strcat("F= ",num2str(F')," :Assymptotically stable"));
    else
        disp(strcat("F= ",num2str(F')," :NOT assymptotically stable"));
    end
    
    if(us >= F(1)*Mf + F(2)*Hf + F(3)*Rf && F(1)>=0 && F(2) >= 0)
        disp("Positive orthant P remains positive invariant");
    else
        disp("Positive orthant P does NOT remain positive invariant");
    end
end

%IC_vect = [1 1 1; 0 1 0; 10 1 0; 5 2 15];
IC_vect = [1 1 1];

for IC = IC_vect'
    for F = F_vect'
        M = zeros(1,length(t));
        H = zeros(1,length(t));
        R = zeros(1,length(t));
        u = zeros(1,length(t));

        M(1) = IC(1); %density Tumor cells
        H(1) = IC(2); %density Active cells
        R(1) = IC(3); %density Resting cells

        for i = 1:length(t)-1
            M_cur = M(i);
            H_cur = H(i);
            R_cur = R(i);

            u(i) = F'*[M_cur - Mf; H_cur - Hf; R_cur - Rf] + us;

            M_dot = 1+a(1)*M_cur*(1-M_cur)-a(2)*M_cur*H_cur;
            H_dot = a(3)*H_cur*R_cur-a(4)*H_cur;
            R_dot = a(5)*R_cur*(1-R_cur)-a(6)*H_cur*R_cur-a(7)*R_cur + u(i);

            M(i+1) = M_cur+M_dot*dt;
            H(i+1) = H_cur+H_dot*dt;
            R(i+1) = R_cur+R_dot*dt;
        end

        figure(1);
        plot(t,M,'r-');
        hold on;
        plot(t,H,'g-');
        plot(t,R,'b-');
        grid on;
        xlabel('Time');
        ylabel('Density');
        legend('Tumor cells', 'Active cells', 'Resting cells');
        title('Time plot condition 3');

        figure(2);
        plot3(R',H',M','k.','MarkerSize',0.3);
        hold on;
        plot3(R(1),H(1),M(1),'g.','MarkerSize',30);
        plot3(R(end),H(end),M(end),'r.','MarkerSize',30);
        grid on;
        xlabel('R');
        ylabel('H');
        zlabel('M');
        legend('Trajectory','Start point','End point','Location','best');
        title('State plot condition 3');
    end
end

zlim([0 inf])

