%% exercise 6 Lyapunov function
clear all; close all; clc;

% Finding a suitable matrix P = P' > 0 such that the quadratic function
% (see exercise 6) serves as a Lyapunov function to prove the stability of
% the fixed point E3 with the parameter values of Condition 3. 

% Condition 3:
a1 = 3;
a2 = 1;
a3 = 4.8;
a4 = 0.4;
a5 = 3.7;
a6 = 1.9;
a7 = 0.1;

% fixed point
R3 = a4/a3; %0.0833
H3 = (a5*(1-R3)-a7)/(a6); %1.7325
M3 = 0.5-((a2)/(2*a1))*H3+sqrt((a1^2+a2^2*H3^2-2*a1*a2*H3+4*a1)/(4*a1^2)); %0.8260

i = 0;
%x = [0:0.1:10, 11:1:100, 110:10:1000];
x = [0:0.3:10];
Vf = zeros(1,length(x));
for m = x
    for h = x
        for r = x
            
            % insert any initial condition >= 0
            M = m;
            H = h;
            R = r;
            
            %
            A = [M-M3; H-H3; R-R3]; 
            
            %
            P = eye(3);
            % P is positive definite because the first value p11 is positive (1) and
            % all eigenvalues are positive (1).
            
            % Lyapunov function
            V = A'*P*A;
                        
            i = i + 1;
            Vf(i) = V;
        end
    end
end

figure(1);
plot(0:length(Vf)-1,Vf);

% check if there is any value for V that is negative or zero
find(Vf <= 0)

% Lyapunov function check:
% * V(x*) = 0, V(x) > 0 (so V is 0 for a fixed point and positive for any
% point): check.
% * V is continuous on S: check
% * V is differentiable along trajectories: check

r = 0:0.3:10;
h = 0:0.3:10;
m = 0:0.3:10;

% figure(2);
% for R = r
%     for H = h
%         Vf_diff_R = 2*(R-R3)*(a5*R*(1-R)-a6*H*R-a7*R);
%         if Vf_diff_R < 0 
%             col = [0,0,1];
%         else
%             col = [1,0,0];
%         end
%         plot3(R, H, Vf_diff_R,'.','color',col);
%         hold on;
%     end
% end
% plot3(R3, H3, 0,'.','color',[0,1,0]);
% grid on;
% xlabel('R');
% ylabel('H');
% zlabel('Vf`')
% 
% figure(3);
% for R = r
%     for H = h
%         Vf_diff_H = 2*(H-H3)*(a3*H*R-a4*H);
%         if Vf_diff_H < 0 
%             col = [0,0,1];
%         else
%             col = [1,0,0];
%         end
%         plot3(R, H, Vf_diff_H,'.','color',col);
%         hold on;
%     end
% end
% plot3(R3, H3, 0,'.','color',[0,1,0]);
% grid on;
% xlabel('R');
% ylabel('H');
% zlabel('Vf`')
% 
% figure(4);
% for M = m
%     for H = h
%         Vf_diff_M = 2*(M-M3)*(1+a1*M*(1-M)-a2*M*H);
%         if Vf_diff_M < 0 
%             col = [0,0,1];
%         else
%             col = [1,0,0];
%         end
%         plot3(M, H, Vf_diff_M,'.','color',col);
%         hold on;
%     end
% end
% plot3(M3, H3, 0,'.','color',[0,1,0],'Markersize',30);
% grid on;
% xlabel('M');
% ylabel('H');
% zlabel('Vf`')

figure(5);
for M = m
    for R = r
        for H = h
            Vf_diff = ...
                1*(R-R3)*(a5*R*(1-R)-a6*H*R-a7*R) +...
                1*(H-H3)*(a3*H*R-a4*H) + ...
                1*(M-M3)*(1+a1*M*(1-M)-a2*M*H);
            if Vf_diff <= 0 
                col = [0,0,1];
            else
                col = [1,0,0];
                plot3(R, H, M,'.','color',col);
                hold on;
            end
            %plot3(R, H, M,'.','color',col);
            %hold on;
        end
    end
end
plot3(R3, H3, M3,'.','color',[0,1,0],'Markersize',30);
grid on;
xlabel('R');
ylabel('H');
zlabel('M');
xlim([0 max(r)]);
ylim([0 max(h)]);
zlim([0 max(m)]);

%% NON LINEARIZED LYAPUNOV
clear all; close all; clc;

% Condition 3:
a1 = 3;
a2 = 1;
a3 = 4.8;
a4 = 0.4;
a5 = 3.7;
a6 = 1.9;
a7 = 0.1;

% fixed point
R3 = a4/a3; %0.0833
H3 = (a5*(1-R3)-a7)/(a6); %1.7325
M3 = 0.5-((a2)/(2*a1))*H3+sqrt((a1^2+a2^2*H3^2-2*a1*a2*H3+4*a1)/(4*a1^2)); %0.8260

% Linearised A matrix
A = [a1-2*a1*M3-a2*H3       -a2*M3      0;
     0                      a3*R3-a4    a3*H3;
     0                      -a6*R3      a5-2*a5*R3-a6*H3-a7];

Q = [1 0 0;
     0 1 0;
     0 0 1];

P = lyap(A',Q);
 
p11 = P(1,1);
p22 = P(2,2);
p33 = P(3,3);
p12 = P(1,2);
p13 = P(1,3);
p23 = P(2,3);

radius = 0.002;
r_vect = 0:radius/1000:radius;
angle_vect = 0:2*pi/120:2*pi;
n_points = length(r_vect)*length(angle_vect)^2;

V = zeros(n_points,1);
V_diff = zeros(n_points,1);

M_vect_valid = zeros(length(angle_vect)^2,1);
H_vect_valid = zeros(length(angle_vect)^2,1);
R_vect_valid = zeros(length(angle_vect)^2,1);
M_vect_invalid = zeros(length(angle_vect)^2,1);
H_vect_invalid = zeros(length(angle_vect)^2,1);
R_vect_invalid = zeros(length(angle_vect)^2,1);

M_dot_vect = zeros(n_points,1);
H_dot_vect = zeros(n_points,1);
R_dot_vect = zeros(n_points,1);

M_diff_vect = zeros(n_points,1);
H_diff_vect = zeros(n_points,1);
R_diff_vect = zeros(n_points,1);

i = 0;
rmax = 0;
boolmax = false;
for r = r_vect
    j = 0;
    for phi = angle_vect
        for theta = angle_vect
            i = i + 1;
            j = j + 1;

            R = R3 - r*cos(phi)*cos(theta);
            H = H3 - r*sin(phi)*cos(theta);
            M = M3 - r*sin(theta);
            
            M_dot = 1 + a1*M*(1-M) - a2*M*H;
            H_dot = a3*H*R - a4*H;
            R_dot = a5*R*(1-R) - a6*H*R-a7*R;

            V(i) = p11*(M-M3)^2 + p22*(H-H3)^2 + p33*(R-R3)^2 ...
                 + 2*p12*((H-H3)*(M-M3)) + 2*p13*((R-R3)*(M-M3)) + 2*p23*((R-R3)*(H-H3));
            
            V_diff(i) = M_dot*(2*p11*(M-M3) + 2*p12*(H-H3) + 2*p13*(R-R3)) + ...
                        H_dot*(2*p12*(M-M3) + 2*p22*(H-H3) + 2*p23*(R-R3)) + ...
                        R_dot*(2*p13*(M-M3) + 2*p23*(H-H3) + 2*p33*(R-R3));
                    
            M_dot_vect(i) = M_dot;
            H_dot_vect(i) = H_dot;
            R_dot_vect(i) = R_dot;
            
            M_diff_vect(i) = M_dot*(2*p11*(M-M3) + 2*p12*(H-H3) + 2*p13*(R-R3));
            H_diff_vect(i) = H_dot*(2*p12*(M-M3) + 2*p22*(H-H3) + 2*p23*(R-R3));
            R_diff_vect(i) = R_dot*(2*p13*(M-M3) + 2*p23*(H-H3) + 2*p33*(R-R3));
            
            cond_fp = (V(i) >= 0 && V_diff(i) <= 0) && (R == R3 && H == H3 && M == M3);
            cond_notfp = (V(i) > 0 && V_diff(i) <= 0);

            if((cond_fp || cond_notfp))
                M_vect_valid(j) = M;
                H_vect_valid(j) = H;
                R_vect_valid(j) = R;
            else
                M_vect_invalid(j) = M;
                H_vect_invalid(j) = H;
                R_vect_invalid(j) = R;
                boolmax = true;
            end
        end
    end
    if(~boolmax)
        rmax = r;
    else
        break;
    end
end

disp(strcat("Maximum radius = ",num2str(rmax)));

M_vect_valid(M_vect_valid==0) = [];
H_vect_valid(H_vect_valid==0) = [];
R_vect_valid(R_vect_valid==0) = [];
M_vect_invalid(M_vect_invalid==0) = [];
H_vect_invalid(H_vect_invalid==0) = [];
R_vect_invalid(R_vect_invalid==0) = [];

figure(1);
plot3(R_vect_valid,H_vect_valid,M_vect_valid,'g.');
hold on;
plot3(R_vect_invalid,H_vect_invalid,M_vect_invalid,'r.','Markersize',30);
plot3(R3,H3,M3,'k.','Markersize',30);
grid on;
xlabel('R');
ylabel('H');
zlabel('M');
title("Lyapunov validation");
xlim([R3-radius R3+radius]);
ylim([H3-radius H3+radius]);
zlim([M3-radius M3+radius]);
legend("Valid Lyapunov points","Invalid Lyapunov point(s)","Fixed point",'Location','best');

% figure(2);
% plot(0:n_points-1,V);
% hold on;
% grid on;
% title("V");

% figure(3);
% plot(0:n_points-1,V_diff);
% hold on;
% grid on;
% title("V'");

% figure(4);
% plot(0:n_points-1,M_diff_vect);
% hold on;
% plot(0:n_points-1,H_diff_vect);
% plot(0:n_points-1,R_diff_vect);
% grid on;
% legend('Mdiff','Hdiff','Rdiff');

% figure(5);
% plot(0:n_points-1,M_dot_vect);
% hold on;
% plot(0:n_points-1,H_dot_vect);
% plot(0:n_points-1,R_dot_vect);
% grid on;
% legend('Mdot','Hdot','Rdot');


%% LINEARIZED LYAPUNOV
clear all; close all; clc;

% Condition 3:
a1 = 3;
a2 = 1;
a3 = 4.8;
a4 = 0.4;
a5 = 3.7;
a6 = 1.9;
a7 = 0.1;

% fixed point
R3 = a4/a3; %0.0833
H3 = (a5*(1-R3)-a7)/(a6); %1.7325
M3 = 0.5-((a2)/(2*a1))*H3+sqrt((a1^2+a2^2*H3^2-2*a1*a2*H3+4*a1)/(4*a1^2)); %0.8260

A = [a1-2*a1*M3-a2*H3       -a2*M3      0;
     0                      a3*R3-a4    a3*H3;
     0                      -a6*R3      a5-2*a5*R3-a6*H3-a7];

Q = [1 0 0;
     0 1 0;
     0 0 1];
 
eigenQ = eig(Q);
posdefQ = true;
for i = 1:length(eigenQ)
    if(real(eigenQ(i)) < 0)
        posdefQ = false;
    end
end
if(posdefQ)
    disp("Q is positive definite");
else
    disp("Q is not positive definite !!!!!!!!!!!!!!!");
    disp(eig(Q));
end

if(posdefQ)
P = lyap(A',Q);

% p11 = P(1,1);
% p22 = P(2,2);
% p33 = P(3,3);
% p12 = P(1,2);
% p13 = P(1,3);
% p23 = P(2,3);

eigenP = eig(P);
posdefP = true;
for i = 1:length(eigenP)
    if(real(eigenP(i)) < 0)
        posdefP = false;
    end
end
if(posdefP)
    disp("P is positive definite");
else
    disp("P is not positive definite !!!!!!!!!!!!!!!");
    disp(eig(P));
end

if(posdefP)
radius = 10;
r_vect = 0:radius/1000:radius;
angle_vect = 0:2*pi/60:2*pi;
n_points = length(r_vect)*length(angle_vect)^2;

V = zeros(n_points,1);
V_diff = zeros(n_points,1);

M_vect = zeros(n_points,1);
H_vect = zeros(n_points,1);
R_vect = zeros(n_points,1);

i = 0;
rmax = 0;
boolmax = false;
for r = r_vect
    for phi = angle_vect
        for theta = angle_vect
            i = i + 1;

            R = R3 - r*cos(phi)*cos(theta);
            H = H3 - r*sin(phi)*cos(theta);
            M = M3 - r*sin(theta);

            V(i) = [M H R]*P*[M H R]';
            
            V_diff(i) = [M H R]*(A'*P+P*A)*[H H R]';
            
            cond_fp = (V(i) >= 0 && V_diff(i) <= 0) && (R == R3 && H == H3 && M == M3);
            cond_notfp = (V(i) > 0 && V_diff(i) <= 0);

            if(~(cond_fp || cond_notfp))
                M_vect(i) = M;
                H_vect(i) = H;
                R_vect(i) = R;
                boolmax = true;
            end 
        end
    end
    if(~boolmax)
        rmax = r;
    end
end

disp(strcat("Maximum radius = ",num2str(rmax)));

M_vect(M_vect==0) = [];
H_vect(H_vect==0) = [];
R_vect(R_vect==0) = [];

figure(1);
plot3(R_vect,H_vect,M_vect,'r.');
hold on;
plot3(R3,H3,M3,'k.','Markersize',30);
grid on;
xlabel('R');
ylabel('H');
zlabel('M');
title("Lyapunov points");
xlim([R3-radius R3+radius]);
ylim([H3-radius H3+radius]);
zlim([M3-radius M3+radius]);

figure(2);
plot(0:n_points-1,V);
hold on;
grid on;
title("V");

figure(3);
plot(0:n_points-1,V_diff);
hold on;
grid on;
title("V'");
end
end


