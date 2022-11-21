clear all;
close all hidden;
clc;

condition = 3;
switch condition
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

var = 3:0.05:4;

Mv = [];
Hv = [];
Rv = [];

for i = 1:length(var)
    
    a(7) = var(i);
    
    syms M H R;
    f1 = 1+a(1)*M*(1-M)-a(2)*M*H;
    f2 = a(3)*H*R - a(4)*H;
    f3 = a(5)*R*(1-R)-a(6)*H*R-a(7)*R;
    assume([M>=0 H>=0 R>=0]);
    S = vpasolve([f1==0 f2==0 f3==0],[M H R]);
    S.M;
    S.H;
    S.R;
    
    Fixedpoints = [double(S.M) double(S.H) double(S.R)];
    
    for s = 1:size(Fixedpoints,1)
        M0 = Fixedpoints(s,1);
        H0 = Fixedpoints(s,2);
        R0 = Fixedpoints(s,3);
        
        A = [a(1)-2*a(1)*M0-a(2)*H0 -a(2)*M0     0;
            0                      a(3)*R0-a(4) a(3)*H0;
            0                      -a(6)*R0     a(5)-2*a(5)*R0-a(6)*H0-a(7)];
        
        E = eig(A);
        if(max(real(E)) < 0)
            stability = 'g';
        else
            stability = 'r';
        end
        
        figure(1);
        plot3(R0,H0,M0,'.','MarkerSize',10,'MarkerEdgeColor',stability);
        hold on;
        grid on;
        xlabel('R');
        ylabel('H');
        zlabel('M');
        
        figure(2);
        subplot(311);
        plot(var(i),R0,'.','MarkerSize',10,'MarkerEdgeColor',stability);
        hold on;
        ylabel('R*');
        xlabel('bifurcation var');
        subplot(312);
        plot(var(i),H0,'.','MarkerSize',10,'MarkerEdgeColor',stability);
        hold on;
        ylabel('H*');
        xlabel('bifurcation var');
        subplot(313);
        plot(var(i),M0,'.','MarkerSize',10,'MarkerEdgeColor',stability);
        hold on;
        ylabel('M*');
        xlabel('bifurcation var');
    end
end

