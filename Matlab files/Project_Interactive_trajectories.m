clear all;
close all hidden;
clc;

f2 = figure(2);
global a;
Condition = 1;
switch Condition
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
global v;
v = [150 30];
global rendertime;
rendertime = 0;
global states;
states = [1 1 1];

f = figure(1);
f.Position = [900 100 350 500];
p = uipanel(f,'Position',[0.1 0.1 0.8 0.8]);
c = [];
t = [];
e = [];
r = [];
b = [];
tic;

for i = 1:7
    h = 350-30*i;
    name = strcat('a',num2str(i));
    t = [t; uicontrol(p,'Style','text','String',name,'Position',[10 h 20 20])];
    e = [e; uicontrol(p,'Style','edit','Value',i,'Position',[180 h 80 20])];
    e(i).Max = 10;
    e(i).Min = 0;
    e(i).String = num2str(a(i));
    c = [c; uicontrol(p,'Style','slider','String',name,'Position',[30 h 150 20])];
    c(i).Max = 10;
    c(i).Min = 0;
    c(i).Value = a(i);
    c(i).Callback = @(src,event)function1(src,event,e(i));
    e(i).Callback = @(src,event)function1(src,event,c(i));
end

r = [r uicontrol(p,'Style','radiobutton','String','Fixed points','Callback',@function1)];
r = [r uicontrol(p,'Style','radiobutton','String','Nullclines','Callback',@function1)];
r = [r uicontrol(p,'Style','radiobutton','String','Trajectories','Callback',@function1)];
for i = 1:length(r)
    h = 140-30*i;
    r(i).Position = [10 h 100 20];
    r(i).Value = 1;
end

b = [b uicontrol(p,'Style','pushbutton','String','Re-plot','Callback',@function1)];

function function1(src, ~, h)
global a;
global v;
global rendertime;
global states;
switch src.Style
    case 'edit'
        var = strcat('a',num2str(src.Value));
        val = str2double(src.String);
        h.Value = val;
    case 'slider'
        var = src.String;
        val = src.Value;
        h.String = num2str(val);
    case 'pushbutton'
        var = '';
    case 'radiobutton'
        var = src.String;
        
end
if(toc > rendertime+0.1)
    tic;
    switch var
        case 'a1'
            a(1) = val;
        case 'a2'
            a(2) = val;
        case 'a3'
            a(3) = val;
        case 'a4'
            a(4) = val;
        case 'a5'
            a(5) = val;
        case 'a6'
            a(6) = val;
        case 'a7'
            a(7) = val;
        case 'Fixed points'
            states(1) = ~states(1);
        case 'Nullclines'
            states(2) = ~states(2);
        case 'Trajectories'
            states(3) = ~states(3);
    end
    
    f2 = figure(2);
    if(~isempty(f2.Children))
        v = f2.Children.View;
    end
    hold off;
    
    if states(1)
        syms Ms Hs Rs;
        f_1 = 1 + a(1)*Ms*(1-Ms) - a(2)*Ms*Hs;
        f_2 = a(3)*Hs*Rs - a(4)*Hs;
        f_3 = a(5)*Rs*(1-Rs) - a(6)*Hs*Rs - a(7)*Rs;
        assume([Ms>=0 Hs>=0 Rs>=0]);
        S = vpasolve([f_1==0 f_2==0 f_3==0],[Ms Hs Rs]);
        Fixedpoints = [double(S.Ms) double(S.Hs) double(S.Rs)];
        
        if ~isempty(Fixedpoints)
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
                plot3(R0,H0,M0,'o','MarkerSize',10,'MarkerFaceColor',stability,'MarkerEdgeColor','k');
                view(v(1), v(2));
                hold on;
            end
        end
    end
    Xmax = 3;
    Ymax = 3;
    Zmax = 3;
    if states(3)
        % Initial conditions
        m0 = 0:1:2;
        h0 = m0;
        r0 = m0;
        
        for M0 = m0
            for H0 = h0
                for R0 = r0
                    % Simulation time
                    t0 = 0;
                    te = 1000;
                    
                    % Insert DE
                    f = @(t,X) [1+a(1)*X(1)*(1-X(1))-a(2)*X(1)*X(2); a(3)*X(2)*X(3)-a(4)*X(2); a(5)*X(3)*(1-X(3))-a(6)*X(2)*X(3)-a(7)*X(3)];
                    
                    % Solve DE
                    [t,xa] = ode45(f,[t0 te],[M0 H0 R0]);
                    
                    % Plot DE solution
                    figure(2);
                    plot3(xa(:,3),xa(:,2),xa(:,1),'k');
                    hold on;
                    % Mark first and last value
                    plot3(xa(1,3),xa(1,2),xa(1,1),'.c','MarkerSize',15);
                    plot3(xa(end,3),xa(end,2),xa(end,1),'.b','MarkerSize',15);
                    Xmax = max(xa(:,1));
                    Ymax = max(xa(:,2));
                    Zmax = max(xa(:,3));
                    view(v(1), v(2));
                end
            end
        end
    end
    if states(2)
        R = 0:0.1:Xmax; %Resting cells
        H = 0:0.1:Ymax; %Active hunting cells
        M = 0:0.1:Zmax; %Tumor cells
        
        [R1,H1] = meshgrid(R,H);
        [R2,M1] = meshgrid(R,M);
        [H2,M2] = meshgrid(H,M);
        
        NM1 = (0.5-(a(2).*H1)/(2*a(1)))+sqrt((a(1)^2+a(2)^2.*H1.^2-2*a(1)*a(2).*H1+4*a(1))/(4*a(1)^2));
        NM2 = (0.5-(a(2).*H1)/(2*a(1)))-sqrt((a(1)^2+a(2)^2.*H1.^2-2*a(1)*a(2).*H1+4*a(1))/(4*a(1)^2));
        NH1 = ones(size(R2)).*0;
        NH2 = (a(5)*(1-R2)-a(7))/(a(6));
        NH3 = (-a(1).*M1.*(M1-1)-1)./(a(2).*M1);
        NR1 = ones(size(H2)).*0;
        NR2 = (-a(6).*H2-a(7)+a(5))/(a(5));
        NR3 = ones(size(H2))*a(4)/a(3);
        
        figure(2);
        alpha = 0.2;
        mesh(R1,H1,NM1,'EdgeColor','r','FaceColor','r','FaceAlpha',alpha,'LineStyle','none');
        hold on;
        mesh(R1,H1,NM2,'EdgeColor','r','FaceColor','r','FaceAlpha',alpha,'LineStyle','none');
        mesh(R2,NH1,M1,'EdgeColor','g','FaceColor','g','FaceAlpha',alpha,'LineStyle','none');
        mesh(R2,NH2,M1,'EdgeColor','g','FaceColor','g','FaceAlpha',alpha,'LineStyle','none');
        mesh(R2,NH3,M1,'EdgeColor','g','FaceColor','g','FaceAlpha',alpha,'LineStyle','none');
        mesh(NR1,H2,M2,'EdgeColor','b','FaceColor','b','FaceAlpha',alpha,'LineStyle','none');
        mesh(NR2,H2,M2,'EdgeColor','b','FaceColor','b','FaceAlpha',alpha,'LineStyle','none');
        mesh(NR3,H2,M2,'EdgeColor','b','FaceColor','b','FaceAlpha',alpha,'LineStyle','none');
        view(v(1), v(2));
    end
    if(~isempty(f2.Children))
        f2.Children.XGrid = 'on';
        f2.Children.YGrid = 'on';
        f2.Children.ZGrid = 'on';
        f2.Children.XLabel.String = 'R';
        f2.Children.YLabel.String = 'H';
        f2.Children.ZLabel.String = 'M';
        f2.Children.XLim = [0 Xmax];
        f2.Children.YLim = [0 Ymax];
        f2.Children.ZLim = [0 Zmax];
    end
    rendertime = toc;
end
end