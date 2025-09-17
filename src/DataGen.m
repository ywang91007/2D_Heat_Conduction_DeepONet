clear all; close all; clc;
% Two Dimension Heat Conduction, Central Finite Difference Method, Unsteady, Thermal Conductivity, Curvature
% Method: https://core.ac.uk/download/pdf/11034604.pdf
%% Control Center (Settings)

fn = {'None','Plot Each Iteration','Save Training Data','Save Testing/Validation Data'};
response=listdlg('PromptString',{'Select which features you would like enabled.','Multiple can be selected.',''},'SelectionMode','multiple','ListString',fn);

if ismember(0,response)==1
    response=[0];
end
%% Set-Up
%%%%%%%%%%%%%%%%%
% Initilization %
%%%%%%%%%%%%%%%%%

%Grid Size
Nx = 100;
Ny = 100; 
%Bolt Parameters: Square
x_k2 = 20; %X-cord
y_k2 = 20; %Y-Cord
l_k2 = 10; %Distance from center to edges
%K-Value [50 < k2 < 200];
k1 = 100;
k2 = 50;
K = k1*ones(Ny,Nx);
K(y_k2-l_k2:y_k2+l_k2-1, x_k2-l_k2:x_k2+l_k2-1) = k2; %K-Value Map
%Temperature Profile: Polynomial or Linear
a = -(0.05:0.01:0.08);
b = 4:8;
c = 300:50:750;
y = 1:Ny;
N = 1;
s = 200;
T_in = zeros(s,Ny);
for i = 1:length(a)
    for j = 1:length(b)
        for k = 1:length(c)
            T_in(N,:) =  a(i)*y.^2+b(j)*y+c(k);
            N = N+1;
        end
    end
end
%T_in = 300;         %Linear Profile
T = zeros(Ny,Nx);
T(:,1) = T_in(1,:)';
T_plus=T;
T_out = zeros(Ny,s);
%Delta Field - Change of temperature per iteration
dx = zeros(Ny,Nx);
tolerance = 5e-2;

%% Code

%%%%%%%%%%%%%%%
% Calculation %
%%%%%%%%%%%%%%%

Max_step = 1e5; %Maximum iteration length
iter = 1000; %Step amount: Change to lengthen simulation [iter < Max_step]
delta_r = 4; %Represents delta*r^2 term

% Training Data
train_bolt = zeros(s,2);
for m = 1:s
    T = zeros(Ny,Nx);
    T(:,1) = T_in(m,:)';
    T_plus=T;
    for n = 1:iter
        for i = 1:Ny
            for j = 2:Nx-1
                % X-Direction:
                x1 = T(i,j)*[K(i,j-1)+2*K(i,j)+K(i,j+1)]/4; %Center Node
                x0 = T(i,j-1)*[K(i,j-1)+K(i,j)]/2; %Left Node
                x2 = T(i,j+1)*[K(i,j+1)+K(i,j)]/2; %Right Node
                dx = (x0-2*x1+x2)/(K(i,j)*(delta_r)); %Differential
            
                T_plus(i,j) = T(i,j)+dx; %Update T+1 model
          
            end
            error(i,Nx) = abs(T_plus(i,Nx) - T(i,Nx)); %Useful for steady state calculations
        end
        T = T_plus; %Save T+1 model as T model
        if n == Max_step
            disp('Max # of iterations reached. Calculation not converging.')
            break
        end
    end
    T_out(:,m) = T(:,x_k2+20);
    train_bolt(m,:) = [x_k2,y_k2]/100;
end

% Test/Validation Data

a1 = [-0.03 -0.04];
b1 = 6;
c1 = [200 250 800 850 900];
N = 1;
t = 10;
T_in_test = zeros(t,Ny);
T_out_test = zeros(Ny,t);
for i = 1:length(a1)
    for j = 1:length(c1)
        T_in_test(N,:) = a1(i)*y.^2+b1*y+c1(j);
        N = N + 1;
    end
end

for m = 1:t
    T = zeros(Ny,Nx);
    T(:,1) = T_in_test(m,:)';
    T_plus=T;
    for n = 1:iter
        for i = 1:Ny
            for j = 2:Nx-1
                % X-Direction:
                x1 = T(i,j)*[K(i,j-1)+2*K(i,j)+K(i,j+1)]/4; %Center Node
                x0 = T(i,j-1)*[K(i,j-1)+K(i,j)]/2; %Left Node
                x2 = T(i,j+1)*[K(i,j+1)+K(i,j)]/2; %Right Node
                dx = (x0-2*x1+x2)/(K(i,j)*(delta_r)); %Differential
            
                T_plus(i,j) = T(i,j)+dx; %Update T+1 model
          
            end
            error(i,Nx) = abs(T_plus(i,Nx) - T(i,Nx)); %Useful for steady state calculations
        end
        T = T_plus; %Save T+1 model as T model
        if n == Max_step
            disp('Max # of iterations reached. Calculation not converging.')
            break
        end
    end
    T_out_test(:,m) = T(:,x_k2+20);
end

%% Plotting
if ismember(2,response) ==1
    figure(1)
    imagesc(T_plus)
    title('Temperature Field');
    xlabel('x');
    ylabel('y');

    figure(2)
    plot(y,T_in(100,:))
    hold on
    plot(y,T_out(:,100))
end

%% Data Export

if ismember(3,response) == 1
    train_branch_input = zeros(8000, 100);
    train_branch_input_norm = train_branch_input;
    train_bolt_input = zeros(8000,2);
    train_trunk_input = zeros(8000, 1);
    train_output = zeros(8000, 1);
    n = 1;
    for i = 1:s
        train_branch_input(n:n+39, 1:100) = repmat(T_in(i,:),40,1);
        %branch_input(n:n+39,101:102) = repmat([x_k2,y_k2],40,1);
        train_bolt_input(n:n+39,:) = repmat(train_bolt(s,:),40,1);
        r = randi(100,40,1);
        train_trunk_input(n:n+39,:) = r;
        train_output(n:n+39,:) = T_out(r,i);
        n = n + 40;
        if i == s
            branch_input_norm_p1 = train_branch_input(:,1:100)/max(train_branch_input(:,1:100),[],'all');
            train_branch_input_norm(:,1:100) = branch_input_norm_p1;
            %branch_input_norm(:,101:102) = branch_input(:,101:102);
        end
    end
    writematrix(train_branch_input_norm, 'train_branch_input_norm.txt');
    writematrix(train_bolt_input,'train_bolt_input.txt')
    writematrix(train_trunk_input, 'train_trunk_input.txt');
    writematrix(train_output, 'train_output.txt');
end

if ismember(4,response) == 1
    
    test_branch_input = T_in_test/max(T_in_test,[],'all');
    %test_branch_input(:,101:102) = repmat([x_k2,y_k2],10,1);
    test_bolt_input = repmat([0.2 0.2],10,1); 
    test_trunk_input = y';
    
    test_output = T_out_test;
    writematrix(test_branch_input, 'test_branch_input_norm.txt');
    writematrix(test_bolt_input, 'test_bolt_input.txt')
    writematrix(test_trunk_input, 'test_trunk_input.txt');
    writematrix(test_output, 'test_output.txt');
end
