%% -----------------------------------------ASSIGNING VALUES----------------------------------------------------
clear; % Cleaning the previous workspace so that we don't have any conflicts
load('DatasetFourier.mat');

t_id = time(1:88);  % Time used for the first 80% of the values
t_val = time(89:109); % Time used for the last 20% of the values 

y_id = y(1:88); % Product quantity for the first 80% of the values
y_val = y(89:109); % Product quantity for the last 20% of the values

%% ---------------------------------------PLOTTING FOR M = 1:7--------------------------------------------------
for m = 1:7
    phi_id = matrix_id(m); % Function used to calculate the identification matrix
    phi_val = matrix_val(m); % Function used to calculate the validation matrix
    
    theta=linsolve(phi_id,y_id); % Theta is our aproximator,calculated only for the first 80% of the values
  
    yhat_id=phi_id*theta; % Approximation for the first 80% of the values
    yhat_val=phi_val*theta; % Approximation for the last 20% of the values
    
    MSE_id(m)=immse(y_id,yhat_id); % Mean square value for the first 80% of the values
    MSE_val(m)=immse(y_val,yhat_val); % Mean square value for the last 20% of the values
    
    % Plotting the graphics
    % figure;
    % plot(time,y),legend('Real values'),
    % title('m = ',m);xlabel('Time');ylabel('Y');axis([1 109 500 3250]);grid;
    % hold on;
    % plot(t_id,yhat_id,'g');
    % plot(t_val,yhat_val,'r'),legend('Y','Yhat_{id}','Yhat_{val}');
end

%% -----------------------------------------PLOTTING THE MSE----------------------------------------------------
figure
m=1:7;
subplot(211);plot(m,MSE_id);title("MSE_{ID}");xlabel('M');ylabel('MSE');grid;
hold on
subplot(212);plot(m,MSE_val);title("MSE_{VAL}");xlabel('M');ylabel('MSE');grid;
[MSE_val_min,m_optim] = min(MSE_val(m)); % We only care about m_optim, which will be used for the last plot

%% ----------------------------------------------CLEANING-------------------------------------------------------------
clear phi_id phi_val theta yhat_id yhat_val; 
% Cleaning the values of the matrices and the approximator theta so that we can make those for our optim m 

%% ----------------------------------------PLOTTING FOR M OPTIM----------------------------------------------
phi_id = matrix_id(m_optim); % Function used to calculate the identification matrix
phi_val = matrix_val(m_optim); % Function used to calculate the validation matrix

theta=phi_id\y_id; % Theta is our approximator, calculated only for the first 80% of the values

yhat_id=phi_id*theta; % Approximation for the first 80% of the values
yhat_val=phi_val*theta; % Approximation for the last 20% of the values

% Plotting the graphics for m_optim
figure;
plot(time,y),legend('Real values');
title("m_{optim} = 3");xlabel('Time');ylabel('Y');axis([1 109 500 3250]);grid;
hold on;
plot(t_id,yhat_id,'g');
plot(t_val,yhat_val,'r'),legend('Y','Yhat_{id}','Yhat_{val}');

% Plotting the validation graphics m_optim
figure;
plot(time(89:109),y(89:109));
title('m_{optim} = 3');xlabel('Time');ylabel('Y');axis([89 109 600 2500]);grid;
hold on;
plot(t_val,yhat_val,'r'),legend('Y','Yhat_{val}');
%% ---------------------------------------------FUNCTIONS------------------------------------------------------------
function phi_id = matrix_id(m) % Function used to calculate the identification matrix
    load('product_16.mat');
    t_id = time(1:88);  % Time for the first 80% of the values 
    p=12; % Used only in creating the matrix, it's the number of months in a year for our data

    for i=1:length(t_id) % More precisely i = 1->88
        for j=1:2+2*m
            phi_id(i,1)=1;        
            phi_id(i,2)=i;   
            r = mod(j,2);
                if r==1
                    phi_id(i,j)=cos(2*pi*i*(j-1)/2/p);
                else
                    phi_id(i,j)=sin(2*pi*i*(j-2)/2/p);
                end
        end
    end
end

function phi_val = matrix_val(m) % Function used to calculate the validation matrix 
    load('product_16.mat');
    t_val = time(89:109); % Time for the last 20% of the values 
    p=12; % Used only in creating the matrix, it's the number of months in a year for our data

     for i=1:length(t_val) % More precisely i = 1->21
        for j=1:2+2*m
            phi_val(i,1)=1;        
            phi_val(i,2)=t_val(i);   
            r = mod(j,2);
                if r==1
                    phi_val(i,j)=cos(2*pi*t_val(i)*(j-1)/2/p);
                else
                    phi_val(i,j)=sin(2*pi*t_val(i)*(j-2)/2/p);
                end
        end
    end
end
        