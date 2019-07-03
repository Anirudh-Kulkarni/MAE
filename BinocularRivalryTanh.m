% This is a program to implement Binocular rivarly (rivarly of two populations in general)
% It is based on the article by Chaing and Low on binocular rivalry!

% From 50 seconds to 500 seconds, only one population is excited 2 seconds in every consecutrive 10 seconds.
% The idea is to check how adaptation affects the activity of the two populations after this stimulus.

%%% The parameters might have to be adjusted in this program !!!

% Author: Anirudh KULKARNI

function BinocularRivalryTanh % Tanh because the nonlinear function used here is tanh!
    global  dt T  Inext1 Inext2 % Defune global variables which can be used outside this function

    T =1000; % total time of simulation in seconds
    dt = 5 *(10^(-3)); % time step of the simulation
    t = 0:dt:T; % total time duration 

    u1 = zeros(1,length(t)); % Predefine four populations ; two populations u1 and u2 represent firing rates
    u2 = zeros(1,length(t));
    % two other populations a1 and a2 represent the adaptation variables of these populations
    a1 = zeros(1,length(t));
    a2 = zeros(1,length(t));
    
    % adjust the amount of noise in the input current to the populations by changing sigma
    sigma = 1;
    
    % Pulse a is the total duration of the current pulse: on and off
    % pulse b is the duration of the current pulse when the current increases: on
    pulsea = 10; % in seconds!
    pulseb = 2;
        
    % run a for loop 
    for i = 1:length(t)-1
          Inext1 = -3  +  sigma*randn;
          Inext2 = -3  +  sigma*randn; % this is the baseline value of the currents to the two populations
          
          % From 50 seconds to 500 seconds, excite only one population 2 seconds in every 10 seconds
          if i*dt < 500 && i*dt > 50
            if mod(i*dt, pulsea) < pulseb
             Inext1 = Inext1  +  9;
            end
          else
            % Otherwise, excite both populations
              if mod(i*dt, pulsea) < pulseb
                Inext1 = Inext1  +  9;
                Inext2 = Inext2  +  9;
              end
          end
         
        %6 step method of Runge-Kutta Fehlberg to integrate the differential equations in the HH function below
         
          K1 = HH(i,[u1(:,i), u2(:,i), a1(:,i), a2(:,i)]);
          k1 = K1(:,1);
          n1 = K1(:,2);
          h1 = K1(:,3);
          s1 = K1(:,4);% obtain 4 k variables (V,n,h,s) from HH function
          
          K2 = HH(i + 0.2,[u1(:,i) + (0.2*k1*dt),u2(:,i) + (0.2*n1*dt),a1(:,i) + (0.2*h1*dt),a2(:,i) + (0.2*s1*dt)]);
          k2 = K2(:,1);
          n2 = K2(:,2);
          h2 = K2(:,3);
          s2 = K2(:,4);
          
          K3 = HH(i + 0.3,[u1(:,i) + (3*k1*dt/40) + (9*k2*dt/40),u2(:,i) + (3*n1*dt/40) + (9*n2*dt/40),a1(:,i) + (3*h1*dt/40) + (9*h2*dt/40),a2(:,i) + (3*s1*dt/40) + (9*s2*dt/40)]);
          k3 = K3(:,1);
          n3 = K3(:,2);
          h3 = K3(:,3);
          s3 = K3(:,4);
          
          K4 = HH(i + 0.6,[u1(:,i) + 0.3*k1*dt-0.9*k2*dt + 1.2*k3*dt, u2(:,i) + 0.3*n1*dt-0.9*n2*dt + 1.2*n3*dt, a1(:,i) + 0.3*h1*dt-0.9*h2*dt + 1.2*h3*dt, a2(:,i) + 0.3*s1*dt-0.9*s2*dt + 1.2*s3*dt]);
          k4 = K4(:,1);
          n4 = K4(:,2);
          h4 = K4(:,3);
          s4 = K4(:,4);
          
          K5 = HH(i + 1,[u1(:,i)-11*k1*dt/54 + 5*k2*dt/2-70*k3*dt/27 + 35*k4*dt/27, u2(:,i)-11*n1*dt/54 + 5*n2*dt/2-70*n3*dt/27 + 35*n4*dt/27, a1(:,i)-11*h1*dt/54 + 5*h2*dt/2-70*h3*dt/27 + 35*h4*dt/27, a2(:,i)-11*s1*dt/54 + 5*s2*dt/2-70*s3*dt/27 + 35*s4*dt/27]);
          k5 = K5(:,1);
          n5 = K5(:,2);
          h5 = K5(:,3);
          s5 = K5(:,4);
          
          K6 = HH(i + 7/8,[u1(:,i) + 1631*k1*dt/55296 + 175*k2*dt/512 + 575*k3*dt/13824 + 44275*k4*dt/110592 + 253*k5*dt/4096, u2(:,i) + 1631*n1*dt/55296 + 175*n2*dt/512 + 575*n3*dt/13824 + 44275*n4*dt/110592 + 253*n5*dt/4096, a1(:,i) + 1631*h1*dt/55296 + 175*h2*dt/512 + 575*h3*dt/13824 + 44275*h4*dt/110592 + 253*h5*dt/4096, a2(:,i) + 1631*s1*dt/55296 + 175*s2*dt/512 + 575*s3*dt/13824 + 44275*s4*dt/110592 + 253*s5*dt/4096]);
          k6 = K6(:,1);
          n6 = K6(:,2);
          h6 = K6(:,3);
          s6 = K6(:,4);
          
    %create next step for each variable
    u1(:,i + 1) = u1(:,i) + (37*k1/378 + 250*k3/621 + 125*k4/594 + 512*k6/1771)*dt;
    u2(:,i + 1) = u2(:,i) + (37*n1/378 + 250*n3/621 + 125*n4/594 + 512*n6/1771)*dt;
    a1(:,i + 1) = a1(:,i) + (37*h1/378 + 250*h3/621 + 125*h4/594 + 512*h6/1771)*dt;
    a2(:,i + 1) = a2(:,i) + (37*s1/378 + 250*s3/621 + 125*s4/594 + 512*s6/1771)*dt;

    end

    %% Plots of the firing rate activities and the adaptation variables of the two populations

    figure
    plot(t,u1);
    hold
    plot(t,-u2,'r');
     
    figure
    plot(t,a1);
    hold
    plot(t,-a2,'r');
    
end


function outf = f(v)%Alpha for variable K_n

    outf = (1 + tanh(v))/2;  
     
end


function dydt = HH(~,y)
  % Here we define the model for the differential equations
    global Inext1 Inext2;
    u1 = y(:,1);
    u2 = y(:,2);
    a1 = y(:,3);
    a2 = y(:,4);

    beta = 10;
    TauA = 100;
    Tau = 200*(10^(-3));
    PhiA = 1;

    dydt = [(1/Tau)*(-u1  +  f(-beta*u2-a1  +  Inext1)),(1/Tau)*(-u2  +  f(-beta*u1-a2  +  Inext2)),((1/TauA)*(-a1 + PhiA*f(-beta*u2-a1 + Inext1))),((1/TauA)*(-a2 + PhiA*f(-beta*u1-a2 + Inext2)))];
 end
