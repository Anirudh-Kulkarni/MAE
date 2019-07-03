% This program simulates the MAE model that is presented in the paper
% "Sustained Rhythmic Brain Activity Underlies Visual Motion Perception in Zebrafish , 
% Perez-Schuster V., Kulkarni A, et al.
% Cell Reports, Volume 17, Issue 4, 18 October 2016, Pages 1098-1112. "


% Author: Anirudh Kulkarni

function MAE_Simulations
  
   % run two for loops to generate multiple datasets
    
    global  dt T  Ip1 Ip2 IDS1 IDS2 Jd Jc ds1 ds2 

    T =14400; % total time of the simulation
    dt = 0.5 *(10^(-2)); % time step
    t = 0:dt:T; % total time duration variable

    v1 = 0.71*ones(1,length(t)); % the direction selective population 1 firing rate
    v2 = 0.71*ones(1,length(t)); % the direction selective population 2 firing rate
    u1 = 0.1*ones(1,length(t)); % the Comparator population 1 firing rate
    u2 = 0*ones(1,length(t)); % the comparator population 2 firing rate
    
    a1 = 0.92*ones(1,length(t)); % Fast adaptation variable of direction selective population 1
    a2 = 0.92*ones(1,length(t)); % Fast adaptation variable of direction selective population 2
    a3 = 0.92*ones(1,length(t)); % Slow adaptation variable of direction selective population 1
    a4 = 0.92*ones(1,length(t)); % Slow adaptation variable of direction selective population 2
    
    b1 = 0*ones(1,length(t)); % Adaptation variable of Comparator population 1
    b2 = 0.1*ones(1,length(t)); % Adaptation variable of Comparator population 2
    
    Ipro1 = 0*ones(1,length(t)); % the proprioceptive input variable to comparator population 1
    Ipro2 = 0*ones(1,length(t)); % the proprioceptive input variable to comparator population 2
    I_pro_on = 0.64; % the proprioceptive input current value of the comparator population to set it off
    
    Jd = 7.5; % excitatory input strength from the direction selective population to the comparator population
    Jc = 15; % the inhibitory input strength across the comparator populations
    
    sigma = 0.11; % the amount of noise observed in spontaneous activity 
    T_on = 4; % the on time of the pulse of proprioceptive input
    T_off = 20 - T_on; % the off time of the pulse of proprioceptive input
    rng('shuffle'); % pick a random seed
    
    ds1 = 0; % a temporary variable represneting v1 
    ds2 = 0; % a temporary variable represneting v12
    
    temp1 = sigma*randn; % noise input to comparator population 1
    temp2 = sigma*randn; % noise input to comparator population 2

    for i = 1:length(t)-1 % run the simulation
          IDS1 = 1.15; % the baseline input current to direction selective population 1
          IDS2 = 1.15; % the baseline input current to direction selective population 2

          Ipro1(i) = -2; % the baseline propriocetive input current to comparator selective population 1
          Ipro2(i) = -2; % the baseline propriocetive input current to comparator selective population 1

          % inserting noise into pulses of current only in the first time step
          % and then maintaing this noise throughout
          if mod(i*dt,T_on+T_off) < dt/2 
              temp1 = sigma*randn;
              temp2 = sigma*randn;
          elseif mod(i*dt, T_on+T_off) == T_on
              temp1 = 0;
              temp2 = 0;
          end

        % inserting current into the Direction selective population 1 for a duration of T_CS
        if (i*dt < 500) || (i*dt > 500 + T_CS)
         
        % and adding proprioceptive input pulses       
            if mod(i*dt, T_on + T_off) < T_on
                Ipro1(i) = I_pro_on + temp1;
                Ipro2(i) = I_pro_on + temp2;
            end
        else
            IDS1(i) = 3.1;
            if mod(i*dt, T_on + T_off) < T_on
                Ipro1(i) = I_pro_on + temp1;
                Ipro2(i) = I_pro_on + temp2;
            end
        end
        
        Ip1 = Ipro1(i);
        Ip2 = Ipro2(i);
        
        % Running the differential equation for the next time step by calling the function HH
         
        dydt = HH(i,[u1(:,i), u2(:,i), a1(:,i), a2(:,i), a3(:,i), a4(:,i), b1(:,i), b2(:,i)]);
     
          
        % create next step for each variable
        
        u1(:,i+1) = u1(:,i) + dydt(1)*dt;
        u2(:,i+1) = u2(:,i) + dydt(2)*dt;
        a1(:,i+1) = a1(:,i) + dydt(3)*dt;
        a2(:,i+1) = a2(:,i) + dydt(4)*dt;
        a3(:,i+1) = a3(:,i) + dydt(5)*dt;
        a4(:,i+1) = a4(:,i) + dydt(6)*dt;
        b1(:,i+1) = b1(:,i) + dydt(7)*dt;
        b2(:,i+1) = b2(:,i) + dydt(8)*dt;

        ds1 = f(IDS1-a1(i+1) - a3(i+1));
        ds2 = f(IDS2-a2(i+1) - a4(i+1));

        v1(i+1) = ds1;
        v2(i+1) = ds2;

    end
    
    % save the variables in a file
    save(strcat('File_to_save'),'u1','u2')
    
end

% This is a nonlinear function to saturate the firing rates so that they don't exceed
% a certain threshold and so that they are poisitive

function outf=f(v)

    v_sat=4.0;
    
    if v > 0
        outf = v;
    else
        outf = 0;
    end
    
    if v > v_sat 
        outf = v_sat; 
    end
    
end

% The function of differential equations

function dydt = HH(~,y)

    global Ip1 Ip2 Jd Jc ds1 ds2 

    u1 = y(:,1);
    u2 = y(:,2);
    a1 = y(:,3);
    a2 = y(:,4);
    a3 = y(:,5);
    a4 = y(:,6);
    b1 = y(:,7);
    b2 = y(:,8);

    PhiF = 6; % Fast adaptation variable strength
    PhiS = 9.4; % Slow adaptation variable strength
    PhiB = 11; % Comparator adaptation variable strength

    TauA1 = 24; % Fast adaptation time scale of direction selection variables
    TauA2 = 294; % Slow adaptation time scale of direction selection variables
    TauB = 10; % Adaptation time scale of comparator variables
              
    Tau = 200*(10^(-3)); % Integration time scale of population activity of comparator populations

    dydt = [(1/Tau)*(-u1+f(Jd*ds1-Jc*u2+Ip1-b1)),(1/Tau)*(-u2+f(Jd*ds2-Jc*u1+Ip2-b2)),...
        (1/TauA1)*(-a1+PhiF*ds1),(1/TauA1)*(-a2+PhiF*ds2),...
        (1/TauA2)*(-a3+PhiS*ds1),(1/TauA2)*(-a4+PhiS*ds2),...
        (1/TauB)*(-b1+PhiB*u1),(1/TauB)*(-b2+PhiB*u2)];

end
