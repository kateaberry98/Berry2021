clear all
close all
tic;
%  rng(100) % random number generator for now so can compare plots more easily
%  rng(20)
%% Runs 'Sim' simulations of ABM for panic-buying model with recruitment, defection and lemming effect
% M0r = initial density of panic-buyers
% M0b = initial density of bystanders
% Pmr = movement rate for panic-buyers (unbiased in direction)
% Pmb = movement rate for bystanders (unbiased in direction)
% Rec = individual recruitment rates
% Def = individual defection rates
% Ttot = Overall simulation time
% Sim = Number of identically-prepared simulations to run
% rho = lemming effect parameter 0 < rho < 0.75
% PanicShop = poisson mean products panic-buyers take
% NormalShop = "" bystanders take
% ProdPoisRate = poisson mean products on shelves
% PanicW = weighting of panic-buyers vs products
%% Setup and Initial Conditions

X=61; %100 % Number of x-ordinates
Y=25; %100 % Number of y-ordinates

% Section 2.3 Logistic Growth
% Rec=[0 0.2 0.4 0.6 0.8]; % DEFINE INDIVIDUAL RECRUITMENT RATES HERE
% Def=[0 0.1875 0.375 0.5625 0.75]; % DEFINE INDIVIDUAL DEFECTION RATES HERE

% Section 2.3 Allee effect
% Rec=[0 0 0 0.15 0.6];
% Def=[0 0 0 0 0.2];

% Section 3 spatially-dependent figures & Section 4
% Rec=[0 0 0 0.225 0.9];
% Def=[0 0 0.85/6 0.425 0.85];

% Section 3 spatially-uniform figures
Rec=[0 0.225 0.45 0.675 0.9];
Def=[0 0.2125 0.425 0.6375 0.85];

% Section 3 spatially-uniform figures
% Rec=[0 0.2125 0.425 0.6375 0.85];
% Def=[0 0.225 0.45 0.675 0.9];

% Rec=[0 3/20 3/10 9/20 6/10]; 
% Def=[0 1/20 1/10 3/20 2/10];

Ttot=20;
Sim=40;

M=100*max(max(Def),max(Rec));

Pmr=M; % Moving prob/rate for panic-buyers
Pmb=M; % Moving prob/rate for bystanders

% Initial panic-buyer density
% M0r = 0.1; %Without aisles
% M0r = 0.035; %Busier supermarket
M0r = 0.025; %Aaverage supermarket
% M0r = 0.015; %Quieter supermarket

% Initial bystander density 
% M0b=0.2-M0r; %Without aisles (equal to M0r)
% M0b=0.07-M0r; %Busier supermarket (equal to M0r)
M0b=0.05-M0r; %Average supermarket (equal to M0r)
% M0b=0.03-M0r; %Quieter supermarket (equal to M0r)
% M0b = 0.025; % Fixed

supermarket = input('Including aisles? 1 for Yes or 0 for No: ');
% supermarket = 0;
if supermarket == 0
    uniform = input('Spatially-uniform initial conditions?: 1 for Yes or 0 for No: ');
end

rho = input('Give rho between 0 and 3/4: '); % Lemming effect, rho > ~0.7 = products do not empty
PanicShop = 3;
NormalShop = 1;
ProdPoisRate = 2.5;
PanicW = 3;
alpha = [];
% beta = [];

%% Nearest neighbour index structures & thetas

XX=zeros(Y,X); % To contain row index i for each (i,j)
YY=zeros(Y,X); % To contain column index j for each (i,j)

for i=1:Y
    for j=1:X
        XX(i,j)=i;
        YY(i,j)=j;
    end
end
S=(X*Y+1)*ones(X*Y,4); % Set up matrix for neighbours, (XY+1) where less neighbours

for i=1:Y
    for j=1:X
        Z=(XX-XX(i,j)).^2+(YY-YY(i,j)).^2-1; % Metric for distance from (i,j)
        V=find(Z<=1e-1 & Z>-1); % Chooses sites neighbouring (i,j)
        S(i+(j-1)*Y,1:length(V))=V; % Lists neighbour indexes for (i,j) indexed by J, along row J
    end
end

S(S==0)=X*Y+1; % Buffer

T = (X*Y+1)*ones(2,8,X*Y);
Thetas = zeros(1,8);

for i=1:Y
    for j=1:X
        Zinf=max(abs(XX-XX(i,j)),abs(YY-YY(i,j)))-1; % Infinity norm distance from (i,j) 
        V=find(Zinf<=1e-2 & Zinf>-1); % Chooses (max 8) nearest sites surrounding (i,j)
        T(1,1:length(V),i+(j-1)*Y)=V; % Lists nearest neighbour indexes for 
        for k = 1:length(V)
            Thetas(k) = atan(((YY(V(k))-YY(i,j))./(XX(V(k))-XX(i,j))))+pi.*(XX(V(k))<XX(i,j))-pi/2;
        end
        T(2,1:length(V),i+(j-1)*Y)=Thetas(1:length(V)); %neighbour thetas 
   end
end
%disp(T)

% Supermarket Aisles

A = zeros(Y,X);
if supermarket==1
    for j = 6:20
        for i = [6:8,14:16,22:24,30:32,38:40,46:48,54:56]
            A(j,i)=1;
        end
    end
    a = sum(sum(A))/(X*Y); %density of X*Y
end

% for rho = [0,7/90,7/45,7/30,14/45,7/18,7/15,49/90,28/45,0.7]
for s=1:Sim 
  disp(s)

    % Products

    % Singular Value
    % eta = 5;
    % P = eta*A;

    % Randomly distributed
    P = poissrnd(ProdPoisRate,Y,X).*A; % Poisson distribution (between 0 and 9)

    % Shelves visited
    Visit = zeros(X*Y,1);

    % Clearing unreachable products
    Aind = find(A==1);
        for i = 1:length(Aind)
            if sum(A(Moore(Aind(i),X,Y)))==8
                P(Aind(i)) = 0; 
            end
        end

    % Recording initial number of products    
    Pstart = P;
    Pcount(1) = sum(sum(P));


    CCR=zeros(Y,X); % Empty lattice for tracking panic-buyers (0s will denote empty sites)
    CCB=zeros(Y,X); % Empty lattice for tracking bystanders
    
    % Initial conditions
   
        m0r=M0r; % Initial panic-buyer density (for spatial-uniformity)
        m0b=M0b; % Initial bystander density (for spatial-uniformity)
        
    if supermarket==1 % Spatially uniform with aisles
        II=find(A==0); % Finds indices of where panic-buyers are allowed to be potentially placed (avoiding shelves)
        CCR(II(rand(1,length(II))<m0r/(1-a)))=1; % Places the panic-buyers according to density (and where no shelves)
        II=find((A==0) & (CCR==0)); % Finds indices of where bystanders are allowed to be potentially placed (avoiding shelves)
        CCB(II(rand(1,length(II))<m0b/(1-a-m0r)))=1; % Places the bystanders acc. to density (and where no shelves or panic-buyers)
    elseif supermarket==0 && uniform==1 % Spatially uniform without aisles
        CCR = double(rand(Y,X)<m0r); % Places the panic-buyers according to density
        II=find(CCR==0); % Finds indices of where bystanders are allowed to be potentially placed
        CCB(II(rand(1,length(II))<m0b/(1-m0r)))=1; % Places the bystanders acc. to density (and where no panic-buyers)
    else % Block Initial Conditions
        CCR(:,11:25)=1;
        CCB(:,36:50)=1;
    end
        
    QRstart=sum(sum(CCR)); % Counts the total number of panic-buyers initially
    QBstart=sum(sum(CCB)); % Counts the total number of bystanders initially
    
    
    CRstart=CCR; % Initial profile of panic-buyers
    CBstart=CCB; % Initial profile of bystanders
    
    
    CCR(Y,X+1)=0;
    CCB(Y,X+1)=0;
    
    %% Running Simulations
    
    % Runs a simulation for each s from same initial conditions

    Time=0; % Time starts at zero
    j=1; % Index to be used for steps
    
    tau=0; % Actually, to be a vector of cumulative timesteps
    QR=QRstart; % Actually, to be a vector of total numbers of panic-buyers
    QB=QBstart; % Actually, to be a vector of total numbers of bystanders
    
    CR=CCR; % To manipulate the lattice
    CB=CCB;
    
    JJ=randi([1,X*Y],1,10000); % Vector (length 10000) of [1,XY] randm integers
    Step = randi([1,4],1,10000);
    U1=rand(1,10000); % First vector (length 10000) of Unif(0,1)s
    U2=rand(1,10000); % Second vector (length 10000) of Unif(0,1)s
    U3=rand(1,10000); % Third vector (length 10000) of Unif(0,1)s
    U4=rand(1,10000); % Fourth vector (length 10000) of Unif(0,1)s
    y=1; % Arbitrary index
    
    
   % Agent snapshots
    figure(1)
    if supermarket==1
        spy(A(:,1:X),'ks')
        hold on
    end
    spy(CR(:,1:X),'r')
    hold on
    spy(CB(:,1:X),'b')
    delete(findall(findall(gcf,'Type','axe'),'Type','text'))
    hold off
    title('Before, Time = 0')
%     pause(0.3)  
%     print -dpng 
    

    %% Giant Gillespie Algorithm
   
     while Time<Ttot % While simulation has time left to run
        J=JJ(y);
        while CB(J)==0 && CR(J)==0
            y=y+1;
            if y==10001
                JJ=randi([1,X*Y],1,10000);
                Step = randi([1,4],1,10000);
                U1=rand(1,10000);
                U2=rand(1,10000);
                U3=rand(1,10000);
                U4=rand(1,10000);
                y=1;
            end
            J=JJ(y);
        end
        % First pick agent, then pick agent event later
        if CR(J)==1
            % Agent chosen will be a panic-buyer
            Agent=1;
            %             Ind = find(CR); % Finds lattice sites occupied by panic-buyers
            %             J=Ind(randi(length(Ind))); % Pick a panic-buyer
            v=S(J,:); % Contains indexes of neighbours of J
            Shopmean = PanicShop;
            
            Nr=sum(CR(v)); % Number of panic-buyer neighbours
            Nb=sum(CB(v)); % Number of bystander neighbours
            Pr=Rec(Nr+1); % Picks the relevant individual-level rec. rate
            Pd=Def(Nb+1); % Picks the relevant individual-level def. rate
        elseif CB(J)==1
            % Agent chosen will be a bystander
            Agent=2;
            %             Ind = find(CB); % Finds lattice sites occupied by bystanders
            %             J=Ind(randi(length(Ind))); % Pick a bystander
            v=S(J,:); % Contains indexes of neighbours of J
            Shopmean = NormalShop;
            
            Nr=sum(CR(v)); % Number of panic-buyer neighbours
            Nb=sum(CB(v)); % Number of bystander neighbours
            Pr=Rec(Nr+1); % Picks the relevent individual-level rec. rate
            Pd=Def(Nb+1); % Picks the relevant individual-level def. rate
        else
            error('NULL')
        end
        
        prop=(Pmr+Pd)*QR(j)+(Pmb+Pr)*QB(j); % Computes (approximate) propensity of the whole system
        tau(j+1)=tau(j)-log(U1(y))/prop; % Next entry of tau is time after timestep
        
        QR(j+1)=QR(j); % Creates a new element for defining later
        QB(j+1)=QB(j); % Creates a new element for defining later

 % Panic-buyers                 
        if Agent==1 && U2(y)<Pmr/(Pmr+Pd)  % If panic-buyer chosen, in the event a panic-buyer moves
            
 % Theta calculations                
            Tj = T(:,:,J);
            THETA = Tj(2,CR(Tj(1,:))==1); % Lists only theta values of neighbouring panic-buyers.
             
    % Panic-buyer influence             
            for i=1:PanicW
                THETA = [THETA,THETA];
            end
            
    % Product influence             
            Check = find(Tj(1,:)==X*Y+1, 1);  % Check if agent is at edge of lattice
            if ~isempty(Check)
                ind = Tj(1,A(Tj(1,1:Check-1))==1); % Don't want to check if these are shelves
            else
                ind = Tj(1,A(Tj(1,:))==1);
            end
                    
            if ~isempty(ind) 
                dex = zeros(1,length(ind));
                for k=1:length(ind)
                    dex(k) = find(Tj(1,:)==ind(k));
                    THETA = [THETA, Tj(2,dex(k))*ones(1,P(ind(k)))];
                end
            end
            Dir = ThetaAvgKate(THETA); % Computes average direction influenced by neighbouring panic-buyers & products
            
 % Movement            
            if Dir == 'N'
                if U3(y) < 1/4 + rho
                    Jnext = J-1;
                elseif U3(y) < 1/2 + 2*rho/3
                    Jnext = J+1;
                elseif U3(y) < 3/4 + rho/3
                    Jnext = J+Y; 
                else
                    Jnext = J-Y;
                end
            elseif Dir == 'S'
                if U3(y) < 1/4 + rho
                    Jnext = J+1;
                elseif U3(y) < 1/2 + 2*rho/3
                    Jnext = J+Y;
                elseif U3(y) < 3/4 + rho/3
                    Jnext = J-Y;
                else
                    Jnext = J-1;
                end
            elseif Dir == 'E'
                if U3(y) < 1/4 + rho
                    Jnext = J+Y;
                elseif U3(y) < 1/2 + 2*rho/3
                    Jnext = J-Y;
                elseif U3(y) < 3/4 + rho/3
                    Jnext = J-1;
                else
                    Jnext = J+1;
                end
            elseif Dir == 'W'
                if U3(y) < 1/4 + rho
                    Jnext = J-Y;
                elseif U3(y) < 1/2 + 2*rho/3
                    Jnext = J-1;
                elseif U3(y) < 3/4 + rho/3
                    Jnext = J+1;
                else
                    Jnext = J+Y;
                end
            else % Dir == 'O'
                Jnext = v(Step(y)); % Unbiased movement
            end
                                  
     
            if  Jnext<(X*Y+1) && Jnext>0 && CR(Jnext)==0 && CB(Jnext)==0 && A(Jnext)==0 % Checks site is in matrix and not occupied
                CR(Jnext)=1; % Neighbouring site now contains panic-buyer
                CR(J)=0; % Site originally selected now empties
                Visit(Jnext,:) = Visit(J,:); % Visited shelves moved to neighbouring site
                Visit(J,:) = 0; % Visited shelves for original site now empties
                J = Jnext;
            else % Movement aborted
            end
            
 % Shop phase
            NP = Moore(J,X,Y); % Neighbours in lattice
            Prod = find(P(NP)>0); % Check if any neighbours are shelves with products
            if ~isempty(Prod) 
                for i = 1:length(Prod)
                    if ~ismember(NP(Prod(i)),Visit(J,:)) % If agent hasn't already visited this shelf
                        P(NP(Prod(i))) = P(NP(Prod(i))) - min(poissrnd(Shopmean),P(NP(Prod(i)))); % Agent may take something
                        visitind = find(Visit(J,:) == 0,1); % find end of this agent's list
                        if (isempty(visitind)) % If no more room in Agent's list
                            Visit = [Visit zeros(X*Y,1)]; % Add a column
                            Visit(J,end) = NP(Prod(i)); % Add to new end of row
                        else
                            Visit(J,visitind) = NP(Prod(i)); % Add to end of row
                        end
                    end
                end
            end
            
 % Defection                             
        elseif Agent==1 % If panic-buyer chosen, in the event a panic-buyer defects
            CR(J)=0; % Site originally selected loses panic-buyer...
            CB(J)=1; % ... who turns into a bystander
            QR(j+1)=QR(j+1)-1; % Lose panic-buyer
            QB(j+1)=QB(j+1)+1; % Gain bystander

 % Bystanders                      
        elseif Agent==2 && U2(y)<Pmb/(Pmb+Pr) % If bystander chosen, in the event a bystander moves
            
 % Theta calculations            
            Tj = T(:,:,J);
            THETA = Tj(2,CR(Tj(1,:))==1); % Lists only theta values of neighbouring panic-buyers.
            
    %  Panic-buyer influence     
            for i=1:PanicW 
                THETA = [THETA,THETA];
            end
            
            if ~isempty(THETA)
                THETA = pi - THETA; % Because bystanders are repelled by panic-buyers
            end
            
    % Product influence                        
            Check = find(Tj(1,:)==X*Y+1, 1); % Check if agent is at edge of lattice
            if ~isempty(Check)
                ind = Tj(1,A(Tj(1,1:Check-1))==1); % Don't want to check if these are shelves
            else
                ind = Tj(1,A(Tj(1,:))==1);
            end
            
            if ~isempty(ind)      
                dex = zeros(1,length(ind));
                for k=1:length(ind)
                    dex(k) = find(Tj(1,:)==ind(k));
                    THETA = [THETA, Tj(2,dex(k))*ones(1,P(ind(k)))];
                end
            end
            Dir = ThetaAvgKate(THETA); % Computes average direction influenced by neighbouring panic-buyers & products
            
 % Movement                        
            if Dir == 'N'
                if U3(y) < 1/4 + rho
                    Jnext = J-1;
                elseif U3(y) < 1/2 + 2*rho/3
                    Jnext = J+1;
                elseif U3(y) < 3/4 + rho/3
                    Jnext = J+Y;
                else
                    Jnext = J-Y;
                end
            elseif Dir == 'S'
                if U3(y) < 1/4 + rho
                    Jnext = J+1;
                elseif U3(y) < 1/2 + 2*rho/3
                    Jnext = J+Y;
                elseif U3(y) < 3/4 + rho/3
                    Jnext = J-Y;
                else
                    Jnext = J-1;
                end
            elseif Dir == 'E'
                if U3(y) < 1/4 + rho
                    Jnext = J+Y;
                elseif U3(y) < 1/2 + 2*rho/3
                    Jnext = J-Y;
                elseif U3(y) < 3/4 + rho/3
                    Jnext = J-1;
                else
                    Jnext = J+1;
                end
            elseif Dir == 'W'
                if U3(y) < 1/4 + rho
                    Jnext = J-Y;
                elseif U3(y) < 1/2 + 2*rho/3
                    Jnext = J-1;
                elseif U3(y) < 3/4 + rho/3
                    Jnext = J+1;
                else
                    Jnext = J+Y;
                end
            else % Dir == 'O'
                Jnext = v(Step(y)); % Unbiased movement
            end
         
            if  Jnext<(X*Y+1) && Jnext>0 && CR(Jnext)==0 && CB(Jnext)==0 && A(Jnext)==0 % Checks site is in matrix and not occupied
                CB(Jnext)=1; % Neighbouring site now contains bystander
                CB(J)=0; % Site originally selected now empties
                Visit(Jnext,:) = Visit(J,:); % Visited shelves moved to neighbouring site
                Visit(J,:) = 0; % Visited shelves for original site now empties
                J = Jnext;
            else % Movement aborted
            end
            
 % Shop phase
            NP = Moore(J,X,Y); % Neighbours in lattice
            Prod = find(P(NP)>0); % Check if any neighbours are shelves with products
            if ~isempty(Prod) 
                for i = 1:length(Prod)
                    if ~ismember(NP(Prod(i)),Visit(J,:)) % If agent hasn't already visited this shelf
                        P(NP(Prod(i))) = P(NP(Prod(i))) - min(poissrnd(Shopmean),P(NP(Prod(i)))); % Agent may take something
                        visitind = find(Visit(J,:) == 0,1); % find end of this agent's list
                        if (isempty(visitind)) % If no more room in Agent's list
                            Visit = [Visit zeros(X*Y,1)]; % Add a column
                            Visit(J,end) = NP(Prod(i)); % Add to new end of row
                        else
                            Visit(J,visitind) = NP(Prod(i)); % Add to end of row
                        end
                    end
                end
            end

 % Recruitment
        elseif Agent==2 % If bystander chosen, in the event a bystander recruits
            CB(J)=0; % Site originally selected loses bystander...
            CR(J)=1; % ... who turns into a panic-buyer
            QR(j+1)=QR(j+1)+1; % Gain panic-buyer
            QB(j+1)=QB(j+1)-1; % Lose bystander          
        end
         
        Time=tau(j+1); % Updating time for the next iteration
        QRend=QR(j+1); % Updating number of panic-buyers after step
        QBend=QB(j+1); % Updating number of bystanders after step
        
        y=y+1;
        if y==10001
            U1=rand(1,10000);
            U2=rand(1,10000);
            U3=rand(1,10000);
            
            U4=rand(1,10000);
            Step = randi([1,4],1,10000);
            JJ=randi([1,X*Y],1,10000);
            y=1;
        end     

     Pcount(j+1) = sum(sum(P)); % Product count  
     
     j=j+1; % Increase step counter to proceed to next step
     
     end
    
    % Agent snapshots
    figure(10)
    if supermarket==1
        spy(A(:,1:X),'ks')
        hold on
    end
    spy(CR(:,1:X),'r')
    hold on
    spy(CB(:,1:X),'b')
    delete(findall(findall(gcf,'Type','axe'),'Type','text'))
    hold off
    title(['After, Time = ', num2str(Time)])
%     pause(0.3)
%     print -dpng 
    
    % Product count plot for single simulations
    if Sim == 1
        figure(11)
        plot(tau,Pcount/Pcount(1),'LineWidth',2)
        xlabel('Time','Interpreter','Latex')
        ylabel('$\hat{P} = P/P_0$','Interpreter','Latex')
        xlim([0 max(tau)])
%         print -dpng 
    end


%% Interpolating and averaging profiles for multiple simulations
    if Sim > 1
    %interpolated output: Phat(tau)
    t=linspace(0,Ttot,1000);
    [tau,U]=unique(tau,'first'); %eliminates any duplicate times
    Pcount=Pcount(U);
    G=griddedInterpolant(tau,Pcount/Pcount(1),'nearest'); %linear interpolation function (Interpolate over one simulation)
    Phat(s,:)=G(t); % Store interpolation vector for s-th simulation
    end
end

%% Plotting average profiles and calculating line of best fit
if Sim > 1
t=linspace(0,Ttot,1000);
MeanPhat=mean(Phat,1); % Average over simulations
StdPhat=std(Phat,1); % Average over simulations
figure(50)
plot(t,MeanPhat,'LineWidth',2) % Plot averaged total 
hold on
xlabel('Time','Interpreter','Latex')
ylabel('Average $\hat{P} = P/P_0$','Interpreter','Latex')
title('Proportion of products $\hat{P}$ vs Time','Interpreter','Latex')
ylim([0 1])
xlim([0 Ttot])

figure(51)
semilogy(t,MeanPhat,'LineWidth',2)
xlabel('Time','Interpreter','Latex')
ylabel('Average $\hat{P} = P/P_0$','Interpreter','Latex')
title('Semilogy plot of $\hat{P}$ vs Time ','Interpreter','Latex')
hold on

% Finds index where MeanPhat is zero, or minimum.
[~,minindex] = min(abs(MeanPhat-0.01));
% Get coefficients of a line fit through the data (y = m*t + c).
coeffs = polyfit(t(1:minindex),log(MeanPhat(1:minindex)),1);
m = coeffs(1); % gradient
b = exp(coeffs(2)); % intercept, c = coeffs(2).
if isempty(alpha)
    alpha = -m;
else
    alpha = [alpha,-m];
end
% if isempty(beta)
%     beta = b;
% else
%     beta = [beta,b];
% end

% figure(51)
semilogy(t(1:minindex), b*exp(m.*t(1:minindex)),'LineWidth',2) % Line of best fit
legend({'$\hat{P}$','Best fit $\hat{P}>0.01$'},'Interpreter','Latex')

% Adding exp(mt)*b to Phat plot
figure(50)
hold on
plot(t,exp(m*t)*b,'LineWidth',2)
legend({'$\hat{P}$','$e^{mt+c}$'},'Interpreter','Latex')

figure(52)
loglog(t,MeanPhat,'LineWidth',2)
xlabel('Time','Interpreter','Latex')
ylabel('Average $\hat{P} = P/P_0$','Interpreter','Latex')
title('loglog plot')

end
% end

%% Alpha Plots - if using for loop in line 145
% edit according to parameter being varied
% rho = [0,7/90,7/45,7/30,14/45,7/18,7/15,49/90,28/45,0.7]; 

% figure(53)
 
% plot(rho,alpha,'b-*','LineWidth',2)
% ylabel('$\alpha$','Interpreter','Latex')
% xlabel('$\rho$','Interpreter','Latex')
% print -dpng 

% plot(PanicShop,alpha,'b-*','LineWidth',2)
% ylabel('$\alpha$','Interpreter','Latex')
% xlabel('Poisson mean number of products taken by panic-buyers','Interpreter','Latex')
% print -dpng panicshopvalpha1

% plot(ProdPoisRate,alpha,'b-*','LineWidth',2)
% ylabel('$\alpha$','Interpreter','Latex')
% xlabel('Poisson mean number of products','Interpreter','Latex')
% print -dpng poisratevalpha3

% plot(PanicW,alpha,'b-*','LineWidth',2)
% ylabel('$\alpha$','Interpreter','Latex')
% xlabel('PanicW','Interpreter','Latex')
% print -dpng panicwvalpha3
% toc;

% plot(M0b,alpha,'b-*','LineWidth',2)
% ylabel('$\alpha$','Interpreter','Latex')
% xlabel('$b_0$','Interpreter','Latex')
% print -dpng b0valpha1

% plot(M0r,alpha,'b-*','LineWidth',2)
% ylabel('$\alpha$','Interpreter','Latex')
% xlabel('$r_0$','Interpreter','Latex')
% print -dpng r0valpha1

beep on; beep;