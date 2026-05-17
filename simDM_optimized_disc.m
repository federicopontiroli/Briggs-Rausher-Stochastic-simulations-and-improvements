function [T,Dynamics] = simDM_optimized_disc(vMinus,vPlus,c,initialState,tMax,dT)
    % [T,Dynamics] = simDM_optimized_disc(vMinus,vPlus,c,initialState,tMax,dT)
    % Simulation by the DM with some optimizations:
    % 1) optimization of the computation of propensity values (to limit the usage of nchoosek)
    % 2) pre-generation of random numbers
    % 3) to limit memory allocation, final timeseries are returned with discretization step dT
    % ---
    % vMinus : stoichiometric matrix of reactants
    % vPlus : stoichiometric matrix of products
    % c: array of stochastic reaction rates
    % initialState: abundances of the initial state
    % tMax: max time instant to simulate
    % dT: discretization step for the returned dynamics
    %
    % For simulating the Oregonator:
    % vMinus = [1 1 0; 0 1 0; 1 0 0; 2 0 0; 0 0 1]
    % vPlus  = [0 0 0; 1 0 0; 2 0 1; 0 0 0; 0 1 0]
    % c = [0.1 2 104 0.016 26]
    % initialState = [500 1000 2100]

    tic % this allows to compute the simulation runtime
    
    % T and Dynamics initialization
    T = (0:dT:tMax)';
    Dynamics = nan(length(T),length(initialState)); % for each step we provide the abundance of each variable
    
    % pre-generation of some random numbers
    initialLenght = 1000;
    randV = rand(1,initialLenght); %so we generate the random numbers before, so that we call this function only one time. This reduces time. 
    nRandVResets = 1; % I keep in memory how many times I generate the vector, so that I can compute how many random numbers have been used
    usedRandomNumbers = 0;
    
    % computation of the stoichiometric matrix
    v = vPlus - vMinus;
    
    % setting initial state
    i = 1;
    T(i) = 0;
    Dynamics(i,:) = initialState;
    currentTime = 0;
    currentState = initialState;
    nSimulationSteps = 0;
    while currentTime < tMax
        % computation of reaction propensities
       % computation of reaction propensities
        a = c;
        for j = 1:length(a)
            if sum(vMinus(j,:) > 0)
                for jj = 1:length(currentState)
        
                    if vMinus(j,jj) == 1
                        a(j) = a(j) * currentState(jj);
        
                    elseif vMinus(j,jj) == 2
                        if currentState(jj) >= 2
                            a(j) = a(j) * (currentState(jj) * (currentState(jj) - 1) / 2);
                        else
                            a(j) = 0;
                            break
                        end
        
                    elseif vMinus(j,jj) > 2
                        if currentState(jj) >= vMinus(j,jj)
                            a(j) = a(j) * nchoosek(currentState(jj), vMinus(j,jj));
                        else
                            a(j) = 0;
                            break
                        end
                    end
        
                end
            end
        end

        
        % sum of reaction propensities
        a0 = sum(a);
    
        % extraction of two unused random numbers from randV
        if (usedRandomNumbers + 2 > length(randV))
            % generation of new random numbers if we reached the end of the array...
            randV = rand(1,initialLenght);
            usedRandomNumbers = 0;
            nRandVResets = nRandVResets + 1;
        end
        r1 = randV(usedRandomNumbers+1);
        r2 = randV(usedRandomNumbers+2);
        usedRandomNumbers = usedRandomNumbers + 2;
        
        % selection of the reaction to fire
        mu = 1;
        while sum(a(1:mu)) < r1*a0
            mu = mu + 1;
        end
        % mu now is the smallest reaction index for which sum(a(1:mu)) >= r(1)*a0
       
        % computation of the next tau
        tau = (1/a0)*log(1/r2);

        if (currentTime + tau <= tMax)
            % dynamics update
            currentTime = currentTime + tau;
            currentState = currentState + v(mu,:); % I apply reaction mu by means of its row of the stoichiometric matrix
    
            % saving of the current state to the dynamics timeseries if needed
            if currentTime >= T(i)+dT
                i = i+1;
                T(i) = currentTime;
                Dynamics(i,:) = currentState;
            end
        else
            % in this case we reach the end of the simulation before firing
            % the reaction
            i = i+1;
            T(i) = currentTime;
            Dynamics(i,:) = currentState; % this is the previous state (because the state has been not updated)
            i = i+1;
            currentTime = tMax;
            T(i) = currentTime;
            Dynamics(i,:) = Dynamics(i-1,:);
        end
        
        % update of the number of simulation steps
        nSimulationSteps = nSimulationSteps + 1;
    end
    
    % cut of the residual part of timeseries that remained set to NaN
    T = T(~isnan(Dynamics(:,1)));
    Dynamics = Dynamics(~isnan(Dynamics(:,1)),:);

    disp(['Total number of computed simulation steps: ' num2str(nSimulationSteps)]);
    disp(['Total number of used random numbers: ' num2str((nRandVResets-1)*initialLenght+usedRandomNumbers)]);

    toc % this prints the simulation runtime (time elapsed from tic to toc) 
end

