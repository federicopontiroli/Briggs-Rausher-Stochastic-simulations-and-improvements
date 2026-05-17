function [T,Dynamics] = simRSSA(vMinus,vPlus,c,initialState,delta,tMax,dT)
    % [T,Dynamics] = simRSSA(vMinus,vPlus,c,initialState,delta,tMax,dT)
    % Simulation by a simple implementation of RSSA with some optimizations:
    % 1) optimization of the computation of propensity values (to limit the usage of nchoosek)
    % 2) pre-generation of random numbers
    % 3) to limit memory allocation, final timeseries are returned with discretization step dT
    % ---
    % vMinus : stoichiometric matrix of reactants
    % vPlus : stoichiometric matrix of products
    % c: array of stochastic reaction rates
    % initialState: abundances of the initial state
    % delta: fluctuation rate (0-1 range)
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
    
    % variable used to limit the number of allowed simulated steps
    maxAllowedSteps = tMax*500000; % we allow at most an average of 500,000 reaction events per unit of time
    
    % pre-generation of some random numbers
    initialLenght = 1000;
    randV = rand(1,initialLenght);
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
    currentStateUp = initialState + round(delta*initialState); % state upper bound
    currentStateDown = initialState - round(delta*initialState); % state lower bound
    
    % check for negative state lower bounds (in case we shift the
    % flucutation interval to have the lower bound set to zero)
    currentStateUp(currentStateDown < 0) = currentStateUp(currentStateDown < 0) + abs(currentStateDown(currentStateDown < 0));
    currentStateDown(currentStateDown < 0) = currentStateDown(currentStateDown < 0) + abs(currentStateDown(currentStateDown < 0));
    
    % computation of reaction propensities (up and down)
    aUp = zeros(size(c));
    aDown = zeros(size(c));
    for j = 1:length(c)
        aUp(j) = computeReactionPropensity(vMinus,c,currentStateUp,j);
        aDown(j) = computeReactionPropensity(vMinus,c,currentStateDown,j);
    end
    
    % sum of the upper bound of reaction propensities
    a0Up = sum(aUp);

    % simulation loop
    nSimulationSteps = 0;
    nFastAccept = 0;
    nSlowAccept = 0;
    nRejections = 0;
    nFluctuationIntUpdates = 0;
    while (currentTime < tMax && nSimulationSteps < maxAllowedSteps)
        stateConsistency = true;
        while (stateConsistency && currentTime < tMax && nSimulationSteps < maxAllowedSteps)
            u = 1;
            reactionAccepted = false;
            if (a0Up > 0) % this if clause should avoid searching for a reaction if no reactions can be selected...
                while (~reactionAccepted)
                    % extraction of three unused random numbers from randV
                    if (usedRandomNumbers + 3 > length(randV))
                        % generation of new random numbers if we reached the end of the array...
                        randV = rand(1,initialLenght);
                        usedRandomNumbers = 0;
                        nRandVResets = nRandVResets + 1;
                    end
                    r1 = randV(usedRandomNumbers+1);
                    r2 = randV(usedRandomNumbers+2);
                    r3 = randV(usedRandomNumbers+3);
                    usedRandomNumbers = usedRandomNumbers + 3;

                    % selection of a reaction candidate by using the upper bound of reaction propensities
                    mu = 1;
                    while sum(aUp(1:mu)) < r1*a0Up
                        mu = mu + 1;
                    end

                    % rejection-based strategy to accept/reject the reaction candidate
                    if (r2 <= aDown(mu)/aUp(mu))
                        reactionAccepted = true; % fast acceptance
                        nFastAccept = nFastAccept + 1;
                    else
                        % computation of the "real" propensity
                        a = computeReactionPropensity(vMinus,c,currentState,mu);
                        if (r2 <= a/aUp(mu))
                            reactionAccepted = true; % slow acceptance
                            nSlowAccept = nSlowAccept + 1;
                        else
                            nRejections = nRejections + 1; % if this is executed, the candidate reaction is rejected
                        end
                    end
                    u = u*r3; % in any case I update u to update the final computation of tau
                end

                % computation of the next tau
                tau = (1/a0Up)*log(1/u);

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
                % if nothing can be done, we just go ahead unitl tMax
                % without further updating the current state
                % PS: I save two times the state to be sure of capturing
                % the steady state condition of the model
                i = i+1;
                T(i) = currentTime;
                Dynamics(i,:) = currentState;
                currentTime = tMax;
                i = i+1;
                T(i) = currentTime;
                Dynamics(i,:) = currentState;
            end
            

            % update of the number of simulation steps
            nSimulationSteps = nSimulationSteps + 1;
            
            % check if the updated state is still consistent with its fluctuation interval
            stateConsistency = checkStateConsistency(currentState,currentStateUp,currentStateDown);
        end
        
        % update of the fluctation interval of the model state
        currentStateUp = currentState + round(delta*currentState); % state upper bound
        currentStateDown = currentState - round(delta*currentState); % state lower bound

        % check for negative state lower bounds (in case we shift the
        % flucutation interval to have the lower bound set to zero)
        currentStateUp(currentStateDown < 0) = currentStateUp(currentStateDown < 0) + abs(currentStateDown(currentStateDown < 0));
        currentStateDown(currentStateDown < 0) = currentStateDown(currentStateDown < 0) + abs(currentStateDown(currentStateDown < 0));
        
        % update of reaction propensities (up and down)
        aUp = zeros(size(c));
        aDown = zeros(size(c));
        for j = 1:length(c)
            aUp(j) = computeReactionPropensity(vMinus,c,currentStateUp,j);
            aDown(j) = computeReactionPropensity(vMinus,c,currentStateDown,j);
        end
        
        % update of the sum of the upper bound of reaction propensities
        a0Up = sum(aUp);
        
        nFluctuationIntUpdates = nFluctuationIntUpdates + 1;
    end
    
    % cut of the residual part of timeseries that remained set to NaN
    T = T(~isnan(Dynamics(:,1)));
    Dynamics = Dynamics(~isnan(Dynamics(:,1)),:);
    
    if (nSimulationSteps >= maxAllowedSteps)
        % printing of a warning message to tell to the user that the simulation has been stopped in advance
        disp(' '); % to print an empty line
        disp(['WARNING: the simulation reached the maximum allowed number of simulation steps (' num2str(maxAllowedSteps) ')!']);
        disp(' '); % to print an empty line
    end
    
    disp(['Total number of computed simulation steps: ' num2str(nSimulationSteps)]);
    disp(['Total number of used random numbers: ' num2str((nRandVResets-1)*initialLenght+usedRandomNumbers)]);
    disp(['Number of fluctuation interval updates: ' num2str(nFluctuationIntUpdates) ' (' num2str(nFluctuationIntUpdates/nSimulationSteps*100) '%)']);
    disp(['Number of fast acceptance steps: ' num2str(nFastAccept) ' (' num2str(nFastAccept/(nFastAccept+nSlowAccept+nRejections)*100) '%)']);
    disp(['Number of slow acceptance steps: ' num2str(nSlowAccept) ' (' num2str(nSlowAccept/(nFastAccept+nSlowAccept+nRejections)*100) '%)']);
    disp(['Number of rejection steps: ' num2str(nRejections) ' (' num2str(nRejections/(nFastAccept+nSlowAccept+nRejections)*100) '%)']);
    toc % this prints the simulation runtime (time elapsed from tic to toc) 
end

function a = computeReactionPropensity(vMinus, c, state, reactionIndex)

    a = c(reactionIndex);

    if sum(vMinus(reactionIndex,:) > 0)
        for i = 1:length(state)

            if vMinus(reactionIndex,i) == 1
                a = a * state(i);

            elseif vMinus(reactionIndex,i) == 2
                if state(i) >= 2
                    a = a * (state(i) * (state(i) - 1) / 2);
                else
                    a = 0;
                    return
                end

            elseif vMinus(reactionIndex,i) > 2
                if state(i) >= vMinus(reactionIndex,i)
                    a = a * nchoosek(state(i), vMinus(reactionIndex,i));
                else
                    a = 0;
                    return
                end
            end

        end
    end
end


function test = checkStateConsistency(currentState,currentStateUp,currentStateDown)
    test = true;
    for i = 1:length(currentState)
        test = test && currentState(i) >= currentStateDown(i) && currentState(i) <= currentStateUp(i);
    end
end