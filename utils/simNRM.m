function [T,Dynamics] = simNRM(vMinus,vPlus,c,initialState,tMax,dT)
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
    control=0;
        
    % computation of reaction propensities and tau
    a = zeros(size(c));
    tau = zeros(size(c));
    for j = 1:length(c)
        a(j) = computeReactionPropensity(vMinus,c,currentState,j);
        tau(j) = 1/a(j)*log(1/randV(j)); % in NRM tau is the absolute simulation time, we use the same formula of the FRM because we start from time 0
    end
    usedRandomNumbers = usedRandomNumbers + length(c);
    a0 = sum(a);

    % simulation loop
    nSimulationSteps = 0;
    while (currentTime < tMax && nSimulationSteps < maxAllowedSteps)
        % selection of the first reaction to fire
        [~,mu] = min(tau); 
        if isinf(tau(mu))
            disp(['Infinite firing time reached. Stopping simulation. Current time = ' num2str(currentTime) ', reaction to be fired = ' num2str(mu)]);
            disp(['Current tau values: ' num2str(tau')]); % Print the current tau values
            break;
        end
        
        if (tau(mu) <= tMax) % in NRM tau is the absolute simulation time
            % dynamics update
            currentTime = tau(mu); % in NRM tau is the absolute simulation time
            currentState = currentState + v(mu,:); % I apply reaction mu by means of its row of the stoichiometric matrix
    
            % saving of the current state to the dynamics timeseries if needed
            if (currentTime >= T(i)+dT || currentTime == tMax)
                i = i+1;
                T(i) = currentTime;
                Dynamics(i,:) = currentState;
            end
        else
            % in this case we reach the end of the simulation before firing the reaction
            % PS: I save two times the state to be sure of capturing
            % the steady state condition of the model
            i = i+1;
            T(i) = currentTime;
            Dynamics(i,:) = currentState; % this is the previous state (because the state has been not updated)
            i = i+1;
            currentTime = tMax;
            T(i) = currentTime;
            Dynamics(i,:) = Dynamics(i-1,:);
        end

        % update of reaction propensities and tau (in this simple implementation we don't use the dependency graph)
        for j = 1:length(a)
            if (j == mu)
                % in this case we compute a new tau
                a_old=a(mu);
                a(j) = computeReactionPropensity(vMinus,c,currentState,j);
                if a_old ~= 0 && a(j) == 0
                    disp(['WARNING: Propensity ' num2str(j) ' becomes zero after firing of reaction ' num2str(mu) ', at step ' num2str(i) ', current time ' num2str(currentTime)]);
                end
                if (usedRandomNumbers + 1 > length(randV))
                    % generation of new random numbers if we reached the end of the array...
                    randV = rand(1,initialLenght);
                    usedRandomNumbers = 0;
                    nRandVResets = nRandVResets + 1;
                end
                
                tau(j) = currentTime + 1/a(j)*log(1/randV(usedRandomNumbers+1)); % in NRM tau is the absolute simulation time, we use the same formula of the FRM because we start from time 0
                usedRandomNumbers = usedRandomNumbers + 1; % I use one random number
                
                if isinf(tau(mu))
                    disp(['WARNING: after firing reaction ' num2str(mu) ' at step ' num2str(i) ', at current time ' num2str(currentTime) ', the next tentative time tau of that reaction is infinite.']);
                end
                
                %if isnan(tau(mu))
                %    disp(['WARNING: The tentative firing time for reaction ' num2str(mu) ' is NaN after firing reaction ' num2str(mu) ', at step ' num2str(i-1) ', current time ' num2str(currentTime) '.']);
                %    control = 1;
                %    break;
                %end
            else
                % in this case I scale the tau without using new random numbers
                aNew = computeReactionPropensity(vMinus,c,currentState,j);
                
                if a(j) ~= 0 && aNew == 0
                    disp(['WARNING: Propensity ' num2str(j) ' becomes zero after firing of reaction ' num2str(mu) ', at step ' num2str(i) ', current time ' num2str(currentTime)]);
                end

                tau(j) = currentTime + (a(j)/aNew)*(tau(j)-currentTime);
                if isinf(tau(j))
                    disp(['WARNING: After firing of reaction ' num2str(mu) ', at step ' num2str(i) ', current time ' num2str(currentTime) ', the tentative firing time for reaction ' num2str(j) ' is infinite.']);
                end  
                %if isnan(tau(j))
                %    disp(['WARNING: The tentative firing time for reaction ' num2str(j) ' is NaN after firing reaction ' num2str(mu) ', at step ' num2str(i) ', current time ' num2str(currentTime) '.']);
                %    disp(['Old propensitity of that reaction: a(' num2str(j) ') = ' num2str(a(j)) ', aNew = ' num2str(aNew)]);
                %    control=1;
                %    break;
                %end
                a(j) = aNew;              

            end        
        end
        if control==1
            disp('Old Propensities a(j):'); % Print old propensities
            disp(a(j)); % Display old propensities
            disp('New Propensities aNew:'); % Print new propensities
            disp(aNew); % Display new propensities
            disp(['Index of broken loop: ' num2str(j)]); % Display the index where the loop has broken
            break; % Exit the loop
        end

         if all(isinf(tau))
            disp(['Warning: after firing reaction ' num2str(mu) ', at step ' num2str(i) ', at current time = ' num2str(currentTime) ', all the tentative times are infinite.']);
         end


        % update of the number of simulation steps
        nSimulationSteps = nSimulationSteps + 1;
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