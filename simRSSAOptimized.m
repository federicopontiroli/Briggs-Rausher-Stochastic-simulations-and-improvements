
function [T,Dynamics] = simRSSAOptimized(vMinus,vPlus,c,initialState,delta,tMax,dT)
    tic
    %% === DEPENDENCY GRAPH ===
    [G, ~] = dependencygraph(vMinus, vPlus);
    M = size(vMinus,1);
    dependentReactions = cell(M,1);
    for mu = 1:M
        dependentReactions{mu} = [mu find(G(mu,:) ~= 0)];
    end

    %% Initialization
    T = (0:dT:tMax)';
    Dynamics = nan(length(T),length(initialState));
    maxAllowedSteps = tMax*500000;
    
    randV = rand(1,1000);
    usedRandomNumbers = 0;
    nRandVResets = 1;
    
    v = vPlus - vMinus;
    i = 1;
    T(i) = 0;
    Dynamics(i,:) = initialState;
    
    currentTime  = 0;
    currentState = initialState;
    
    % Setup iniziale degli intervalli
    currentStateUp   = currentState + round(delta*currentState);
    currentStateDown = max(currentState - round(delta*currentState),0);

    %% Initial propensities
    M = length(c);
    aUp   = zeros(M,1);
    aDown = zeros(M,1);
    for j = 1:M
        aUp(j)   = computeReactionPropensity(vMinus,c,currentStateUp,j);
        aDown(j) = computeReactionPropensity(vMinus,c,currentStateDown,j);
    end
    a0Up = sum(aUp);

    %% Statistics
    nSimulationSteps = 0;
    nFastAccept = 0;
    nSlowAccept = 0;
    nRejections = 0;
    nFluctuationIntUpdates = 0; % NUOVO: Conta quante volte aggiorniamo gli intervalli

    %% ================= MAIN LOOP =================
    while currentTime < tMax && nSimulationSteps < maxAllowedSteps
        
        % Nota: Rimosso il while(stateConsistency) esterno perché nel Lazy
        % gestiamo la consistenza puntualmente dentro il loop.
        
        u = 1;
        accepted = false;
        while ~accepted
            if usedRandomNumbers + 3 > length(randV)
                randV = rand(1,1000);
                usedRandomNumbers = 0;
                nRandVResets = nRandVResets + 1;
            end
            r1 = randV(usedRandomNumbers+1);
            r2 = randV(usedRandomNumbers+2);
            r3 = randV(usedRandomNumbers+3);
            usedRandomNumbers = usedRandomNumbers + 3;
            
            % Selezione reazione candidata
            mu = find(cumsum(aUp) >= r1*a0Up,1);
            
            % Se aUp o a0Up sono diventati "stale" (vecchi), potrebbero generare
            % indici non validi se non gestiti, ma matematicamente aUp è sempre
            % un upper bound valido finché lo stato è nell'intervallo.
            if isempty(mu) 
                 mu = M; % Fallback di sicurezza numerica
            end

            if r2 <= aDown(mu)/aUp(mu)
                accepted = true;
                nFastAccept = nFastAccept + 1;
            else
                a = computeReactionPropensity(vMinus,c,currentState,mu);
                if r2 <= a/aUp(mu)
                    accepted = true;
                    nSlowAccept = nSlowAccept + 1;
                else
                    nRejections = nRejections + 1;
                end
            end
            u = u*r3;
        end
        
        %% Advance time and state
        tau = log(1/u)/a0Up;
        currentTime = currentTime + tau;
        currentState = currentState + v(mu,:);
        
        if currentTime >= T(i)+dT
            i = i+1;
            T(i) = currentTime;
            Dynamics(i,:) = currentState;
        end
        
        %% === LAZY UPDATE LOGIC ===
        % 1. Identifichiamo le specie toccate dalla reazione
        affectedSpecies = find(v(mu,:) ~= 0);
        
        % 2. Controlliamo se una qualsiasi di queste è uscita dall'intervallo
        needsUpdate = false;
        for s = affectedSpecies
            if currentState(s) < currentStateDown(s) || currentState(s) > currentStateUp(s)
                needsUpdate = true;
                break; % Basta che una sia fuori per forzare l'aggiornamento
            end
        end
        
        % 3. Se (e solo se) serve aggiornare, ricalcoliamo bounds e propensities
        if needsUpdate
            nFluctuationIntUpdates = nFluctuationIntUpdates + 1;
            
            % A. Aggiorniamo SOLO gli intervalli delle specie toccate
            for s = affectedSpecies
                deltaS = round(delta * currentState(s));
                currentStateDown(s) = max(currentState(s) - deltaS, 0);
                currentStateUp(s)   = currentState(s) + deltaS;
            end
            
            % B. Aggiorniamo le propensities SOLO per le reazioni dipendenti
            % (Poiché il grafo delle dipendenze mappa Reazione->Reazione, e le 
            % reazioni in dependentReactions{mu} sono esattamente quelle che usano 
            % le specie modificate da mu, questo è sufficiente).
            for j = dependentReactions{mu}
                aDown(j) = computeReactionPropensity(vMinus,c,currentStateDown,j);
                aUp(j)   = computeReactionPropensity(vMinus,c,currentStateUp,j);
            end
            
            % C. Aggiorniamo la somma totale
            a0Up = sum(aUp);
        end
        
        nSimulationSteps = nSimulationSteps + 1;
    end
    
    %% Cleanup & Display
    T = T(~isnan(Dynamics(:,1)));
    Dynamics = Dynamics(~isnan(Dynamics(:,1)),:);
    
    disp('--- RISULTATI SIMULAZIONE LAZY ---');
    disp(['Total number of computed simulation steps: ' num2str(nSimulationSteps)]);
    disp(['Total number of used random numbers: ' ...
          num2str((nRandVResets-1)*length(randV) + usedRandomNumbers)]);
    disp(['Number of fluctuation interval updates: ' num2str(nFluctuationIntUpdates) ...
          ' (' num2str(nFluctuationIntUpdates/nSimulationSteps*100) '% of steps)']);
    disp(['Number of fast acceptance steps: ' ...
          num2str(nFastAccept) ' (' ...
          num2str(nFastAccept/(nFastAccept+nSlowAccept+nRejections)*100) '%)']);
    disp(['Number of slow acceptance steps: ' ...
          num2str(nSlowAccept) ' (' ...
          num2str(nSlowAccept/(nFastAccept+nSlowAccept+nRejections)*100) '%)']);
    disp(['Number of rejection steps: ' ...
          num2str(nRejections) ' (' ...
          num2str(nRejections/(nFastAccept+nSlowAccept+nRejections)*100) '%)']);
    toc
end

%% --- Helper Functions rimaste invariate ---
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
% Nota: dependencygraph() deve essere presente nel path o nello script