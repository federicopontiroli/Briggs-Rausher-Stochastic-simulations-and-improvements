function [T,Dynamics] = simOptimized(vMinus,vPlus,c,initialState,tMax)
    % SSA ottimizzato con dependency graph
    % Calcolo propensioni manuale per vMinus=1 o vMinus=2
    % vMinus, vPlus: matrici stechiometriche
    % c: vettore dei rate stocastici
    % initialState: stato iniziale
    % tMax: tempo massimo
    

    initialLength = 1000;
    T = nan(initialLength,1); 
    Dynamics = nan(initialLength,length(initialState));
    %propensities = nan(initialLength, size(c, 1)); % Ensure propensities has the correct number of rows based on the size of c
    
    [G,~] = dependencygraph(vMinus,vPlus); % G = dependencyGraph
    M = size(vMinus,1);
    dependentReactionsList = cell(M,1);
    for mu = 1:M
        dependentReactionsList{mu} = [mu, find(G(mu,:) ~= 0)]; % include se stessa
    end
    
    % Pre-generazione random
    randV = rand(1,initialLength);
    nRandVResets = 1; 
    usedRandomNumbers = 0;
    
    v = vPlus - vMinus; % matrice stechiometrica
    
    % Inizializzazione
    i = 1;
    T(i) = 0;
    Dynamics(i,:) = initialState;
    
    % Calcolo iniziale delle propensioni
    a = zeros(size(c));
    for j = 1:length(c)
        a(j) = c(j);
        for jj = 1:length(initialState)
            if vMinus(j,jj) == 1
                a(j) = a(j) * Dynamics(i,jj);
            elseif vMinus(j,jj) == 2
                a(j) = a(j) * Dynamics(i,jj) * (Dynamics(i,jj)-1)/2;
            elseif vMinus(j,jj) > 2
                % Se ci sono coefficienti maggiori di 2, puoi usare nchoosek
                a(j) = a(j) * nchoosek(Dynamics(i,jj), vMinus(j,jj));
            end
        end
    end
    %propensities(1,:) = a(:).';
    a0 = sum(a);
    
    tic; % inizio cronometro
    
    while T(i) < tMax
        % Estrazione due numeri casuali
        if usedRandomNumbers+2 > length(randV)
            randV = rand(1,initialLength);
            usedRandomNumbers = 0;
            nRandVResets = nRandVResets + 1;
        end
        r1 = randV(usedRandomNumbers+1);
        r2 = randV(usedRandomNumbers+2);
        usedRandomNumbers = usedRandomNumbers + 2;

        % Selezione reazione
        mu = 1;
        while sum(a(1:mu)) < r1*a0
            mu = mu + 1;
        end
        
        % Calcolo tau
        tau = (1/a0)*log(1/r2);
        
        % Aggiornamento dinamica
        i = i+1;
        if i > length(T)
            T = [T; nan(initialLength,1)];
            Dynamics = [Dynamics; nan(initialLength,length(initialState))];
        end
        
        if T(i-1) + tau <= tMax
            T(i) = T(i-1) + tau;
            Dynamics(i,:) = Dynamics(i-1,:) + v(mu,:);
        else
            T(i) = tMax;
            Dynamics(i,:) = Dynamics(i-1,:);
            break; % termina il ciclo se superiamo tMax
        end

        for k = dependentReactionsList{mu}
            a(k) = c(k);
            for jj = 1:length(initialState)
                if vMinus(k,jj) == 1
                    a(k) = a(k) * Dynamics(i,jj);
                elseif vMinus(k,jj) == 2
                    a(k) = a(k) * Dynamics(i,jj) * (Dynamics(i,jj)-1)/2;
                elseif vMinus(k,jj) > 2
                    a(k) = a(k) * nchoosek(Dynamics(i,jj), vMinus(k,jj));
                end
            end
        end
        %propensities(i,:) = a(:).';


        a0 = sum(a); % aggiorna propensione totale
    end
    
    % Rimuove le righe non usate
    T = T(1:i);
    Dynamics = Dynamics(1:i,:);
    %propensities=propensities(1:i,:);
    
    % Stop timing and display elapsed time
    elapsedTime = toc;
    disp(['Total number of computed simulation steps: ' num2str(length(T)-1)]);
    disp(['Total number of used random numbers: ' num2str((nRandVResets-1)*initialLength + usedRandomNumbers)]);
    disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);
end
