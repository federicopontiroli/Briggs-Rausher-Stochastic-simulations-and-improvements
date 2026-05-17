function [T, Dynamics] = simOptimized_disc(vMinus, vPlus, c, initialState, tMax, dT)

    % ==== PREPARAZIONE ====
    [G,~] = dependencygraph(vMinus,vPlus);
    M = size(vMinus,1);

    dependentReactions = cell(M,1);
    for mu = 1:M
        dependentReactions{mu} = [mu find(G(mu,:)~=0)];
    end

    v = vPlus - vMinus;

    % Timeline discreta
    T = (0:dT:tMax)';
    Dynamics = nan(length(T), length(initialState));

    % Stato iniziale
    currentState = initialState;
    currentTime  = 0;
    Dynamics(1,:) = currentState;

    % Random numbers
    randV = rand(1,1000);
    used = 0;

    % Propensioni iniziali
    a = zeros(size(c));
    for j = 1:length(c)
        a(j) = c(j);
        for s = 1:length(initialState)
            if vMinus(j,s)==1
                a(j) = a(j)*currentState(s);
            elseif vMinus(j,s)==2
                a(j) = a(j)*currentState(s)*(currentState(s)-1)/2;
            end
        end
    end
    a0 = sum(a);

    % Indice salvataggio
    saveIdx = 2;

    % Start timer
    tic;
    stepCount = 0; % Initialize step count

    % ==== LOOP SSA ====
    while currentTime < tMax

        % RNG
        if used+2 > length(randV)
            randV = rand(1,1000);
            used = 0;
        end
        r1 = randV(used+1);
        r2 = randV(used+2);
        used = used + 2;

        % Reazione
        mu = find(cumsum(a) >= r1*a0,1);

        % Tempo
        tau = -log(r2)/a0;
        currentTime = currentTime + tau;

        % Evoluzione stato
        currentState = currentState + v(mu,:);

        % Salvataggio discreto
        while saveIdx <= length(T) && currentTime >= T(saveIdx)
            Dynamics(saveIdx,:) = currentState;
            saveIdx = saveIdx + 1;
        end

        % Update propensioni (solo dipendenti)
        for k = dependentReactions{mu}
            a(k) = c(k);
            for s = 1:length(currentState)
                if vMinus(k,s)==1
                    a(k) = a(k)*currentState(s);
                elseif vMinus(k,s)==2
                    a(k) = a(k)*currentState(s)*(currentState(s)-1)/2;
                 elseif vMinus(k,s) > 2
                    a(k) = a(k) * nchoosek(Dynamics(k,s), vMinus(k,s));
                end
            end
        end
        a0 = sum(a);
        
        % Increment step count
        stepCount = stepCount + 1;
    end

    % Display results
    elapsedTime = toc; % Stop timer
    fprintf('Total steps: %d\n', stepCount);
    fprintf('Total random numbers generated: %d\n', used);
    fprintf('Elapsed time: %.4f seconds\n', elapsedTime); % Display elapsed time

end