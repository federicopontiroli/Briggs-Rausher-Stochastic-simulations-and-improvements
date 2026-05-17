function [c, order] = generateStochasticrates(k, vMinus, R)
    % Genera i tassi stocastici a partire dai rate deterministici
    % k: rate deterministici
    % vMinus: matrice stechiometrica dei reagenti
    % R: è tale che V=R/avogadro

    avogadro = 6.02214076e23;
    nReactions = length(k);
    c = zeros(nReactions,1);
    [nRows, ~] = size(vMinus);
    prodFactorials = zeros(nRows,1);
    order = sum(vMinus,2);  % ordine della reazione (somma dei reagenti)
    for i = 1:nRows
        % Prende solo le entrate maggiori di zero
        nonZeroEntries = vMinus(i, vMinus(i,:) > 0);
        % Calcola il prodotto dei fattoriali
        prodFactorials(i) = prod(arrayfun(@factorial, nonZeroEntries));
    end
    for i=1:nRows
        if order(i)==1
            c(i)=k(i)*prodFactorials(i);
        elseif order(i)==2
            c(i) = (k(i)*prodFactorials(i))/R;
        elseif order(i)==3
            c(i) = (k(i) * prodFactorials(i)) / (R^2);
        elseif order(i)==4
            c(i) = (k(i) * prodFactorials(i)) / (R^3);
        else
            c(i) = (k(i) * prodFactorials(i)) / (R^(order(i)-1));
        end
        % Increase precision for the results
        c(i) = round(c(i), 15); % rounding to 10 decimal places for more precision
    end
end

