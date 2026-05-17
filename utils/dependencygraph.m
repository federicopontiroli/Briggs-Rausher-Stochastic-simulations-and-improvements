function [dependencyGraph, affectedSpecies] = dependencygraph(vMinus, vPlus)
    % Define the number of reactions
    M = size(vMinus, 1);

    % Initialize the adjacency matrix for the dependency graph
    dependencyGraph = eye(M);
    
    % Initialize the affected species list
    affectedSpecies = cell(M, 1);
    affectsSpecies = cell(M, 1);

    % Compute the affected species and construct the dependency graph
    for i = 1:M
        affectedSpecies{i} = []; % Initialize each list for species affected by reaction R_i
        for j = 1:size(vMinus, 2)
            if (vPlus(i, j) - vMinus(i, j)) ~= 0
                affectedSpecies{i} = [affectedSpecies{i}, j]; % Add species index j to the list of affected species for reaction R_i
            end
        end
    end
    for i = 1:M
        affectsSpecies{i} = []; % Initialize each list for species affected by reaction R_i
        for j = 1:size(vMinus, 2)
            if vMinus(i, j) ~= 0
                affectsSpecies{i} = [affectsSpecies{i}, j]; % Add species index j to the list of affected species for reaction R_i
            end
        end
    end

    % Construct the dependency graph
    for i = 1:M
        for j = i+1:M
            % Check for common affected species between reaction R_i and R_j
            if ~isempty(intersect(affectedSpecies{i}, affectsSpecies{j}))
                dependencyGraph(i, j) = 1; % Create a link from R_i to R_j     
            end
            if ~isempty(intersect(affectedSpecies{j}, affectsSpecies{i}))
                dependencyGraph(j, i) = 1; % Create a link from R_i to R_j
            end
        end
    end
end
