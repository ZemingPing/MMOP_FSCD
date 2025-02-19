function f  = DSS_of_replace_decision_chromosome_kmeans(intermediate_chromosome, M, V,pop)

[ C,IA,IC ] = unique(intermediate_chromosome(:,1:V),'rows');
intermediate_chromosome = intermediate_chromosome(IA,:);


[N, m] = size(intermediate_chromosome);

% Get the index for the population sort based on the rank
[temp,index] = sort(intermediate_chromosome(:,M + V + 1));

clear temp m

% Now sort the individuals based on the index
for i = 1 : N
    sorted_chromosome(i,:) = intermediate_chromosome(index(i),:);
end

% Find the maximum rank in the current population
max_rank = max(intermediate_chromosome(:,M + V + 1));

% Start adding each front based on rank and crowing distance until the
% whole population is filled.
previous_index = 0;
for i = 1 : max_rank
    % Get the index for current rank i.e the last the last element in the
    % sorted_chromosome with rank i.
    current_index = max(find(sorted_chromosome(:,M + V + 1) == i));
    
    % Check to see if the population is filled if all the individuals with
    % rank i is added to the population.
    if current_index > pop

        remaining = pop - previous_index;
        % Get information about the individuals in the current rank i.
        temp_pop = ...
            sorted_chromosome(previous_index + 1 : current_index, :);

        temp_remain_pop = DSS(temp_pop,remaining,M,V);

        for j = 1 : remaining
            f(previous_index + j,:) = temp_remain_pop(j,:);
        end
        return;
    elseif current_index < pop
        % Add all the individuals with rank i into the population.
        f(previous_index + 1 : current_index, :) = ...
            sorted_chromosome(previous_index + 1 : current_index, :);
    else
        % Add all the individuals with rank i into the population.
        f(previous_index + 1 : current_index, :) = ...
            sorted_chromosome(previous_index + 1 : current_index, :);
        return;
    end
    % Get the index for the last added individual.
    previous_index = current_index;
end
end


