% The City College of New York, City University of New York
% Written by Olga Chsherbakova
% Date: October, 2023

% Genetic Algorithm for Traveling Salesman Problem

solve_tsp('cities.csv');
function solve_tsp(filename)
    % Solve the Traveling Salesman Problem using Genetic Algorithm and 
    % visualize the best path

    % Read cities data from CSV file
    cities_table = readtable(filename);
    % take all columns (City Name, Latitude, Longitude)
    cities_data = table2cell(cities_table(:, 1:3)); 
    
    % Generate a weight matrix where the weights are the distances 
    % between cities
    distance_matrix = generate_distance_matrix(cities_data);

    % Parameters for the Genetic Algorithm
    num_generations = 200;
    pop_size = 100;
    mutate_percent = 0.2;
    tourn_size = 5;

    % Run Genetic Algorithm
    best_path = genetic_algorithm_tsp(distance_matrix, num_generations, ...
        pop_size, mutate_percent, tourn_size);

    % Calculate the total distance of the best path
    total_distance = evaluate_individual(best_path, distance_matrix);

    % Print the best path with distances in the command window
    print_best_path(best_path, cities_data,distance_matrix,total_distance);

    % Plot the best path with the total distance in the title
    plot_best_path(best_path, cities_data, total_distance);
end

% GPS function
function distance = get_distance(latitude1,longitude1,latitude2,longitude2)
    R = 6373.0; % radius of the Earth in km
    % Convert from degrees to radians
    lat1 = deg2rad(latitude1);
    lon1 = deg2rad(longitude1);
    lat2 = deg2rad(latitude2);
    lon2 = deg2rad(longitude2);
    % Get the change in coordinates
    dlon = lon2 - lon1;
    dlat = lat2 - lat1;

    % calculate the distance between the two coordinates
    % using Haversine formula
    a = sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2;
    c = 2 * atan2(sqrt(a), sqrt(1 - a));

    distance = R * c;
    % convert distance to miles
    distance = distance * 0.621371;
end

% Function to generate distance matrix for a list of cities
function distance_matrix = generate_distance_matrix(cities)
    num_cities = size(cities, 1);
    distance_matrix = zeros(num_cities, num_cities);
    for i = 1:num_cities
        for j = 1:num_cities
            distance_matrix(i, j) = get_distance(cities{i, 2}, ...
                cities{i, 3}, cities{j, 2}, cities{j, 3});
        end
    end
end

% Function to evaluate the total distance of the path represented by 
% the chromosome
function total_distance = evaluate_individual(chromosome, distance_matrix)
    total_distance = 0;
    for i = 1:length(chromosome) - 1
        total_distance = total_distance + distance_matrix(chromosome(i),...
            chromosome(i + 1));
    end
    total_distance = total_distance + distance_matrix(chromosome(end),...
        chromosome(1));
end

% Function to perform ordered crossover to produce a child
function child = ordered_crossover(parent1, parent2)
    n = length(parent1);
    child = zeros(1, n);

    start_point = randi(n);
    end_point = mod((start_point + randi(n - 2)), n);

    if start_point <= end_point
        child(start_point:end_point) = parent1(start_point:end_point);
    else
        child(start_point:end) = parent1(start_point:end);
        child(1:end_point) = parent1(1:end_point);
    end
    pointer = 1;
    for i = 1:n
        if ~ismember(parent2(i), child)
            while child(pointer) ~= 0
                pointer = pointer + 1;
            end
            child(pointer) = parent2(i);
        end
    end
end

% Function to perform mutation by swapping two genes in the chromosome
function mutated = swap_mutate(chromosome, mutate_percent)
    if rand() < mutate_percent
        n = length(chromosome);
        gene1 = randi(n);
        gene2 = randi(n);
        while gene2 == gene1
            gene2 = randi(n);
        end
        tmp = chromosome(gene1);
        chromosome(gene1) = chromosome(gene2);
        chromosome(gene2) = tmp;
    end
    mutated = chromosome;
end

% Function to perform tournament selection
function selected = tournament_selection(population, fitnesses, tourn_size)
    n = size(population, 1);
    selected = zeros(size(population));
    for i = 1:n
        candidates = randi(n, [1, tourn_size]);
        [~, idx] = min(fitnesses(candidates));
        selected(i, :) = population(candidates(idx), :);
    end
end

% Genetic Algorithm for TSP
function best_path = genetic_algorithm_tsp(distance_matrix,...
        num_generations, pop_size, mutate_percent, tourn_size)
    num_cities = size(distance_matrix, 1);

    % Initialize population
    population = zeros(pop_size, num_cities);
    for i = 1:pop_size
        population(i, :) = randperm(num_cities);
    end
    for gen = 1:num_generations
        fitnesses = zeros(1, pop_size);
        for i = 1:pop_size
            fitnesses(i) = evaluate_individual(population(i, :), ...
                distance_matrix);
        end
        selected = tournament_selection(population, fitnesses, tourn_size);
        children = zeros(size(population));
        for i = 1:2:pop_size
            children(i, :) = ordered_crossover(selected(i, :), ...
                selected(i + 1, :));
            children(i + 1, :) = ordered_crossover(selected(i + 1, :),...
                selected(i, :));
        end
        for i = 1:pop_size
            children(i, :) = swap_mutate(children(i, :), mutate_percent);
        end

        population = children;
    end
    final_fitnesses = zeros(1, pop_size);
    for i = 1:pop_size
        final_fitnesses(i) = evaluate_individual(population(i, :), ...
            distance_matrix);
    end
    [~, best_idx] = min(final_fitnesses);
    best_path = population(best_idx, :);
end

% --- Function to visualize the best path on a map using geoplot
function plot_best_path(best_path, cities_data, total_distance)
    % Extract latitudes and longitudes for the best path
    lats = cell2mat(cities_data(best_path, 2));
    longs = cell2mat(cities_data(best_path, 3));

    % Close the loop by adding the starting city to the end
    lats = [lats; lats(1)];
    longs = [longs; longs(1)];

    figure;
    geoplot(lats, longs, '-o', 'Color', 'red', 'LineWidth', 2);
    geobasemap('landcover');
    title(['Best Route Found for TSP using GA - Total Cost: ' ...
        num2str(total_distance)]);
    
    % Add city names to the map
    city_names = cities_data(best_path, 1);
    for i = 1:length(best_path)
        text(lats(i), longs(i), city_names{i}, 'VerticalAlignment',...
            'bottom', 'HorizontalAlignment', 'center', 'FontSize', 8,...
            'FontWeight', 'bold');
    end
end

% Function to print the best path with distances
function print_best_path(best_path, cities_data, distance_matrix,...
            total_distance)
    num_cities = length(best_path);
    fprintf('BEST PATH COST:\n');
    fprintf('%-15s  %-15s  %-15s  %-15s\n', 'City', 'Latitude', ...
        'Longitude', 'Distance');
    for i = 1:num_cities
        current_city = best_path(i);
        next_city = best_path(mod(i, num_cities) + 1);
        current_city_info = cities_data(current_city, :);
        next_city_info = cities_data(next_city, :);
        distance = distance_matrix(current_city, next_city);
        fprintf('%-15s  %-15.4f  %-15.4f  %-15.2f\n', ...
            current_city_info{1}, current_city_info{2}, ...
            current_city_info{3}, distance);
    end
    fprintf('Total Distance: %f\n', total_distance);
end
