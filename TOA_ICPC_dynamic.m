function [estimate_position, RMSE] = TOA_ICPC_dynamic(LS_gauss, r_measure, targ, anchor_errors, gauss_distance, targ_initial, time, v_node)
    % TOA_ICPC_dynamic: Dynamic algorithm for target position estimation using Time of Arrival (TOA) measurements.
    % This function estimates the target position using a two-step optimization process, considering anchor position errors,
    % measurement noise, and target movement constraints.
    %
    % Inputs:
    %   LS_gauss      : Matrix of noisy anchor positions (N x 2), where N is the number of anchors.
    %   r_measure     : Vector of noisy distance measurements from the target to each anchor (N x 1).
    %   targ          : True target position (1 x 2). Note: This is only used for RMSE calculation, not for estimation.
    %   anchor_errors : Vector of CRLB values for each anchor's position uncertainty (N x 1).
    %   gauss_distance: Variance of the TOA measurement noise.
    %   targ_initial  : Previous estimate of the target position (1 x 2).
    %   time          : Time elapsed since the previous estimate (scalar).
    %   v_node        : Velocity of the target (scalar).
    %
    % Outputs:
    %   estimate_position : Estimated target position (1 x 2).
    %   RMSE             : Root Mean Square Error between the estimated and true target positions.
    %
    % Description:
    % The function uses a two-step optimization process to estimate the target position. In the first step, it estimates
    % the target position and updates the anchor positions using a Particle Swarm Optimization (PSO) algorithm. In the
    % second step, it refines the target position estimate using the updated anchor positions, their corresponding CRLB values,
    % and a movement constraint based on the target's velocity and elapsed time.

    % Get the number of anchors
    num_anchor = size(LS_gauss, 1);

    % Reshape the noisy anchor positions into a single row vector
    LS_gauss_reshape = reshape(LS_gauss.', 1, []);

    % Define lower and upper bounds for anchor position estimates
    Anchor_estimate_lower_bound = LS_gauss_reshape - 50;
    Anchor_estimate_upper_bound = LS_gauss_reshape + 50;

    % Use the iCHAN algorithm to get an initial estimate of the target position
    targ_iCHAN = TOA_iCHAN(LS_gauss, r_measure, 100, 0.5);

    % Define lower and upper bounds for the target position and anchor positions
    lower_bound = [max(targ_iCHAN(1)-50, 0), max(targ_iCHAN(2)-50, 0), Anchor_estimate_lower_bound]; % Lower bounds
    upper_bound = [min(targ_iCHAN(1)+50, 1000), min(targ_iCHAN(2)+50, 1000), Anchor_estimate_upper_bound]; % Upper bounds

    % First optimization step: Update anchor positions and their CRLB values
    options = optimoptions(@particleswarm, 'SwarmSize', 250, 'MaxIterations', 5000, 'FunctionTolerance', 1e-6, 'Display', 'none');
    fitness_function = @(x) calculate_fitness_dsf(x(1:2), x(3:2*(num_anchor+1)), LS_gauss, r_measure, gauss_distance, anchor_errors, num_anchor);
    [x, ~] = particleswarm(fitness_function, 2*(num_anchor+1), lower_bound, upper_bound, options);

    % Update the anchor positions
    Anchor_estimate_update = [x(3:2:end); x(4:2:end)].';

    % Update the CRLB values for the anchor positions
    [~, FIM_inverse] = TOA_local_CRLB(LS_gauss, targ, gauss_distance, anchor_errors);
    CRLB_anchor_update = zeros(1, num_anchor);
    for i = 1:num_anchor
        index_count = 2*i + 1;
        CRLB_anchor_update(i) = trace(FIM_inverse(index_count:index_count+1, index_count:index_count+1));
    end

    % Second optimization step: Refine the target position estimate with movement constraints
    fitness_function = @(x) calculate_fitness_dsf_dynamic(x(1:2), x(3:2*(num_anchor+1)), Anchor_estimate_update, r_measure, gauss_distance, CRLB_anchor_update, num_anchor, targ_initial, time, v_node);
    [x, ~] = particleswarm(fitness_function, 2*(num_anchor+1), lower_bound, upper_bound, options);

    % Calculate the RMSE between the estimated and true target positions
    RMSE = sqrt((x(1) - targ(1))^2 + (x(2) - targ(2))^2);

    % Extract the estimated target position
    estimate_position = [x(1), x(2)];
end

function fitness = calculate_fitness_dsf(position, true_anchor, wrong_anchor, distances, gauss_distance, CRLB_LS1, num_anchor)
    % calculate_fitness_dsf: Fitness function for the Particle Swarm Optimization (PSO) algorithm.
    % This function computes the fitness value based on the difference between estimated and measured distances,
    % as well as the difference between true and noisy anchor positions.
    %
    % Inputs:
    %   position       : Current estimate of the target position (1 x 2).
    %   true_anchor    : True anchor positions (reshaped into a single row vector).
    %   wrong_anchor   : Noisy anchor positions (N x 2).
    %   distances      : Noisy distance measurements (N x 1).
    %   gauss_distance : Variance of the TOA measurement noise.
    %   CRLB_LS1       : CRLB values for the anchor position uncertainties (N x 1).
    %   num_anchor     : Number of anchors.
    %
    % Output:
    %   fitness        : Fitness value to be minimized by the PSO algorithm.

    % Reshape the true anchor positions into a matrix
    true_anchor = reshape(true_anchor, 2, num_anchor)';

    % Initialize arrays for distance and anchor fitness components
    estimated_distance = zeros(num_anchor, 1);
    fitness_distance = zeros(num_anchor, 1);
    fitness_anchor = zeros(num_anchor, 1);

    % Compute fitness components for each anchor
    for i = 1:num_anchor
        wrong_anchor_1 = wrong_anchor(i, :);
        true_anchor_1 = true_anchor(i, :);
        CRLB = CRLB_LS1(i);

        % Compute the estimated distance between the target and the anchor
        estimated_distance(i) = norm(position - true_anchor_1);

        % Compute the distance fitness component
        fitness_distance(i) = (estimated_distance(i) - distances(i))^2 / (gauss_distance^2);

        % Compute the anchor position fitness component
        fitness_anchor(i) = (norm(wrong_anchor_1 - true_anchor_1))^2 / CRLB;
    end

    % Sum the fitness components to get the total fitness value
    fitness = sum(fitness_distance) + sum(fitness_anchor);
end

function fitness = calculate_fitness_dsf_dynamic(position, true_anchor, wrong_anchor, distances, gauss_d, CRLB_LS1, num_anchor, targ_initial, time, v_node)
    % calculate_fitness_dsf_dynamic: Fitness function for the PSO algorithm with movement constraints.
    % This function computes the fitness value based on the difference between estimated and measured distances,
    % the difference between true and noisy anchor positions, and a movement constraint based on the target's velocity.
    %
    % Inputs:
    %   position       : Current estimate of the target position (1 x 2).
    %   true_anchor    : True anchor positions (reshaped into a single row vector).
    %   wrong_anchor   : Noisy anchor positions (N x 2).
    %   distances      : Noisy distance measurements (N x 1).
    %   gauss_d        : Variance of the TOA measurement noise.
    %   CRLB_LS1       : CRLB values for the anchor position uncertainties (N x 1).
    %   num_anchor     : Number of anchors.
    %   targ_initial   : Previous estimate of the target position (1 x 2).
    %   time           : Time elapsed since the previous estimate (scalar).
    %   v_node         : Velocity of the target (scalar).
    %
    % Output:
    %   fitness        : Fitness value to be minimized by the PSO algorithm.

    % Reshape the true anchor positions into a matrix
    true_anchor = reshape(true_anchor, 2, num_anchor)';

    % Initialize arrays for distance and anchor fitness components
    estimated_distance = zeros(num_anchor, 1);
    fitness_distance = zeros(num_anchor, 1);
    fitness_anchor = zeros(num_anchor, 1);

    % Compute fitness components for each anchor
    for i = 1:num_anchor
        wrong_anchor_1 = wrong_anchor(i, :);
        true_anchor_1 = true_anchor(i, :);
        CRLB = CRLB_LS1(i);

        % Compute the estimated distance between the target and the anchor
        estimated_distance(i) = norm(position - true_anchor_1);

        % Compute the distance fitness component
        fitness_distance(i) = (estimated_distance(i) - distances(i))^2 / (gauss_d^2);

        % Compute the anchor position fitness component
        fitness_anchor(i) = (norm(wrong_anchor_1 - true_anchor_1))^2 / CRLB;
    end

    % Compute the movement constraint fitness component
    fitness_move = 0.1 * norm(position - targ_initial) - (time * v_node)^2; % Movement penalty term

    % Sum the fitness components to get the total fitness value
    fitness = sum(fitness_distance) + sum(fitness_anchor) + fitness_move;
end