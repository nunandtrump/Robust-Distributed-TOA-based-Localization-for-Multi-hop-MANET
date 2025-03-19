function OUT_prev = TOA_iCHAN(LS_gauss, r_measure, Sim_num, step)
    % TOA_iCHAN: Time of Arrival (TOA) based localization using the iCHAN algorithm.
    % This function estimates the mobile station (MS) position using iterative least squares.
    %
    % Inputs:
    %   LS_gauss  : Matrix of known sensor positions (m x 2), where m is the number of sensors.
    %   r_measure : Vector of measured distances from the MS to each sensor (m x 1).
    %   Sim_num   : Maximum number of iterations for the algorithm.
    %   step      : Convergence threshold for the RMSE (Root Mean Square Error).
    %
    % Output:
    %   OUT_prev  : Estimated position of the MS (1 x 2).
    %
    % Description:
    % The function iteratively estimates the MS position using the iCHAN algorithm,
    % which is based on weighted least squares. It updates the position estimate
    % until the RMSE between consecutive estimates is below the specified threshold.

    [m, ~] = size(LS_gauss); % Number of sensors
    OUT_prev = zeros(1, 2); % Initialize the previous output (MS position)

    for iter = 1:Sim_num
        % Estimate MS position
        k = sum(LS_gauss.^2, 2); % Sum of squares of sensor positions
        A1 = [-2*LS_gauss, ones(m, 1)]; % Design matrix for the first step
        p1 = r_measure'.^2 - k; % Vector of measured distances squared minus k
        B1 = diag(2*r_measure'); % Weighting matrix for the first step

        if iter > 1
            % Iteratively optimize r_calculate and Q1
            r_calculate = pdist2(OUT_prev, LS_gauss); % Calculate distances using current estimate
            Q1 = diag(ones(1, m) .* (r_measure - r_calculate).^2); % Update weighting matrix
        else
            % For the first iteration, set Q1 as the identity matrix
            Q1 = eye(m);
        end

        % Compute covariance matrix and estimate theta1
        cov1 = B1 * Q1 * B1;
        theta1 = pinv(A1' * pinv(cov1) * A1) * A1' * pinv(cov1) * p1;
        cov_theta1 = pinv(A1' * pinv(cov1) * A1);

        % Second step of the iCHAN algorithm
        A2 = [1, 0; 0, 1; 1, 1];
        p2 = [theta1(1,1)^2; theta1(2,1)^2; theta1(3,1)];
        B2 = diag([2*theta1(1,1); 2*theta1(2,1); 1]);
        cov2 = B2 * cov_theta1 * B2;
        theta2 = pinv(A2' * pinv(cov2) * A2) * A2' * pinv(cov2) * p2;
        theta2 = abs(theta2);

        % Update current estimate of MS position
        OUT_curr = sign(theta1(1:2,:)) .* theta2.^(1/2);
        OUT_curr = OUT_curr';

        % Calculate RMSE between current and previous estimates
        RMSE_curr = sqrt((OUT_curr(1) - OUT_prev(1))^2 + (OUT_curr(2) - OUT_prev(2))^2);

        % Check for convergence
        if RMSE_curr < step
            OUT_prev = OUT_curr;
            break; % Exit loop if convergence is achieved
        end

        % Update previous estimate
        OUT_prev = OUT_curr;
    end
end