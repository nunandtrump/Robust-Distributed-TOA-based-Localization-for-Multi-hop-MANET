function collinearity = calculateCollinearity(anchor_node, anchor_node_CRLB, gauss_distance)
    % calculateCollinearity: Determine if anchor nodes are collinear based on their error circles.
    % This function checks whether the anchor nodes are collinear by examining if their error circles
    % overlap in such a way that a common intersection point exists.
    %
    % Inputs:
    %   anchor_node      : Matrix of anchor positions (N x 2), where N is the number of anchors.
    %   anchor_node_CRLB  : CRLB for each anchor's position uncertainty (N x 1).
    %   gauss_distance   : Variance of the TOA measurement noise (scalar).
    %
    % Output:
    %   collinearity     : A binary value indicating collinearity:
    %                       - 0: Anchor nodes are not collinear (a common intersection point exists).
    %                       - 1: Anchor nodes are collinear (no common intersection point exists).
    %
    % Description:
    % The function calculates the error circle radius for each anchor node based on the measurement noise
    % and anchor position uncertainty. It then checks if there exists a common intersection point for all
    % error circles. If such a point exists, the anchor nodes are considered non-collinear; otherwise, they
    % are considered collinear.

    % Number of anchor nodes
    num_anchor = size(anchor_node, 1);

    % Calculate the radius of the error circle for each anchor node
    radius_anchor_error = zeros(1, num_anchor);
    for i = 1:num_anchor
        radius_anchor_error(i) = 2 * (gauss_distance + sqrt(anchor_node_CRLB(i)/2));
    end

    % Symbolic variable for theta (angle)
    syms theta;

    % Initialize counter and cell array to store inequalities
    count = 0;
    inequalities = cell(1, 0.5 * num_anchor * (num_anchor - 1));

    % Generate inequalities for each pair of anchor nodes
    for m = 1:(num_anchor - 1)
        for n = (m + 1):num_anchor
            count = count + 1;
            c = anchor_node(m, 1) - anchor_node(n, 1); % xm - xn
            d = anchor_node(m, 2) - anchor_node(n, 2); % ym - yn
            e = radius_anchor_error(m) + radius_anchor_error(n); % rm + rn
            inequalities{count} = abs(c * cos(theta) + d * sin(theta)) <= e;
        end
    end

    % Define range and step size for theta
    lower_bound = 0;
    upper_bound = 3.14;
    step_size = 0.1;

    % Search for a value of theta that satisfies all inequalities
    found_point = false;
    for val = lower_bound:step_size:upper_bound
        is_satisfied = true;
        for i = 1:0.5 * num_anchor * (num_anchor - 1)
            ineq_val = subs(inequalities{i}, theta, val); % Substitute theta with the current value
            is_satisfied = is_satisfied && isAlways(ineq_val); % Check if the inequality holds
        end
        if is_satisfied
            found_point = true;
            collinearity = 0; % Anchor nodes are not collinear
            break; % Exit loop if a satisfying point is found
        end
    end

    % If no satisfying point is found, anchor nodes are collinear
    if ~found_point
        collinearity = 1; % Anchor nodes are collinear
    end
end
