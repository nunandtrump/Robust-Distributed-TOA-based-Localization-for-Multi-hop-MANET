function [CRLB, FIM_inverse] = TOA_local_CRLB(LS, targ, gauss_distance, gauss_anchor)
    % TOA_local_CRLB: Compute the local Cramér-Rao Lower Bound (CRLB) for Time of Arrival (TOA) based localization.
    % This function calculates the CRLB for the position estimation of a target using TOA measurements.
    %
    % Inputs:
    %   LS            : Matrix of anchor positions (N x 2), where N is the number of anchors.
    %   targ          : Target position (1 x 2).
    %   gauss_distance: Variance of the TOA measurement noise.
    %   gauss_anchor  : Variance of the anchor position uncertainty (N x 1).
    %
    % Outputs:
    %   CRLB          : Cramér-Rao Lower Bound for the target position estimation.
    %   FIM_inverse   : Inverse of the Fisher Information Matrix (FIM).
    %
    % Description:
    % The function computes the CRLB for TOA-based localization, considering both measurement noise
    % and anchor position uncertainty. The CRLB provides a lower bound on the variance of any unbiased
    % estimator of the target position.

    N = size(LS, 1);  % Number of anchors
    FIM = zeros(2*(N+1), 2*(N+1));  % Initialize the Fisher Information Matrix (FIM)

    % Precompute distances and differences between target and anchors
    r = zeros(N, 1);  % Distances between target and anchors
    dx = zeros(N, 1); % Differences in x-coordinates
    dy = zeros(N, 1); % Differences in y-coordinates

    for i = 1:N
        r(i) = sqrt((LS(i,:) - targ) * (LS(i,:) - targ)');  % Euclidean distance
        dx(i) = targ(1) - LS(i, 1);  % Difference in x-coordinates
        dy(i) = targ(2) - LS(i, 2);  % Difference in y-coordinates
    end

    % Fill the upper-left 2x2 block of the FIM (related to the target position)
    FIM(1,1) = sum((dx.^2 ./ r.^2) / (gauss_distance^2));  % FIM(1,1)
    FIM(1,2) = sum((dx .* dy ./ r.^2) / (gauss_distance^2));  % FIM(1,2)
    FIM(2,1) = FIM(1,2);  % FIM(2,1) is symmetric to FIM(1,2)
    FIM(2,2) = sum((dy.^2 ./ r.^2) / (gauss_distance^2));  % FIM(2,2)

    % Fill the remaining blocks of the FIM (related to anchor positions and their uncertainties)
    for i = 3:(2*(N+1))
        if mod(i, 2) == 0
            % Even indices correspond to y-coordinates of anchors
            num = i/2 - 1;  % Anchor index
            FIM(1,i) = -(dx(num) * dy(num) / r(num)^2) / (gauss_distance^2);  % FIM(1,i)
            FIM(i,1) = FIM(1,i);  % Symmetric entry
            FIM(2,i) = -(dy(num) * dy(num) / r(num)^2) / (gauss_distance^2);  % FIM(2,i)
            FIM(i,2) = FIM(2,i);  % Symmetric entry
            FIM(i,i) = 2/(gauss_anchor(num)) + (dy(num) * dy(num) / r(num)^2) / (gauss_distance^2);  % FIM(i,i)
        else
            % Odd indices correspond to x-coordinates of anchors
            num = (i - 1)/2;  % Anchor index
            FIM(i,i+1) = (dx(num) * dy(num) / r(num)^2) / (gauss_distance^2);  % FIM(i,i+1)
            FIM(i+1,i) = FIM(i,i+1);  % Symmetric entry
            FIM(1,i) = -(dx(num) * dx(num) / r(num)^2) / (gauss_distance^2);  % FIM(1,i)
            FIM(i,1) = FIM(1,i);  % Symmetric entry
            FIM(2,i) = -(dy(num) * dx(num) / r(num)^2) / (gauss_distance^2);  % FIM(2,i)
            FIM(i,2) = FIM(2,i);  % Symmetric entry
            FIM(i,i) = 2/(gauss_anchor(num)) + (dx(num) * dx(num) / r(num)^2) / (gauss_distance^2);  % FIM(i,i)
        end
    end

    % Compute the inverse of the FIM
    FIM_inverse = inv(FIM);

    % Extract the top-left 2x2 block of the inverse FIM (related to the target position)
    top_left_2x2 = FIM_inverse(1:2, 1:2);

    % Compute the CRLB as the trace of the top-left 2x2 block
    CRLB = trace(top_left_2x2);  % Cramér-Rao Lower Bound
end