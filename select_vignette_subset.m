function select_vignette_subset(infile, outname, limits, breakdown, nsurvey, heuristic, seed)
%
%   infile : string
%   outname : string
%   limits : 1 x 4 array of floats
%   breakdown : 1 x 3 array of integers
%   nsurvey : integer or 'inf'
%   heuristic : string {'random', 'heuristic'}
%   seed : integer to seed random number generator
%
%   Copyright 2014, Erin K. Molloy
%   mailto://emolloy2@illinois.edu

    % Set seed for random number generator
    rng(seed)

    % Import full survey data
    Q = importdata(infile);
    X = Q.data;
    [m, n] = size(X);
    header = Q.colheaders;

    % Define bins: low [a,b], neutral (b,d), high [d,e]
    l = size(limits);
    if (l(2) == 4)
        a = repmat(limits(1), 1, n);
        b = repmat(limits(2), 1, n);
        d = repmat(limits(3), 1, n);
        e = repmat(limits(4), 1, n);
    else
        a = min(X);  e = max(X);
        c = mean(X); s = std(X);
        b = c - s;   d = c + s;
    end

    if strcmp('random', heuristic)
        order = randperm(m);
        nq = (n - 1) * size(breakdown,2);
    else
        % Create selection masks
        low_mask  = zeros(m, n);
        mid_mask  = zeros(m, n);
        high_mask = zeros(m, n);
        for i = 2:n
            low_mask(:,i)  = [X(:,i) >= a(i)] .* [X(:,i) <= b(i)];
            mid_mask(:,i)  = [X(:,i) >  b(i)] .* [X(:,i) <  d(i)];
            high_mask(:,i) = [X(:,i) >= d(i)] .* [X(:,i) <= e(i)];
        end
        masks = cat(3, low_mask, mid_mask, high_mask);
    end

    % Create vignette subsets
    i = 1;
    while (i <= nsurvey)
        Y = [];

        if strcmp('random', heuristic)
            [X, Y] = create_subset_random(X, Y, i, nq, order);
            outnamei = strcat(outname, '_', int2str(i), '_random');
        else
            [X, Y] = create_subset_heuristic(X, Y, breakdown, masks);
            outnamei = strcat(outname, '_', int2str(i), '_heuristic');
        end

        save_subset(Y, header, outnamei, a, b, d, e);
        i = i + 1;
    end
end


function [A, B] = create_subset_random(A, B, i, nq, order)
    st = (i-1) * nq + 1;
    en = st + nq - 1;
    if (en < size(A, 1))
        B = A(order(st:en),:);
    else
        error('Error: No elements in bin... quiting!')
    end
end


function [A, B] = create_subset_heuristic(A, B, breakdown, masks)
    [m, n] = size(A);

    levels = cat(2, ones(1,breakdown(1)), ...
                    ones(1,breakdown(2))*2, ...
                    ones(1,breakdown(3))*3);
    order = randperm(size(levels,2));
    levels = levels(order);

    for i = levels
        if (i == 1)
            mask = masks(:,:,1);
        elseif (i == 2)
            mask = masks(:,:,2);
        else
            mask = masks(:,:,3);
        end

        order = randperm(n-1);
        factors = order + 1;
        for j = factors
            subset = A(:,1) .* mask(:,j);
            subset(subset==0) = [];
            
            % Identifies and adds question to subset
            l = size(subset,1);
            if (l ~= 0)
                randsel = randi([1,l]);
                question = subset(randsel);
            
                % Add question to B
                B = cat(1, B, A(question,:));
            
                % Remove question from A
                A(question, 1) = 0;
            else
                error('Error: No elements in bin... quiting!')
            end
        end
    end
end


function save_subset(B, header, outname, a, b, d, e)
    % Save subset of vignettes
    outfile = strcat(outname, '.csv');
    T = array2table(B, 'VariableNames', header);
    writetable(T, outfile, 'Delimiter', ',');
            
    % Compute and save associated metrics
    n = size(B, 2);
    average = mean(B,1);
    stddev = std(B,1);
    minimum = zeros(1,n);
    maximum = zeros(1,n);
    total = zeros(3,n);
    for i = 2:n
        minimum(i) = min(B(:,i));
        maximum(i) = max(B(:,i));
        total(1,i) = sum([B(:,i) >= a(i)] .* [B(:,i) <= b(i)]);
        total(2,i) = sum([B(:,i) >  b(i)] .* [B(:,i) <  d(i)]);
        total(3,i) = sum([B(:,i) >= d(i)] .* [B(:,i) <= e(i)]);
    end
    ancillary = cat(1, average(:,2:n), stddev(:,2:n), ...
                       minimum(:,2:n), maximum(:,2:n), total(:,2:n), ...
                       a(:,2:n), b(:,2:n), d(:,2:n), e(:,2:n));

    outfile = strcat(outname, '_info.csv');
    rownames = {'mean'; 'stddev'; 'min'; 'max'; ...
                'nlow'; 'nmid'; 'nhigh'; ...
                'a'; 'b'; 'd'; 'e'};
    T = array2table(ancillary, 'VariableNames', header(2:n), ...
                               'RowNames', rownames);
    writetable(T, outfile, 'Delimiter', ',', 'WriteRowNames', true);
end
