function [err] = meg_PLVclassify(Class_Vec, Mahal_Vec, meg_tw)
   
    feat = [];
    feat = Mahal_Vec(:,1);
    feat = feat';

    for zz = 1:size(meg_tw,1)
        CL{zz} = meg_tw{zz,2};
    end
    % size(CL) needs to be 4996 x 1
    % size(x_pow_str.data) = 4996 x 34
    % so we have 34 dim data

    % FIND CVPARTITION1 AND ADD TO PATH
    % cvpartition1 -> boolean value seperates train & test data 80/20 and
    % creates 5-fold data to test on
    [ testIdx, trainIdx, nTest] = cvpartition1(size(feat',1), 5);

   C = cell(5,1); err = zeros(5,1); P = cell(5,1); logp = cell(5,1); coeff = cell(5,1);
        for i = 1:5 
            % Classify the data
            % classify( test data, training data, group labels );
            % size(group) -> 59952 x 1
            % size(g) -> 5 x 59952 (4996*12=59952);    size(find(g(1,:))) -> 1 x 47952

          [C{i, 1},err(i, 1),P{i, 1},logp{i, 1},coeff{i, 1}] = ...
                                                                    classify( feat(testIdx(i,:))', ...
                                                                    feat(trainIdx(i,:))' , ...
                                                                    Class_Vec(trainIdx(i,:)), 'linear' );

        end
        
end

