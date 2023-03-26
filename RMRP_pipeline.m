% Matlab script for illustrating the secondary structure of RMRP's 5'UTR
% Define sequence
RMRP_seq = 'GUUCGUGCUGAAGGCCUGUAUCCUAGGCUACACACUGAGGACUCUGUUCCUCCCCUUUCCGCCUAGGGGAAAGUCCCCGGACCUCGGGCAGAGAGUGCCACGUGCAUACGCACGUAGACAUUCCCCGCUUCCCACUCCAAAGUCCGCCAAGAAGCGUAUCCCGCUGAGCGGCGUGGCGCGGGGGCGUCAUCCGUCAGCUCCCUCUAGUUACGCAGGCAGUGCGUGUCCGCGCACCAACCACACGGGGCUCAUUCUCAGCGCGGC';

%3 known pathogenic SNP variations induced on reference sequence
String1 = RMRP_seq;
String1(252) = 'G';
String2 = RMRP_seq;
String2(40) = 'A';
String3 = RMRP_seq;
String3(182) = 'U';

%Their sequences are saved as strings and structures saved as WT's
wild_type_1 = rnafold(String1);
wild_type_2 = rnafold(String2);
wild_type_3 = rnafold(String3);

% Use rnafold to predict optimal secondary structure
h = waitbar(0,'Calculating RNA Secondary Structure...');
RNAbracket = rnafold(RMRP_seq, 'MinLoopSize', 7);
waitbar(1,h,'Calculating RNA Secondary Structure...');
close(h)
% Use rnaplot to generate structural diagram
% Use rnaplot to generate the dotdiagram of the raw sequence
h = waitbar(0, 'Generating Structural Diagrams...');
ha1 = rnaplot(RNAbracket, 'Format', 'Dotdiagram');
title(ha1, 'Structural Dot-Diagram of Raw Sequence')
waitbar(0.25,h,'Generating Structural Diagrams...');
% Use rnaplot to generate the circle diagram of the raw sequence
ha2 = rnaplot(RNAbracket);
title(ha2, 'Structural Circle diagram of Raw Sequence')
waitbar(0.5,h,'Generating Structural Diagrams...');
% Use rnafold to find the minimum free-energy secondary structure
smalldeltaG = rnafold(RMRP_seq);
% Use rnaplot to plot the minimum secondary structure dotdiagram
ha3 = rnaplot(smalldeltaG, 'Format', 'Dotdiagram');
title(ha3, 'Structural Dot-Diagram of MFE Sequence')
waitbar(0.75,h,'Generating Structural Diagrams...');
% Use rnaplot to plot the minimum secondary structure circle diagram
ha4 = rnaplot(smalldeltaG);
title(ha4, 'Structural Circle Diagram of MFE Sequence')
waitbar(1,h,'Generating Structural Diagrams...');
close(h)
% Use rnaplot to generate the dotdiagram of the wildtype 1
ha5 = rnaplot(wild_type_1, 'Format', 'Dotdiagram', 'selection', 14);
title(ha5, 'Structural Dot-Diagram of Wild Type 1')
% Use rnaplot to generate the circle diagram of the wildtype 1
ha6 = rnaplot(wild_type_1);
title(ha6, 'Structural Circle diagram of Wild Type 1')
% Use rnaplot to generate the dotdiagram of the wildtype 2
ha7 = rnaplot(wild_type_2, 'Format', 'Dotdiagram', 'selection', 56);
title(ha7, 'Structural Dot-Diagram of Wild Type 2')
% Use rnaplot to generate the circle diagram of the wildtype 2
ha8 = rnaplot(wild_type_2);
title(ha8, 'Structural Circle diagram of Wild Type 2')
% Use rnaplot to generate the dotdiagram of the wildtype 3
ha9 = rnaplot(wild_type_3, 'Format', 'Dotdiagram', 'selection', 22);
title(ha9, 'Structural Dot-Diagram of Wild Type 3')
% Use rnaplot to generate the circle diagram of the wildtype 3
ha10 = rnaplot(wild_type_3);
title(ha10, 'Structural Circle diagram of Wild Type 3')

%Print similarity scores to command window
fprintf('Predicted Structure v WT1 = \n\n');
score_comp(RNAbracket,wild_type_1)
fprintf('Predicted Structure v WT2 = \n\n');
score_comp(RNAbracket,wild_type_2)
fprintf('Predicted Structure v WT3 = \n\n');
score_comp(RNAbracket,wild_type_3)
fprintf('WT1 v WT2 = \n\n');
score_comp(wild_type_1,wild_type_2)
fprintf('WT2 v WT3 = \n\n');
score_comp(wild_type_2,wild_type_3)
fprintf('WT1 v WT3 = \n\n');
score_comp(wild_type_1, wild_type_3)
%%
% Define structure
RMRP_struct = rnafold(RMRP_seq);
% Mutate every position in the RNA with every possible nucleotide
mut_table = cell(length(RMRP_seq)*4, 6);
h = waitbar(0, 'Generating mutation table...');
for i = 1:length(RMRP_seq)
    waitbar(i/length(RMRP_seq), h, sprintf('Mutating position %d of %d...', i, length(RMRP_seq)));
    for n = 1:4
        mutated_seq = RMRP_seq;
        if n == 1
            mutated_seq(i) = 'A';
            nucleotide = 'A';
        elseif n == 2
            mutated_seq(i) = 'C';
            nucleotide = 'C';
        elseif n == 3
            mutated_seq(i) = 'G';
            nucleotide = 'G';
        elseif n == 4
            mutated_seq(i) = 'U';
            nucleotide = 'U';
        end
        mutated_struct = rnafold(mutated_seq);
        score = score_comp(RMRP_struct, mutated_struct);
        % Determine pairdness from dot-bracket notation changes from reference structure
        original_pairdness = 'Unpaired';
        if RMRP_struct(i) == '('
            original_pairdness = 'Paired';
        end
        mutated_pairdness = 'Unpaired';
        if mutated_struct(i) == '('
            mutated_pairdness = 'Paired';
        end
        row = {i, RMRP_seq(i), nucleotide, score, original_pairdness, mutated_pairdness};
        mut_table(i*4+n, :) = row; % Note the change in this line
    end
end
close(h);
% Save table to a tab-separated file
fid = fopen('RMRP_gsgeorge.txt', 'w');
fprintf(fid, 'Mutated Position\tOriginal Nucleotide at this position\tNucleotide was mutated to this nucleotide\tComparison score (higher is more similar)\tOriginal Nucleotide Pairdness\tMutated Nucleotide Pairdness\n');
for i = 1:length(mut_table)
    fprintf(fid, '%d\t%s\t%s\t%f\t%s\t%s\n', mut_table{i,:});
end
fclose(fid);
%%
% Read in the data
fid = fopen('RMRP_gsgeorge.txt');
data = textscan(fid, '%d %s %s %f %s %s', 'Delimiter', '\t', 'HeaderLines', 1);

h1 = waitbar(0.2, 'The data is being read...');

% Calculate the minimum comparison score
min_score = min(data{4});

% Create a waitbar to inform the user that the minimum score is being calculated
h2 = waitbar(0.4, 'The minimum score is being calculated...');

% Determine the location of mutation with the lowest comparison score
lowest_score_pos = data{1}(data{4} == min_score);

% Create a waitbar to inform the user that the lowest score position is being determined
h3 = waitbar(0.6, 'The lowest score position is being determined...');

% Output the nucleotide change
fprintf('The mutation at position %d has the lowest comparison score, with a score of %f. The nucleotide was changed from %s to %s.\n', ...
    lowest_score_pos, min_score, data{2}{data{4} == min_score}, data{3}{data{4} == min_score});
%Output the most mutated structure with most mutated position highlighted
fclose(fid);

% Create a waitbar to inform the user that the raw sequence is being plotted
h4 = waitbar(0.8, 'The raw sequence is being plotted...');

% Use rnaplot to generate the circle diagram of the raw sequence
ha2 = rnaplot(RNAbracket, 'Selection', lowest_score_pos);
title(ha2, 'Structural Circle diagram of Raw Sequence with Lowest Score Position Highlighted')

% Use rnaplot to generate the dotdiagram of the raw sequence
ha1 = rnaplot(RNAbracket, 'Format', 'Dotdiagram', 'Selection', lowest_score_pos);
title(ha1, 'Structural Dot-Diagram of Raw Sequence with Lowest Score Position Highlighted')

% Close the waitbar when the script is finished
close(h4);

% Create a waitbar to inform the user that the lowest score is being mutated
h5 = waitbar(0.9, 'The lowest score is being mutated...');

% Mutate the lowest score position to the nucleotide given
mutated_seq = RMRP_seq;
mutated_seq(lowest_score_pos) = data{3}{data{4} == min_score};

% Create a waitbar to inform the user that the mutated sequence is being plotted
h6 = waitbar(0.95, 'The mutated sequence is being plotted...');

% Use rnaplot to generate the circle diagram of the mutated sequence
ha3 = rnaplot(rnafold(mutated_seq), 'Selection', lowest_score_pos);
title(ha3, 'Structural Circle diagram of Mutated Sequence with Lowest Score Position Highlighted')

% Use rnaplot to generate the dotdiagram of the mutated sequence
ha4 = rnaplot(rnafold(mutated_seq), 'Format', 'Dotdiagram', 'Selection', lowest_score_pos);
title(ha4, 'Structural Dot-Diagram of Mutated Sequence with Lowest Score Position Highlighted')

% Close the waitbar when the script is finished
close(h6);
close(h5);
close(h3);
close(h2);
close(h1);

%%
% Create a waitbar for user to track progress
h = waitbar(0, 'Calculating RMRP structure...');

RMRP_struct = rnafold(RMRP_seq);
waitbar(0.25, h);

% Create a histogram of the similarity scores
figure;
histogram(data{4});
xlabel('Score');
ylabel('Number of Mutations');
title('Histogram of Mutation Scores');
waitbar(0.5, h);

% Determine where the wild-type comparison scores fall on the distribution
wild_type_scores = [score_comp(RMRP_struct, wild_type_1), score_comp(RMRP_struct, wild_type_2), score_comp(RMRP_struct, wild_type_3)];
hold on;
scatter(wild_type_scores(1), ones(1,1)*max(data{4})/2, 'filled', 'r');
scatter(wild_type_scores(2), ones(1,1)*max(data{4})/2, 'filled', 'g');
scatter(wild_type_scores(3), ones(1,1)*max(data{4})/2, 'filled', 'b');
legend('Mutations', 'Wild-Type 1', 'Wild-Type 2', 'Wild-Type 3');
hold off;
waitbar(1, h);

% Close the waitbar
close(h);
%%
% Create a waitbar for user to track progress
h = waitbar(0, 'Calculating RMRP structure...');

RMRP_struct = rnafold(RMRP_seq);
waitbar(0.25, h);

% Create a kernel distribution of the similarity scores
figure;
histfit(data{4}); %I might need to change this to better fit my data in the future, possibly implementing beta distribution
xlabel('Score');
ylabel('Number of Mutations');
title('Kernel Distribution of Mutation Scores');
waitbar(0.5, h);

% Determine where the wild-type comparison scores fall on the distribution
wild_type_scores = [score_comp(RMRP_struct, wild_type_1), score_comp(RMRP_struct, wild_type_2), score_comp(RMRP_struct, wild_type_3)];
hold on;
scatter(wild_type_scores(1), ones(1,1)*max(data{4})/2, 'filled', 'r');
scatter(wild_type_scores(2), ones(1,1)*max(data{4})/2, 'filled', 'g');
scatter(wild_type_scores(3), ones(1,1)*max(data{4})/2, 'filled', 'b');
legend('Mutations', 'Standard Distribution','Wild-Type 1', 'Wild-Type 2', 'Wild-Type 3');
hold off;
waitbar(1, h);

% Close the waitbar
close(h);

%%
fid = fopen('RMRP_gsgeorge.txt', 'r');
file_data = textscan(fid, '%d %s %s %f %s %s', 'HeaderLines', 1);
fclose(fid);

[sorted_scores, sorted_idx] = sort(file_data{4});

for i=1:10
    fprintf('%d\t%s\t%s\t%f\t%s\t%s\n', file_data{1}(sorted_idx(i)), ...
        file_data{2}{sorted_idx(i)}, file_data{3}{sorted_idx(i)}, ...
        file_data{4}(sorted_idx(i)), file_data{5}{sorted_idx(i)}, ...
        file_data{6}{sorted_idx(i)});
end
%%
% read in data
fid = fopen('RMRP_gsgeorge.txt', 'r');
file_data = textscan(fid, '%d %s %s %f %s %s', 'HeaderLines', 1);
fclose(fid);

% calculate averages
paired_to_unpaired_score = mean(file_data{4}(strcmp(file_data{5}, 'Paired') & strcmp(file_data{6}, 'Unpaired')));
unpaired_to_paired_score = mean(file_data{4}(strcmp(file_data{5}, 'Unpaired') & strcmp(file_data{6}, 'Paired')));
no_change_score = mean(file_data{4}(strcmp(file_data{5}, file_data{6})));
overall_score = mean(file_data{4});

% print out results
fprintf('Average Similarity Score for Mutations that went from Paired to Unpaired: %f\n', paired_to_unpaired_score);
fprintf('Average Similarity Score for Mutations that went from Unpaired to Paired: %f\n', unpaired_to_paired_score);
fprintf('Average Similarity Score for Mutations that did not result in a change in Pairedness: %f\n', no_change_score);
fprintf('Average Similarity Score for all Mutations: %f\n', overall_score);
%%
% read in data
fid = fopen('RMRP_gsgeorge.txt', 'r');
file_data = textscan(fid, '%d %s %s %f %s %s', 'HeaderLines', 1);
fclose(fid);

% calculate averages
paired_to_unpaired_score = mean(file_data{4}(strcmp(file_data{5}, 'Paired') & strcmp(file_data{6}, 'Unpaired')));
unpaired_to_paired_score = mean(file_data{4}(strcmp(file_data{5}, 'Unpaired') & strcmp(file_data{6}, 'Paired')));
no_change_score = mean(file_data{4}(strcmp(file_data{5}, file_data{6})));
overall_score = mean(file_data{4});

% create array of data points
data_points = [paired_to_unpaired_score, unpaired_to_paired_score, no_change_score, overall_score];

% create bar graph
figure;
bar(data_points, 'FaceColor', [0 0.5 0.9]);
xlabel('Data Type');
ylabel('Average Similarity Score');
title('Average Similarity Score for Mutations');

set(gca, 'xticklabel', {'Paired to Unpaired', 'Unpaired to Paired', 'No Change', 'Overall'});

% %% THIS IS COMMENTED OUT BECAUSE IT INVOLVES AMINO ACID STRUCTURE (Too
% pretty to delete!)
% %The script will now analize the amino acid structure of these sequences
% % Define sequence for RMRP
% RMRP_seq = 'GUUCGUGCUGAAGGCCUGUAUCCUAGGCUACACACUGAGGACUCUGUUCCUCCCCUUUCCGCCUAGGGGAAAGUCCCCGGACCUCGGGCAGAGAGUGCCACGUGCAUACGCACGUAGACAUUCCCCGCUUCCCACUCCAAAGUCCGCCAAGAAGCGUAUCCCGCUGAGCGGCGUGGCGCGGGGGCGUCAUCCGUCAGCUCCCUCUAGUUACGCAGGCAGUGCGUGUCCGCGCACCAACCACACGGGGCUCAUUCUCAGCGCGGC';
% 
% %Read the reference sequence
% REF_seq = 'GUUCGUGCUGAAGGCCUGUAUCCUAGGCUACACACUGAGGACUCUGUUCCUCCCCUUUCCGCCUAGGGGAAAGUCCCCGGACCUCGGGCAGAGAGUGCCACGUGCAUACGCACGUAGACAUUCCCCGCUUCCCACUCCAAAGUCCGCCAAGAAGCGUAUCCCGCUGAGCGGCGUGGCGCGGGGGCGUCAUCCGUCAGCUCCCUCUAGUUACGCAGGCAGUGCGUGUCCGCGCACCAACCACACGGGGCUCAUUCUCAGCGCGGC';
% 
% String1 = RMRP_seq;
% String1(252) = 'G';
% String2 = RMRP_seq;
% String2(40) = 'A';
% String3 = RMRP_seq;
% String3(182) = 'U';
% 
% % Analyze protein translation
% % Translate the sequences into their corresponding amino acid sequences
% wild_type_1_AA = nt2aa(String1);
% wild_type_2_AA = nt2aa(String2);
% wild_type_3_AA = nt2aa(String3);
% ref_type_AA = nt2aa(REF_seq);
% 
% % Compare the sequences to determine if there are any differences
% wild_type_1_diff = ntdiff(wild_type_1_AA, ref_type_AA);
% wild_type_2_diff = ntdiff(wild_type_2_AA, ref_type_AA);
% wild_type_3_diff = ntdiff(wild_type_3_AA, ref_type_AA);
% 
% % Analyze secondary structure
% % Predict the secondary structure of the FTL sequence
% wild_type_1 = rnafold(String1);
% wild_type_2 = rnafold(String2);
% wild_type_3 = rnafold(String3);
% ref_type = rnafold(REF_seq);
% 
% % Compare the secondary structures for each sequence
% wt1_diff = ntdiff(wild_type_1, ref_type);
% wt2_diff = ntdiff(wild_type_2, ref_type);
% wt3_diff = ntdiff(wild_type_3, ref_type);
% 
% % Analyze mRNA stability
% % Calculate the free energy of the sequences
% [wt1_energy, wt1_structure] = rnafoldenergy(String1);
% [wt2_energy, wt2_structure] = rnafoldenergy(String2);
% [wt3_energy, wt3_structure] = rnafoldenergy(String3);
% [ref_energy, ref_structure] = rnafoldenergy(REF_seq);
% 
% % Compare the free energies
% wt1_energy_diff = wt1_energy - ref_energy;
% wt2_energy_diff = wt2_energy - ref_energy;
% wt3_energy_diff = wt3_energy - ref_energy;
% 
% % Analyze binding of mRNA to regulatory proteins
% % Determine the regulatory proteins that bind to the FTL sequence
% wt1_regulatory_proteins = findRegulatoryProteins(String1);
% wt2_regulatory_proteins = findRegulatoryProteins(String2);
% wt3_regulatory_proteins = findRegulatoryProteins(String3);
% ref_regulatory_proteins = findRegulatoryProteins(REF_seq);
% 
% % Compare the regulatory proteins
% wt1_regulatory_diff = setdiff(wt1_regulatory_proteins, ref_regulatory_proteins);
% wt2_regulatory_diff = setdiff(wt2_regulatory_proteins, ref_regulatory_proteins);
% 
% % Plot the protein structures
% wt1_struct = proteinplot(wild_type_1_AA);
% wt2_struct = proteinplot(wild_type_2_AA);
% wt3_struct = proteinplot(wild_type_3_AA);
% ref_struct = proteinplot(ref_type_AA);
%% Prior to reconstruction of mut_table in parallel, clear all variables from matlab enviornment
clearvars;

% Define sequence
RMRP_seq = 'GUUCGUGCUGAAGGCCUGUAUCCUAGGCUACACACUGAGGACUCUGUUCCUCCCCUUUCCGCCUAGGGGAAAGUCCCCGGACCUCGGGCAGAGAGUGCCACGUGCAUACGCACGUAGACAUUCCCCGCUUCCCACUCCAAAGUCCGCCAAGAAGCGUAUCCCGCUGAGCGGCGUGGCGCGGGGGCGUCAUCCGUCAGCUCCCUCUAGUUACGCAGGCAGUGCGUGUCCGCGCACCAACCACACGGGGCUCAUUCUCAGCGCGGC';

% Define structure
RMRP_struct = rnafold(RMRP_seq);

% Calculate similarity score for every single nucleotide change to RMRP_seq
score_matrix = zeros(4,length(RMRP_seq));
parfor i = 1:length(RMRP_seq)
    for j = 1:4
        temp_seq = RMRP_seq;
        if j == 1
            temp_seq(i) = 'A';
        elseif j == 2
            temp_seq(i) = 'C';
        elseif j == 3
            temp_seq(i) = 'G';
        elseif j == 4
            temp_seq(i) = 'U';
        end
        temp_struct = rnafold(temp_seq);
        score_matrix(j,i) = score_comp(temp_struct, RMRP_struct);
    end
end

% Generate heatmap of single nucleotide change scores
nucs = {'A', 'C', 'G', 'U'};
figure;
imagesc(score_matrix);
colormap(jet);
colorbar;
title('Single Nucleotide Change Scores');
xlabel('Position');
ylabel('Nucleotide Change');
set(gca, 'YTick', 1:4);
set(gca, 'YTickLabel', nucs);