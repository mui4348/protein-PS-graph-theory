clear
% Define file names for the two cases
fileNames = {'wt', 'wt_unpack'};
for kk = 1:numel(fileNames)
for ii = 270:10:340 % temperatures

temp = int2str(ii);
%mt_alpha0.925_60chain_3C_protNoPBC_1k_lst500frms_290K.pdb
fileName1 = strcat(fileNames{kk},'_1C_protNoPBC_1k_',temp,'K.pdb');
fileName2 = strcat(fileNames{kk},'_2C_protNoPBC_1k_',temp,'K.pdb');
fileName3 = strcat(fileNames{kk},'_3C_protNoPBC_1k_',temp,'K.pdb');

threshold = 8;
Nodes_with_linkss1 = chainLinkss(fileName1, threshold);
[coords1, residues1] = calculatecoords(fileName1);
Nodes_with_linkss2 = chainLinkss(fileName2, threshold);
[coords2, residues2] = calculatecoords(fileName2);
Nodes_with_linkss3 = chainLinkss(fileName3, threshold);
[coords3, residues3] = calculatecoords(fileName3);

%residues;

Nodes_with_linkss = [Nodes_with_linkss1,Nodes_with_linkss2,Nodes_with_linkss3];
coords = [coords1,coords2,coords3];


for i = 1:numel(Nodes_with_linkss)
    for j = 1:numel(Nodes_with_linkss{i})
        chainj = Nodes_with_linkss{i}{j}(1);
        for k = 2:numel(Nodes_with_linkss{i}{j})
            chaink = Nodes_with_linkss{i}{j}(k);
            [Indices_row, Indices_col] = calculatePairDataidices(coords{i}, threshold, chainj, chaink);
            interating_nodes_framei_chainjk{i}{j}{k} = unique([transpose(Indices_col), transpose(Indices_row)]);
            Indices_chainj_node{i}{j}{k} = num2cell(transpose(unique(Indices_row)));
            Indices_chaink_node{i}{j}{k} = num2cell(transpose(unique(Indices_col)));
        end
    end
end

% Initialize an empty array to store all numbers
all_numbers = [];

% Loop through the cell array to extract all numbers
for i = 1:length(interating_nodes_framei_chainjk)
    for j = 1:length(interating_nodes_framei_chainjk{i})
        for k = 1:length(interating_nodes_framei_chainjk{i}{j})
            % Extract the numeric matrix
            current_matrix = interating_nodes_framei_chainjk{i}{j}{k};
            
            % Append the numbers to the all_numbers array
            all_numbers = [all_numbers; current_matrix(:)];
        end
    end
end
%%
% Compute the histogram counts
num_bins = 33; % Number of bins corresponding to your 33 labels (residues per chain)
[counts, edges] = histcounts(all_numbers, num_bins);

% Convert counts to fractions
total_count = sum(counts);
fractions = counts / total_count;

% Normalize fractions to the range 0 to 1
max_fraction = max(fractions);
normalized_fractions = fractions / max_fraction;
%%
% Write the results to a CSV file
csv_filename = strcat('chain_residue_connections_',fileNames{kk},'_',temp,'K.csv');
%writetable(variable,strcat(Variable_table{i},'_mt_3k_',temp,'K.csv'));

fid = fopen(csv_filename, 'w');
fprintf(fid, 'Residue,Count\n');
for i = 1:num_bins
    fprintf(fid, '%s,%d\n', residues1{i}, normalized_fractions(i));
end
fclose(fid);
disp(['Data has been written to ' csv_filename]);

end
end

