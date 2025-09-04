clear
% Define file names for the two cases
Names = { 'wt', 'wt_unpack' };
%Names = { 'MaSpI1' };
num_trajectories = 3;
for ii =  270:10:300
%ii = 320;
for k = 1:numel(Names)
    Name = Names{k};
% Loop through the 3 trajectories
for j = 1:num_trajectories
    % Define file parameters
    temp = int2str(ii);
    traj = int2str(j);
    fileName = strcat(Name, '_', traj, 'C_protNoPBC_1k_', temp, 'K.pdb');
	display(['Processing file: ', fileName]);
threshold_1 = 8;
chainNum = 60;
[ChainPair,sig_node_count] = ChainPairDFF(fileName,threshold_1);
Chain_Adjacency_matrices = {};

for i=1:numel(ChainPair)
    logical_indices{i} = ChainPair{i} < threshold_1 & ChainPair{i} ~= 0;
    Chain_Adjacency_matrices{i} = double(logical_indices{i});
end 
% % 
for i = 1:numel(Chain_Adjacency_matrices)
    longest_path_arr(i,1) = i;
    [longest_path_arr(i,2),longest_path_arr(i,3)] = count_unique_pathss(Chain_Adjacency_matrices{i});
disp(i);
disp(i/numel(Chain_Adjacency_matrices)*100)
end
longest_path_ = max(longest_path_arr(:,3));
%longest_path_arr_normalized = longest_path_arr(:,2)/longest_path_;
T = array2table(longest_path_arr);
column_names = {'Network Number', 'No. of paths (NumOfClusts)', 'Longest path size (MaxClustSize)'};
T.Properties.VariableNames = column_names;
chnNum = int2str(chainNum);
excel_file = strcat(Name,'_',traj,'C_',chnNum,'chains_clust_size_Num_1k_',temp,'K.xlsx');
writetable(T, excel_file);
end
end
end


% % Functions

% will return 4 things. 
% 1: no of connections, 
% 2: list of 1, 
% 3: number of components including isolated nodes and connections, 
% 4: List of 3.
function [paths,coll,saved_paths,conns,saved_conns] = count_unique_pathss(adjacency_matrix)
    num_nodes = size(adjacency_matrix, 1);
    node_list = 1:num_nodes;
    Links = {};
    % Compute unique paths using dynamic programming
    paths = 0;
    k = 1;
    paths_local = zeros(1,num_nodes);
    cancel = zeros(1,num_nodes);
      for k = 1:num_nodes
        nodes = 0;
        for i = 1:num_nodes
            if adjacency_matrix(k,i) == 1
                    nodes = nodes + 1;
                    Links{k}(1,nodes) = i;
            end
        end
      end
    Connections = {};
    
    for i = 1:numel(Links)
        Connections{i}(1,1) = i;
        [~,col] = size(Links{i});
        for j = 1:col
            Connections{i}(1,j+1)=Links{i}(1,j);
        end
    end


for ddd = 1:numel(Connections)
    for i = 1:numel(Connections)
        for j = i+1:numel(Connections)
        combined = intersect(Connections{i},Connections{j});
        if ~isempty(combined)
            Connections{i} = union(Connections{i},Connections{j});
        end
        end
    end
end

for i = 1:numel(Connections)
    for j = i+1:numel(Connections)
        combined = intersect(Connections{i},Connections{j});
        if ~isempty(combined)
            Connections{i} = union(Connections{i},Connections{j});
            Connections{j} = [];
        end
    end
end
coll = 0;
conns = 0;
for i = 1:numel(Connections)
    if ~isempty(Connections{i})
        conns = conns + 1;
        saved_conns{conns} = Connections{i};
        [~,col] = size(Connections{i});
        if col>1
            if col>coll
            coll = col;
            end
        paths = paths + 1;
        saved_paths{paths} = Connections{i};
        end       
    end

end
end
