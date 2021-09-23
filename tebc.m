function [EBC,BC,mark] = tebc(M,mask,type)
%% Calculate the node and edges betweenness centrality for targeted inter-node communications. 
%---------------------------------------------------------------------------------------------------------------------------------------------------%
% - Z.K.X. 2021/07/03
%---------------------------------------------------------------------------------------------------------------------------------------------------%
%% Input
%  (1) M: Binary/weighted connection matrix.
%  (2) mask: 0/1 matrix, that indicates the targeted paths between pair of nodes. 
%  (3) type: 1 - weighted matrix (default)/ 0 - distance matrix 
%% Output
%  (1) EBC: Edge betweenness centrality matrix.
%  (2) BC: Nodal betweenness centrality vector.
%  (3) mark
%      [1] mark.mark_edge: Path identification for edge.
%      [2] mark.mark_node1: Path identification for node (with direct path).
%      [3] mark.mark_node2: Path identification for node (without direct path).
%---------------------------------------------------------------------------------------------------------------------------------------------------%
%%
if nargin < 2
    mask = ones(length(M));
end
if nargin < 3
    type = 1;
end

if length(unique(M)) > 2
    if type == 1
        E = find(M); M(E)=1./M(E);
    end
end

EBC = zeros(length(M),length(M));
BC = zeros(length(M),1);

mark.mark_edge = repmat({[]},length(M),length(M));
mark.mark_node1 = repmat({[]},length(M),1);
mark.mark_node2 = repmat({[]},length(M),1);

for i = 1:length(M) - 1
    M(i,i) = 0;  
    if mean(M(i,:)) ~= 0
        [d,spath] = BellmanFord(M,i);        
        for j = i+1:length(M) 
            if mask(i,j) ~= 0 & d(j) ~= 0                
                pt = spath{j};                           
                for s = 1:length(pt)
                    path = pt{s};
                    for k = 1:length(path) - 1
                        EBC(path(k),path(k+1)) = EBC(path(k),path(k+1)) + 1/length(pt);
                        mark.mark_edge{path(k),path(k+1)} = [mark.mark_edge{path(k),path(k+1)};[i,j]];  
                        mark.mark_edge{path(k+1),path(k)} = [mark.mark_edge{path(k+1),path(k)};[i,j]];  
                        mark.mark_node1{path(k)} = [mark.mark_node1{path(k)};[i,j]];
                        if length(path) > 2 & k >= 2  
                            BC(path(k)) = BC(path(k)) + 1/length(pt);        
                            mark.mark_node2{path(k)} = [mark.mark_node2{path(k)};[i,j]];
                        end
                    end                                                
                end                                      
           end        
        end                
    end
end

EBC = tril(EBC) + flipud(rot90(triu(EBC)));
EBC = flipud(rot90(tril(EBC))) + EBC;
BC = BC*2;
end

%% Author: Bruno Luong
function [d, spath] = BellmanFord(arg1, arg2)
 
if isa(arg1,'graph') || isa(arg1,'digraph')
    G = arg1;
    if ismultigraph(G)
        edges = table2array(G.Edges);
        ij = edges(:,1:2);
        w = edges(:,3);
        n = max(ij,[],'all');
        A = accumarray(ij, w, [n,n], @min, 0, true);
        if isa(arg1,'graph')
            A = triu(A,1).'+A;
        end
    else
        A = G.adjacency('weight');
    end
else
    A = arg1;
end
if nargin < 2
    start = 1; % starting node
else
    start = arg2;
end
A = A.';
n = size(A,1);
start = max(min(start,n),1);
% initialize distance matrix
d = inf(n,1);
du = 0;
i = start;
if nargout < 2
    % only distance is requested 
    while true
        d(i) = du;
        [i, j, dij] = find(A(:,i));
        if isempty(i)
            break
        end
        [du, p] = sort(du(j)+dij);
        [i, p]  = sort(i(p)); % here we requires stable sorting, which is the case with MATLAB
        b = [true; diff(i,1,1)>0];
        i = i(b);
        du = du(p(b));
        b = du < d(i);
        i = i(b);
        du = du(b);
    end
else
    % shortest paths requested
    prev = cell(n,1);
    neig = inf(1,n);
    for c = 0:n-1
        d(i) = du;
        neig(i) = c;
        jprev = i;
        [i, j, dij] = find(A(:,jprev));
        if isempty(i)
            break
        end
        I = [i, du(j)+dij, jprev(j)];
        I = sortrows(I, [1 2]);
        dI = diff([0 0; I(:,1:2)],1,1) ~= 0;
        change = find(dI(:,1)|dI(:,2));
        Ic = I(change,:);
        b = dI(change,1) & (Ic(:,2) <= d(Ic(:,1)));
        lgt = diff([change; size(dI,1)+1],1,1);
        I = I(repelem(b, lgt),:);
        if isempty(I)
            break
        end
        lgtk = lgt(b);
        istart = cumsum([1; lgtk]);
        istop = istart(2:end)-1;
        istart(end) = [];
        % keep track the previous visit nodes
        % doing like jprev = mat2cell(I(:,3), lgtk, 1); without overhead
        jprev = cell(size(istart));
        j = I(:,3);
        for k=1:size(jprev)
            jprev{k} = j(istart(k):istop(k));
        end
        i = I(istart,1);
        du = I(istart,2);
        neq = du < d(i);
        if all(neq) % optimization by not bother with CELLFUN
            prev(i) = jprev;
        else
            prev(i(neq)) = jprev(neq);
            eq = ~neq;
            ieq = i(eq);
            prev(ieq) = cellfun(@union, prev(ieq), jprev(eq), 'unif', 0);
        end
    end
    % Second past to build shortest paths
    spath = cell(size(prev));
    spath{start} = {start};
    [neig, k] = sort(neig);
    k = k(2:find(isfinite(neig),1,'last'));
    for n = k
        spath{n} = cellfun(@(p) [p,n], cat(1, spath{prev{n}}), 'unif', 0);
    end % for-loop
end  
end % BellmanFord