function gridconvert(fin,bot)

% display a hbrook format mesh 
%
% usage:  show_hbrook_mesh(fin,dtype,ftype)
%
% inputs:
% fin  (req): mesh file (string) 
% dtype(opt): display type
%             = 0 show everything 
%             = 1 show mesh only 
% ftype(opt): format type 
%             = 'new' no strange floats for boundary nodes and edges
%             = 'old' two strange floats after boundary nodes and edges
%
% G. Cowles, SMAST

% % process input and set defaults
%   if(nargin < 1)
%     msg = 'error using show_hbrook_mesh: must supply at least the mesh';
%     error(msg);
%   end;
% 
%   if isstr(varargin{1})
%     fin = varargin{1};
%   else
%     msg = 'error using show_hbrook_mesh: must supply at least the mesh';
%     error(msg);
%   end;
% 
% % check for dtype argument and process
%   if(nargin > 1) & (~isstr(varargin{2}))
%     dtype = varargin{2};
%   else
%     dtype = 0;
%   end;
% 
% % check for ftype argument and process 
%   if isstr(varargin{end})
%     ftype = varargin{end};
%   else
%     ftype = 'new';
%   end;
% 
%   clf;

dtype = 0; ftype = 'new';
  cnt = 0;

% read in header info for number of vertices, edges, and tris
  [dum,nvrtx,dum,nedge,dum,ntri] = textread(fin,'%s %d %s %d %s %d',1);
  cnt = cnt + 1;

% read vertex info
  [dum,x,y]   = textread(fin,'%s %f %f',nvrtx,'headerlines',cnt);
  coords      = zeros(nvrtx,2);
  coords(:,1) = x;
  coords(:,2) = y;
  clear x;
  clear y;
  cnt         = cnt + nvrtx;

% read edge info and shift to start from 1
  [dum,e1,e2] = textread(fin,'%s %d %d',nedge,'headerlines',cnt);
  enode       = zeros(nedge,2);
  enode(:,1)  = e1 + 1;
  enode(:,2)  = e2 + 1;
  cnt         = cnt + nedge;

% read triangle info and shift to start from 1
  [dum,n1,n2,n3] = textread(fin,'%s %d %d %d',ntri,'headerlines',cnt);
  nv             = zeros(ntri,3);
  nv(:,1)        = n1 + 1;
  nv(:,2)        = n2 + 1;
  nv(:,3)        = n3 + 1;
  cnt            = cnt + ntri;

% read number of side boundaries
  [dum,nsbd] = textread(fin,'%s %d',1,'headerlines',cnt);
  cnt        = cnt + 1;
  if(nsbd > 0)
  sbndry_cnt = zeros(nsbd,1);
  sbndry_lst = zeros(nedge,nsbd);
  sbndry_id  = zeros(nsbd,1);

% loop through side boundary types and read side boundary lists
  for i=1:nsbd

     % read side boundary id 
     [dum,sbndry_id(i)] = textread(fin,'%s %d',1,'headerlines',cnt);
     cnt                = cnt + 1;

     % read side boundary count for this type
     [dum,sbndry_cnt(i)] = textread(fin,'%s %d',1,'headerlines',cnt);
     cnt                = cnt + 1;

     % read side boundary for this type and shift to start from 1
     if(ftype == 0)
     [dum,tmp] = textread(fin,'%s %d',sbndry_cnt(i),'headerlines',cnt);
     else
     [dum,tmp,jnk,jnk] = textread(fin,'%s %d %f %f',sbndry_cnt(i),'headerlines',cnt);
     end;
     cnt                = cnt + sbndry_cnt(i);
     sbndry_lst(1:sbndry_cnt(i),i) = tmp(1:sbndry_cnt(i)) + 1;

  end;
  end; %nsbd > 0


% read number of boundary vertices
  [dum,nvbd] = textread(fin,'%s %d',1,'headerlines',cnt);
  cnt        = cnt + 1;
  if(nvbd > 0)
  vbndry_cnt = zeros(nvbd,1);
  vbndry_lst = zeros(nvrtx,nvbd);
  vbndry_id  = zeros(nvbd,1);

% loop through vertex boundary types and read vertex boundary lists
  for i=1:nvbd

     % read vertex boundary id 
     [dum,vbndry_id(1)] = textread(fin,'%s %d',1,'headerlines',cnt);
     cnt                = cnt + 1;

     % % read vertex boundary count for this type
     % [dum,vbndry_cnt(1)] = textread(fin,'%s %d',1,'headerlines',cnt);
     % cnt                = cnt + 1;
     vbndry_cnt(1)=1;

     % read vertex boundary for this type and shift to start from 1
     if(ftype == 0)
     [dum,tmp] = textread(fin,'%s %d',vbndry_cnt(1),'headerlines',cnt);
     else
     [dum,tmp,jnk,jnk] = textread(fin,'%s %d %f %f',vbndry_cnt(1),'headerlines',cnt);
     end;
     cnt                = cnt + vbndry_cnt(1);
     vbndry_lst(1:vbndry_cnt(1),i) = tmp(1:vbndry_cnt(1)) + 1;

  end;
  end; %nvbd > 0

% convert coordinates to polar

tol = 1e-14;  % tolerance to consider two points coinciding

r_eps= exp(bot)/(1-exp(bot));

r_phys = (1+ r_eps)*exp(coords(:,2)) - r_eps;
x_phys = r_phys.*cos(-coords(:,1));
y_phys = r_phys.*sin(-coords(:,1));
coords_phys = [x_phys y_phys];

keep = true(size(coords_phys,1),1); 

for i = 1:(size(coords_phys,1) -1)
    for j = (i+1):size(coords_phys,1)
        dist = norm(coords_phys(i,:) - coords_phys(j,:));
        if dist < tol
            keep(j) = false;
        end
    end
end

coords_clean = coords_phys(keep,:);

% find the mapping between old and new nodes
map = ones(size(coords_clean,1),1);
count = 1;
for i=1:size(coords_phys,1)
    if coords_phys(i,:)==coords_clean(count,:)
        map(i)=count;
        count = count+1;
    end
end

mapped_edges = map(enode);
edges_clean = mapped_edges(mapped_edges(:,1) ~= mapped_edges(:,2), :);

edges_sorted = sort(edges_clean, 2);
[edges_unique, ia] = unique(edges_sorted, 'rows','stable');
edges_clean2 = edges_clean(ia,:);

edge_map = NaN(size(enode,1),1);
count = 1;
% find the mapping between old and new edges
for i=1:size(mapped_edges,1)
    if mapped_edges(i,:)==edges_clean2(count,:)
        edge_map(i)=count;
        count = count+1;
    end
end


% update edge boundaries based on new numbers
count = 1;
for i = 1:length(sbndry_id)
    tmp_1 = edge_map(sbndry_lst(1:sbndry_cnt(i),i));
    tmp = tmp_1(~isnan(tmp_1));
    if ~isempty(tmp)
        sbndry_lst_new(1:length(tmp),count) = tmp;
        count = count+1;      
    end
end

% Update triangles
T_mapped = map(nv);



% T_sorted = sort(T_mapped, 2);
% [T_unique, ia] = unique(T_sorted, 'rows','stable');

T_clean = T_mapped( ...
    (T_mapped(:,1) ~= T_mapped(:,2)) & ...
    (T_mapped(:,1) ~= T_mapped(:,3)) & ...
    (T_mapped(:,2) ~= T_mapped(:,3)), :);

% Output the grd file
filename = './xy_rstrt1_b0.grd';
fid = fopen(filename, 'w');

% --- 1. Write Header ---
fprintf(fid, 'npnt: %d nseg: %d ntri: %d\n', size(coords_clean,1), size(edges_clean2,1), size(T_clean,1));

% --- 2. Write Points ---
for i = 1:size(coords_clean,1)
    fprintf(fid, '%d: %.10e %.10e\n', i-1, coords_clean(i,1), coords_clean(i,2));
end

% --- 3. Write Edges (Segments) ---
for i = 1:size(edges_clean2,1)
    fprintf(fid, '%d: %d %d\n', i-1, edges_clean2(i,1)-1, edges_clean2(i,2)-1);
end

% --- 4. Write Triangles ---
for i = 1:size(T_clean,1)
    fprintf(fid, '%d: %d %d %d\n', i-1, T_clean(i,1)-1, T_clean(i,2)-1, T_clean(i,3)-1);
end

% --- 5. Write Boundary Edge Groups ---
fprintf(fid, 'nebd: %d\n', numel(sbndry_lst_new(1,:)));

count = 1;
for k = 1:numel(sbndry_lst_new(1,:))
    fprintf(fid, 'idnum: %d\n', count);
    count = count+1;
    fprintf(fid, 'number: %d\n', nnz(sbndry_lst_new(:,k)));
    for j = 1:nnz(sbndry_lst_new(:,k))
        fprintf(fid, '%d: %d\n', j-1, sbndry_lst_new(j,k)-1);
    end
end

% % --- 6. Write Boundary Points (Optional) ---
% if exist('boundary_points', 'var') && ~isempty(boundary_points)
%     fprintf(fid, 'nvbd: %d\n', numel(boundary_points));
%     for k = 1:numel(boundary_points)
%         fprintf(fid, 'idnum: %d\n', boundary_points(k).idnum);
%         fprintf(fid, 'point: %d\n', boundary_points(k).point);
%     end
% end

fclose(fid);



% test= edge_map(sbndry_lst(1:sbndry_cnt(1),1));
% test_1 = test(~isnan(test));
% for i = 1:size(sbndry_id)
%     for j = 1:sbndry_cnt(i)
%         count = 0;



% keep = true(size(coords,1),1); % start by keeping all rows
% map = (1:size(coords,1))';  % start by mapping each vertex to itself

% for i = 1:size(coords,1)
%     if ~keep(i)
%         continue; % skip if already marked for removal
%     end
%     for j = i+1:size(coords,1)
%         if keep(j)
%             dist = norm(coords(i,:) - coords(j,:));
%             if dist < tol
%                 keep(j) = false; % remove the duplicate/close one
%                 map(j) = i;
%             end
%         end
%     end
% end

%Now I need to reorder the points based on points that are removed

% 
% coords_clean = coords(keep,:);
% mapped_edges = map(enode);
% edges_sorted = sort(mapped_edges, 2); 
% [edges_unique, ia] = unique(edges_sorted, 'rows','stable');
% edges_clean = edges_unique(edges_unique(:,1) ~= edges_unique(:,2), :);
% 
% % Update triangles
% T_mapped = map(nv);
% % Remove degenerate triangles
% T_clean = T_mapped( ...
%     (T_mapped(:,1) ~= T_mapped(:,2)) & ...
%     (T_mapped(:,1) ~= T_mapped(:,3)) & ...
%     (T_mapped(:,2) ~= T_mapped(:,3)), :);
% 
% 
% edge_map = zeros(size(e1,1),1); % preallocate
% % enode_sorted = sort(enode, 2);
% 
% for i = 1:size(e1)
%     % Find where this old edge is in the new list
%     if mapped_edges(i,1) == mapped_edges(i,2)
%         edge_map(i) = NaN;
%     else
%         idx = find(ismember(i, ia));
%         if ~isempty(idx)
%             edge_map(i) = idx;
%         else
%             edge_map(i) = NaN; % mark missing edges (e.g., if it was removed)
%         end
%     end
% end

% plot grid
%   patch('Vertices',coords,'Faces',nv,...
%           'edgecolor','black','facecolor','white');
%   hold on
% 
% 
%   if(dtype == 0)
% % plot boundary edges (class 1 = red, class 2 = blue, class 3 = green) 
%   color = ['r-','b-','g-','r-','b-','g-','r-','b-','g-'];
%   for i=1:nsbd
%      for j=1:sbndry_cnt(i)
%         n1 = enode(sbndry_lst(j,i),1);
%         n2 = enode(sbndry_lst(j,i),2);
%         x(1) = coords(n1,1); 
%         x(2) = coords(n2,1); 
%         y(1) = coords(n1,2); 
%         y(2) = coords(n2,2); 
%         plot(x,y,color(i))
%      end;
%   end;
%   clear x;
%   clear y;
% 
% % plot boundary verts (class 1 = blue, class 2 = green, class 3 = red)
%   clear color;
%   color = ['b';'g';'r'];
%   for i=1:nvbd
%      for j=1:vbndry_cnt(i)
%         n1 = vbndry_lst(j,i);
%         x  = coords(n1,1);
%         y  = coords(n1,2);
%         plot(x,y,[color(i),'+'])
%      end;
%   end;
% 
%   end; %dtype == 0
% 
%   axis equal

end
