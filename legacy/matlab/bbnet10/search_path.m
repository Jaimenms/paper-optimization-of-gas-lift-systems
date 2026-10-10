function path = search_path(M)

[N,S] = size(M);

path = [];
for j=1:S
    iy = j;
    path{j}=j;
    while ~isempty(iy)
        ix = find(M(:,iy)==1);
        iy = find(M(ix,:)==-1);
        path{j}=[path{j} iy];
    end
end   