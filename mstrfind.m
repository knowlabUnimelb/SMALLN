function idx = mstrfind(C, S)

idx = cellfun(@(s)find(strcmp(C, s)), S, 'UniformOutput', false);
idx(cellfun(@isempty, idx)) = [];
idx = cell2mat(idx);