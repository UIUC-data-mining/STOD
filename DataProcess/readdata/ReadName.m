function name = ReadName(filename)
f = fopen(filename);
name = textscan(f,'%[^\n]');
fclose(f);
for i=1:length(name{1})    
    if ~isempty(name{1}{i}) && name{1}{i}(end)<32
        name{1}{i}(end) = [];
    end
end
