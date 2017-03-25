function [arrayOfLines] = fileImport(dataFile)

filename = fopen(dataFile);
line = fgetl(filename);
n=1;
arrayOfLines = cell(30,1);
while ischar(line)
    splitLine = strread(line, '%s');
    if (isempty(splitLine) == 0 && strcmp(splitLine(1), '#') == 0 && strcmp(splitLine(1), ' ') == 0)
        arrayOfLines{n} = splitLine;
        n = n + 1;
    end
    line = fgetl(filename);
end
fclose(filename);