close all; clear;
dirFolder='';

folder='test*';

folderStruc=dir(folder);
folderName={};

for i = 1:max(size(folderStruc))
    folderName{i} = folderStruc(i).name;
end

data={};
X={};
for i = 1:max(size(folderName))
    data{i} = csvread(strcat(dirFolder,folderName{i},'/basisFunction.csv'));
    X{i} = csvread(strcat(dirFolder,folderName{i},'/X.csv'));
end

% Number of Basis function or Test
N = max(size(data));


for i = 1:N
    
    % Number supports
    S = size(data{i},1);
    leg={};
    
    subplot(2,round(N/2+.4),i);
    for j=1:S
        plot(X{i},data{i}(j,:));
        hold on;
       leg{end+1}=sprintf('(%d)',j);
    end
    legend(leg);
    title(sprintf('Basis Function (%d)',i));
    xlabel('X');
    ylabel('Value');
    grid minor
end