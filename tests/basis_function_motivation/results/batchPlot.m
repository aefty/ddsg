
close all; clear;
dirFolder='';

folderStruc=dir('test*');
folderName={};

for i = 1:max(size(folderStruc))
    folderName{i} = folderStruc(i).name;
end

exact={};
gridpoints={};
gridpoints_write={};
sgInterpolate={};
hdmrInterpolate={};
x={};
w={};
c={};

pointsSG=[];
pointsHDMR=[];

errorSGavg=[];
errorHDMRavg=[];

errorSGmax=[];
errorHDMRmax=[];

errorSGavgL2=[];
errorHDMRavgL2=[];

for i = 1:max(size(folderName))
    exact{i} = csvread(strcat(dirFolder,folderName{i},'/exact.csv'));
    p_sg = csvread(strcat(dirFolder,folderName{i},'/gridpoints_sg.csv'));
    p_hdmr = csvread(strcat(dirFolder,folderName{i},'/gridpoints_hdmr.csv'));
    
    gridpoints{i} =[p_hdmr;p_sg];
    
    sgInterpolate{i} = csvread(strcat(dirFolder,folderName{i},'/sgInterpolate.csv'));
    hdmrInterpolate{i} = csvread(strcat(dirFolder,folderName{i},'/hdmrInterpolate.csv'));
    x{i} = csvread(strcat(dirFolder,folderName{i},'/x.csv'));
    w{i} = csvread(strcat(dirFolder,folderName{i},'/w.csv'));
    c{i} = csvread(strcat(dirFolder,folderName{i},'/c.csv'));
    
    delSG   = abs(exact{i} - sgInterpolate{i});
    delHDMR = abs(exact{i} - hdmrInterpolate{i});

    pointsSG(i,:) = gridpoints{i}(2,:);
    pointsHDMR(i,:) = gridpoints{i}(1,:);

    errorSGavg(i,:)   = sum(delSG)./max(size(exact{i}));
    errorHDMRavg(i,:) = sum(delHDMR)./max(size(exact{i}));
    
    errorSGmax(i,:)  = max(delSG);
    errorHDMRmax(i,:)= max(delHDMR);
    
    errorSGmin(i,:)  = min(delSG);
    errorHDMRmin(i,:)= min(delHDMR);
    
    errorSGavgL2(i,:)  = sqrt(sum(delSG.^2))/max(size(delSG));
    errorHDMRavgL2(i,:)= sqrt(sum(delHDMR.^2))/max(size(delSG));
    
    
end

tit = {'F1-OSCILLATORY','F2-PRODUCT PEAK','F3-CORNER PEAK','F4-LOG GAUSSIAN','F5-CONTINUOUS','F6-MAX'};
%leg = {'HDMR Max Error','HDMR Mean Error','HDMR Min Error','SG Max Error','SG Avg Error','SG Min Error'};
leg = {'HDMR \infty Error','HDMR L2 Error','SG \infty Error','SG L2 Error'};

xlab={'Grid Points'};
ylab={'Relative Error'};

for f = 1:6
figure 
loglog(pointsHDMR(:,f),errorHDMRmax(:,f),'--or'); hold on;
loglog(pointsHDMR(:,f),errorHDMRavgL2(:,f),'-or'); hold on;
%loglog(pointsHDMR(:,f),errorHDMRmin(:,f),'-.or'); hold on;

loglog(pointsSG(:,f),errorSGmax(:,f),'--*b'); hold on;
loglog(pointsSG(:,f),errorSGavgL2(:,f),'-*b'); hold on;
%loglog(pointsSG(:,f),errorSGmin(:,f),'-.*b');

legend(leg);
title(tit{f});
xlabel(xlab{1});
ylabel(ylab{1});
grid minor

end