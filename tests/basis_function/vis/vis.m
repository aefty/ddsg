

data={};

data.default=[
    65	3.73E-03
    257 2.27E-04
    1025	1.45E-05
    4097	9.30E-07
    16385	5.65E-08];


data.quad=[
    127	3.66E-04
    511	3.82E-06
    2047	1.20E-08
    8191	1.88E-10];


data.flipup=[
    127	1.08E-03
    127	5.62E-05
    2047	3.57E-06
    8191	2.25E-07];


semilogy(data.default(:,1),data.default(:,2),':'); hold on;
semilogy(data.quad(:,1),data.quad(:,2),'-'); hold on;
semilogy(data.flipup(:,1),data.flipup(:,2),'--');
legend({'Defualt','Quad','Flipup'});
