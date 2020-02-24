close all; clear all; clc

rng(0);
n = 200;
K = 2; 
e = [.7 .3]; 
w = [-2 1]; 
b = [.5 -.5]; 
v = [.2 .1]; 
for i=1:n
    x(i) = rand;
    if rand < e(1);
        y(i) = w(1)*x(i) + b(1) + randn*sqrt(v(1));
    else
        y(i) = w(2)*x(i) + b(2) + randn*sqrt(v(2));
    end
end
figure
plot(x,y,'bo')
hold on
t=0:0.01:1;
plot(t,w(1)*t+b(1),'k')
plot(t,w(2)*t+b(2),'k')

ek = [.5 .5];
wk = [1 -1];
bk = [0 0];
sk = repmat(var(y),1,2);

etak = zeros(2,K);

rik =zeros(n,K);

cont = true;

stopThreshold = 10 ^ (-4);
stepCount = 0;


pk = zeros(n,K);

phi = zeros(n,1);

X = [ones(n,1) x(:)];


thetaik = zeros(n,K);

theta = 0;
thetaOld = 0;
thetaHistory = zeros(10000,1);

L = 0;
LOld = 0;
Li = zeros(n,1);
LHistory = zeros(10000,1);

firstFlag = 1;

while cont == true
    stepCount = stepCount + 1;
    
    % E-step
    for sIndex = 1 : n
        for kIndex = 1 : K
            gPdf = normpdf(y(sIndex), wk(kIndex) * x(sIndex) + bk(kIndex),...
            sqrt(sk(kIndex)));
            rik(sIndex, kIndex) = ek(kIndex) * gPdf;
        end
        if(sum(rik(sIndex,:)) ~= 0)
            rik(sIndex,:) = rik(sIndex,:) / sum(rik(sIndex,:));
        else
            keyboard
        end
    end
    
    % M-step
    for kIndex = 1 : K
        ek(kIndex) = sum(rik(:,kIndex))/n;
        ck = diag(rik(:,kIndex));
        etak(:,kIndex) = (X' * ck * X) \ X' * ck * y';
        bk(kIndex) = etak(1, kIndex);
        wk(kIndex) = etak(2, kIndex);
        for sIndex = 1 : n
            pk(sIndex,kIndex) = (y(sIndex) - wk(kIndex) * x(sIndex) - bk(kIndex)) ^ 2;
        end
        sk(kIndex) = sum(rik(:,kIndex) .* pk(:,kIndex)) / sum(rik(:,kIndex));
    end
    
    Li = zeros(n,1);
    for sIndex = 1 : n
        for kIndex = 1 : K
            Li(sIndex, 1) =Li(sIndex, 1) + ek(kIndex) * normpdf(y(sIndex),wk(kIndex) *...
            x(sIndex) + bk(kIndex), sqrt(sk(kIndex)));
        end
    end
    
    L = sum(log(Li(:,1)));
    
    if abs(LOld - L) < stopThreshold
        if firstFlag == 1
            LOld = L;
            LHistory(stepCount) = L;
            firstFlag =0;
        else
            LHistory(stepCount) = L;
            cont = false;
        end
        
    else
        LOld = L;
        LHistory(stepCount) = L;
        firstFlag =0;
    end
end

plot(t,wk(1)*t+bk(1),'r--')
plot(t,wk(2)*t+bk(2),'r--')


figure
plot(LHistory(1:stepCount))
xlabel('Iterations');
ylabel('Log liklihood')