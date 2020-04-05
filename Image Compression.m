% load image
close all; 
clear all; 
clc
load mandrill

k=1e2;
y = double(imread('mandrill.tiff'));

% partition image to M*M blocks and loop through blocks
M = 2; 
n = prod(size(y))/(3*M*M); 
d = size(y,1);
c = 0;
x=zeros(n,3*M*M);
for i=1:M:d
    for j=1:M:d
        c = c+1;
        % reshape each block to vector of dimension 3*M^2
        x(c,:) = reshape(y(i:i+M-1,j:j+M-1,:),[1,M*M*3]);
    end
end

% initialize centroid location
dprime = size(x,2);
rng(0);
perm = randperm(n);
m = x(perm(1:k),:); 

map_old = zeros(n,1);
map = ones(n,1);

% update centroid location
iter = 1;
while norm(map - map_old)>0
    iter
    map_old = map;
    
    dist = dist2(x,m);
    [Min,map] = min(dist,[],2);
    
    w(iter) = 0;
    for l=1:k
        a = map-l;
        zeroInd = find(a==0);
        m(l,:) = mean(x(zeroInd,:));
        w(iter) = w(iter)+sum(sum((x(zeroInd,:)-repmat(m(l,:),[size(x(zeroInd,:),1) 1])).^2,2));
    end
    iter = iter+1;
end

% construct compressed image
for l=1:k
    a = map-l;
    zeroInd = find(a==0);
    x1(zeroInd,:) = repmat(m(l,:),[numel(zeroInd) 1]);
end

c=0;
for i=1:M:d
    for j=1:M:d
        c = c+1;
        y1(i:i+M-1,j:j+M-1,:) = reshape(x1(c,:),[M M 3]);
    end
end

% objective function
figure
plot(w)
xlabel('Iterations');
ylabel('W(c)')

% compare compressed & original
figure; 
subplot(121)
imagesc(y1/256)
title('Compressed')
subplot(122)
imagesc(y/256)
title('Original')

% difference image
figure;
dy=y1-y;
imagesc(dy)