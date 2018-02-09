clear all
clc
close all

gtCols = {'sub1','x1','y1','xy1','int1',...
    'sub2','x2','y2','xy2','int2',...
    'sub3','x3','y3','xy3','int3',...
    'sub4','x4','y4','xy4','int4'};

boCols = {'sub1','wlsEst1','wlsStd1','pCountWls1','fitEst1','fitStd1','pCountFit1',...
    'sub2','wlsEst2','wlsStd2','pCountWls2','fitEst2','fitStd2','pCountFit2',...
    'sub3','wlsEst3','wlsStd3','pCountWls3','fitEst3','fitStd3','pCountFit3',...
    'sub4','wlsEst4','wlsStd4','pCountWls4','fitEst4','fitStd4','pCountFit4'};

goCols = 'pCountgo';

load('smallNsimulation.mat');
GT = groundtruth(1:size(bo,1),:);
BO = bo;
GO = go;

%% Sort stuff
[GTs, idx] = sortrows([sum(GT(:,strncmp(gtCols, 'int', 3)),2), GT(:,strncmp(gtCols, 'xy', 2))], [1 2 3 4 5]);
BO = BO(idx,:);
GO = GO(idx,:);

%%
for i = 1:size(BO,1)
    [~, jdx] = sort(GTs(i,2:end));
    y = BO(i,strncmp(boCols, 'pCountWls', 9));
    pcount(i,:) = y;
end
avg = aggregate([GTs(:,1), mean(GTs(:,2:end), 2)], 1, 2);
groupavg = aggregate([GTs(:,1), GO], 1,2 );
groupstd = aggregate([GTs(:,1), GO], 1,2, @std,1 );
groupn = aggregate([GTs(:,1), GO], 1,2, @count, 1 );

%% 
subplot(1,2,1)
edges = -20:5:70;
int = GTs(:,2:end);
pc = pcount(:);
[dist, hidx] = histc(int(:), edges);
avgpc = [unique(hidx), edges(unique(hidx))', aggregate([hidx, pc], 1, 2, @mean, 1),...
    aggregate([hidx, pc], 1, 2, @std, 1)./sqrt(dist(dist~=0))];

cnt = 1; 
outavg = [];
for i = 1:numel(edges)
    if cnt <= size(avgpc, 1) && i == avgpc(cnt,1) 
        outavg = [outavg; avgpc(cnt,:)];
        cnt = cnt + 1;
    else
        outavg = [outavg; i, edges(i), nan, nan];
    end
end

%% Plot filled plot for the individual analysis
x = [outavg(:,2)', fliplr(outavg(:,2)')];
m = [outavg(:,3)', fliplr(outavg(:,3)')];
s = [2.58 * outavg(:,4)', - 2.58 * fliplr(outavg(:,4)')];
lm = numel(m);

noNanIdx = find(~isnan(m(1:lm./2)));      % Find the indexes in the first half which aren't nan
breakpoint = find(diff(noNanIdx) ~= 1);   % Find where the gap is

cidx = [(1:lm./2)', (lm:-1:(lm./2 +1))']; % Create vectors for fill
p1 = cidx(noNanIdx(1:breakpoint),:);      % Find indexes where the first half of m does not equal nan (up to first break)
p1 = [p1(:,1); flipud(p1(:,2))];
fill(x(p1), m(p1) + s(p1), [0 .5 .5], 'FaceAlpha', .5)
hold on
p2 = cidx(noNanIdx(breakpoint+1:end),:);  % Find indexes where the first half of m does not equal nan (after first break)
p2 = [p2(:,1); flipud(p2(:,2))]; 
h1 = fill(x(p2), m(p2) + s(p2), [0 .5 .5], 'FaceAlpha', .5);
plot(x(p1(1:numel(p1)/2)), m(p1(1:numel(p1)/2)), '--k')
plot(x(p2(1:numel(p2)/2)), m(p2(1:numel(p2)/2)), '--k')

%% Plot filled plot for the group analysis
x = [avg(:,2)', fliplr(avg(:,2)')];
m = [groupavg(:,2)', fliplr(groupavg(:,2)')];
s = [2.58 * (groupstd./sqrt(groupn))', -2.58 * fliplr((groupstd./sqrt(groupn))')];
h2 = fill(x, m + s, [.75 0 .5], 'FaceAlpha', .5);hold on
plot(x, m, '--k')
set(gca,'YLim', [0 1.05])

xlabel('\beta_{xy}(ms)')
ylabel('p(Significant Interaction)')
title('Weighted Least Squares')
legend([h1, h2], 'Individual', 'Group')

%%
for i = 1:size(BO,1)
    [~, jdx] = sort(GTs(i,2:end));
    y = BO(i,strncmp(boCols, 'pCountFit', 9));
    pcount(i,:) = y;
end
avg = aggregate([GTs(:,1), mean(GTs(:,2:end), 2)], 1, 2);
groupavg = aggregate([GTs(:,1), GO], 1,2 );
groupstd = aggregate([GTs(:,1), GO], 1,2, @std,1 );
groupn = aggregate([GTs(:,1), GO], 1,2, @count, 1 );

%% 
subplot(1,2,2)
edges = -20:5:70;
int = GTs(:,2:end);
pc = pcount(:);
[dist, hidx] = histc(int(:), edges);
avgpc = [unique(hidx), edges(unique(hidx))', aggregate([hidx, pc], 1, 2, @mean, 1),...
    aggregate([hidx, pc], 1, 2, @std, 1)./sqrt(dist(dist~=0))];

cnt = 1; 
outavg = [];
for i = 1:numel(edges)
    if cnt <= size(avgpc, 1) && i == avgpc(cnt,1) 
        outavg = [outavg; avgpc(cnt,:)];
        cnt = cnt + 1;
    else
        outavg = [outavg; i, edges(i), nan, nan];
    end
end


%% Plot filled plot for the individual analysis
x = [outavg(:,2)', fliplr(outavg(:,2)')];
m = [outavg(:,3)', fliplr(outavg(:,3)')];
s = [2.58 * outavg(:,4)', - 2.58 * fliplr(outavg(:,4)')];

lm = numel(m);

noNanIdx = find(~isnan(m(1:lm./2)));      % Find the indexes in the first half which aren't nan
breakpoint = find(diff(noNanIdx) ~= 1);   % Find where the gap is

cidx = [(1:lm./2)', (lm:-1:(lm./2 +1))']; % Create vectors for fill
p1 = cidx(noNanIdx(1:breakpoint),:);      % Find indexes where the first half of m does not equal nan (up to first break)
p1 = [p1(:,1); flipud(p1(:,2))];
fill(x(p1), m(p1) + s(p1), [0 .5 .5], 'FaceAlpha', .5)
hold on
p2 = cidx(noNanIdx(breakpoint+1:end),:);  % Find indexes where the first half of m does not equal nan (after first break)
p2 = [p2(:,1); flipud(p2(:,2))]; 
h1 = fill(x(p2), m(p2) + s(p2), [0 .5 .5], 'FaceAlpha', .5);
plot(x(p1(1:numel(p1)/2)), m(p1(1:numel(p1)/2)), '--k')
plot(x(p2(1:numel(p2)/2)), m(p2(1:numel(p2)/2)), '--k')

%% Plot filled plot for the group analysis
x = [avg(:,2)', fliplr(avg(:,2)')];
m = [groupavg(:,2)', fliplr(groupavg(:,2)')];
s = [2.58 * (groupstd./sqrt(groupn))', -2.58 * fliplr((groupstd./sqrt(groupn))')];
h2 = fill(x, m + s, [.75 0 .5], 'FaceAlpha', .5);
set(gca,'YLim', [0 1.05])
plot(x, m, '--k')

xlabel('\beta_{xy}(ms)')
ylabel('p(Significant Interaction)')
title('Maximum Likelihood')