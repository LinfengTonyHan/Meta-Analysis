function ES = metaCorr(studyList, dataSheet, corrType)
%% This function is used for a meta-analysis of correlations (i.e., Pearson's r). It computes the overall effect sizes and generates the forest plot.
%
%% Developed by Linfeng Han (he/him). For questions, please contact Linfeng at hanlf@sas.upenn.edu
%
%% studyList:
% A list of studies incorporated in the meta-analysis, the data should be formatted in a cell array.
% Example: {'Epstein & Kanwisher, 1998','Epstein et al., 2017','Julian et al., 2018'}
%
%% dataSheet:
% An M by 2  numerical matrix, where M equals to the number of studies included. It includes the correlation coefficients
% and sample sizes of all studies (two necessary parameters used in the
% meta-analysis). Each row should consist of one correlation coefficient and one sample size.
% Example: [0.42, 302; 0.24, 120; 0.55, 59]
%
%% corrType:
% It should be 'r' or 'z', where 'r' refers to Pearson's r, and 'z' refers to Fisher's z;
% The default corrType is Pearson's r.

if nargin == 2
    corrType = 'r';
end

numStudy = length(studyList); %number of studies
numObs = size(dataSheet, 1); %number of observations

if numStudy ~= numObs
    err('Invalid input. The number of studies does not match the number of observations.');
end

Ns = numObs;

%% forest plot parameters
digits = 2;
textLoc.CI_X = 1.3;
textLoc.negCorr_X = -0.25;
textLoc.posCorr_X = 0.25;
textLoc.label_Y = 1.75;
textLoc.weight_X = 1;

%% Process the data: transforming from Pearson's r (Pr) to Fisher's z (Fz) and averaging across studies
rawCorrCoefs = dataSheet(:, 1);
sampleSize = dataSheet(:, 2);

if strcmp(corrType, 'r')
    Fz = 0.5 * (log(1 + rawCorrCoefs) - log(1 - rawCorrCoefs));
    Pr = rawCorrCoefs;
elseif strcmp(corrType, 'z')
    Fz = rawCorrCoefs;
    Pr = (exp(2 * Fz) - 1)/(exp(2 * Fz) + 1); %transforming Fisher's z back to Pearson's r
end
FzSEM = 1./sqrt(sampleSize - 3);
CorrSEM = (1 - Pr.^2) .* FzSEM;

Range = icdf('normal',[0.025 0.975], 0, 1);
lowerRange = Range(1);
upperRange = Range(2);

lowerLim = Pr + lowerRange * CorrSEM;
upperLim = Pr + upperRange * CorrSEM;

%% Calculate the overall effect sizes and adding them to the first row
weight = sampleSize * 100/sum(sampleSize);
weightedFz = sum(Fz .* sampleSize)/sum(sampleSize);
weightedPr = (exp(2 * weightedFz) - 1)/(exp(2 * weightedFz) + 1);
weightedFzSEM = 1/sqrt(sum(sampleSize) - 3);
weightedCorrSEM = (1- weightedPr^2) * weightedFzSEM;
Fz(2:(end + 1)) = Fz; Fz(1) = weightedFz;
Pr(2:(end + 1)) = Pr; Pr(1) = weightedPr;
FzSEM(2:(end + 1)) = FzSEM; FzSEM(1) = weightedFzSEM;
CorrSEM(2:(end + 1)) = CorrSEM; CorrSEM(1) = weightedCorrSEM;
sampleSize(2:(end + 1)) = sampleSize; sampleSize(1) = sum(sampleSize);
lowerLim(2:(end + 1)) = lowerLim; lowerLim(1) = weightedPr + lowerRange * weightedCorrSEM;
upperLim(2:(end + 1)) = upperLim; upperLim(1) = weightedPr + upperRange * weightedCorrSEM;
weight(2:(end + 1)) = weight; weight(1) = 100;

%% Formatting of ES
%Pearson's r, %Sample Size, %Fisher's z, %Fisher's z-SEM, %Corr-SEM, %Lower Limit, %Upper Limit
ES = [Pr, sampleSize, Fz, FzSEM, CorrSEM, lowerLim, upperLim];
Pr = roundn(Pr, -digits); Fz = roundn(Fz, -digits); FzSEM = roundn(FzSEM, -digits);
CorrSEM = roundn(CorrSEM, -digits); lowerLim = roundn(lowerLim, -digits); upperLim = roundn(upperLim, -digits);
weight = roundn(weight, -digits);

figure;
hold on
xlim([-1, 2]);
ylim([0, Ns + 2]);
Fplot = plot([0, 0], [0, Ns + 1.4], '--k', 'LineWidth', 1.8);
FplotAxis = ancestor(Fplot, 'axes');
xrule = FplotAxis.XAxis;
yrule = FplotAxis.YAxis;
xticks([-1:0.2:1, 2]);
xticklabels({'-1', '-0.8', '-0.6', '-0.4', '-0.2', '0', '0.2', '0.4', '0.6', '0.8', '1', 'Correlation Coefficient'});
xrule.FontSize = 15;
xrule.FontName = 'Palatino';

yticks(1:(Ns+2));
yticklabels(['Overall Effect Size', studyList, 'Studies']);
yrule.FontSize = 13;
yrule.FontName = 'Times';

for ite = 1:(Ns + 1)
    textContent.weightValue{ite} = num2str(weight(ite));
    text(textLoc.weight_X, ite, textContent.weightValue{ite}, 'HorizontalAlignment', 'left', 'fontsize', 14);
    
    textContent.CI{ite} = [num2str(Pr(ite)), ' [', num2str(lowerLim(ite)), ' ', num2str(upperLim(ite)), ']'];
    text(textLoc.CI_X, ite, textContent.CI{ite}, 'HorizontalAlignment', 'left', 'fontsize', 14);
    if ite > 1
        plot([lowerLim(ite), upperLim(ite)], [ite, ite], '-k', 'LineWidth',1.8);
        plot(Pr(ite), ite, 's','markerSize', sqrt(weight(ite)), 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'); %Based on weight
    elseif ite == 1
        plot([lowerLim(ite), upperLim(ite)], [ite, ite], '-b', 'LineWidth',3);
        plot(Pr(ite), ite, 'd','markerSize', 15, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
    end
end

textContent.negCorr = 'Negatively Correlated';
textContent.posCorr = 'Positively Correlated';
textContent.confInt = 'Mean [95% Confidence Interval]';
textContent.weight = 'Weight';
text(textLoc.negCorr_X, Ns + textLoc.label_Y, textContent.negCorr, 'HorizontalAlignment', 'right', 'fontsize', 14);
text(textLoc.posCorr_X, Ns + textLoc.label_Y, textContent.posCorr, 'HorizontalAlignment', 'left', 'fontsize', 14);
text(textLoc.CI_X, Ns + textLoc.label_Y, textContent.confInt, 'HorizontalAlignment', 'left', 'fontsize', 14);
text(textLoc.weight_X, Ns + textLoc.label_Y, textContent.weight, 'HorizontalAlignment', 'left', 'fontsize', 14);
end
