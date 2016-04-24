%% Set Sampling Rate, Smoothing Window Size & Baseline and Response Period
% Sampling Rate
fs=1000;
% Smoothing Window (ms)
smval=100;
% Beginning of Baseline
baseBeg=-0.506;
% End of Plot Time
plotEnd=0.505;
plotTim=baseBeg:(1/fs):plotEnd;
% Required viewing time during 1st Presentation (ms)
viewThresh=500;
% Minimum time required for timelock.time
viewReq=1;
% All 45 of Mike's VPLT Recordings
filelist={  'IW0604173'; 'IW0604183'; 'IW0604212'; 'IW0606203'; 'IW0606212';...
    'IW0606223'; 'IW0606272'; 'IW0606283'; 'IW0606292'; 'IW0606302';...
    'IW0607052'; 'IW0607062'; 'IW0607072'; 'IW0607122'; 'IW0607133';...
    'IW0607172'; 'IW0607202'; 'IW0607212'; 'IW0607253'; 'IW0607262';...
    'IW0607272'; 'IW0607312'; 'IW0608022'; 'IW0608032'; 'IW0608042';...
    'MP0702015'; 'MP0702084'; 'MP0702093'; 'MP0702123'; 'MP0702133';...
    'MP0702224'; 'MP0702233'; 'MP0702273'; 'MP0704113'; 'MP0704132';...
    'MP0704162'; 'MP0704182'; 'MP0704192'; 'MP0704232'; 'MP0705022';...
    'MP0705102'; 'MP0705142'; 'MP0803202'; 'MP0803212'; 'MP0803252'};
%
facewaveforms = cell(1,size(filelist,1));
nonfacewaveforms = cell(1,size(filelist,1));
facewaveforms2 = cell(1,size(filelist,1));
nonfacewaveforms2 = cell(1,size(filelist,1));
for fileloop=19:20%:size(filelist,1)
    
    % Get Filename
    fid=filelist{fileloop}(1:9);
    % Get Trial Indices for Images of Faces and Non-Faces
    % and Name of Stimulus Set (setnum) and logical index of face trials (sel)
    [faces, nonfaces, setnum, sel] = evalface(fid(1:9));
    % Load timelock structure for 1st image presentation
    f1=load(['C:\Users\seth.koenig\Documents\MATLAB\Face Cells\' fid 'fratenov']);
    % Load timelock structure for 2nd image presentation
    f2=load(['C:\Users\seth.koenig\Documents\MATLAB\Face Cells\' fid 'fraterep']);
    
    resp1dat=squeeze(f1.timelock.trial(:,1,ft_nearest(f1.timelock.time,0):ft_nearest(f1.timelock.time,(viewThresh/1000))));
    
    % Find All Images Viewed for at Least viewThresh ms During Encoding
    stimView{1,1}=find((sum(~isnan(resp1dat),2)>=viewThresh));
    % All FACE Images Viewed for at Least viewThresh ms During Encoding
    stimView{1,2}=intersect(faces,stimView{1,1});
    % All NON-FACE Images Viewed for at Least viewThresh ms During Encoding
    stimView{1,3}=intersect(nonfaces,stimView{1,1});
    clear resp1dat
    
    for cellloop=1:size(f1.timelock.label,1)
        % Get Name of Cell
        cellID=[fid '.' f1.timelock.label{cellloop}(end-1:end)];
        % Index the cell
        cellind=strmatch(['sig00' f1.timelock.label{cellloop}(end-1:end)],f1.timelock.label);
        
        if f1.timelock.time(end)>=viewReq && f2.timelock.time(end)>=viewReq
            
            % 1st Presentation:
            % Faces
            d1=squeeze(f1.timelock.trial(stimView{1,2},cellind,ft_nearest(f1.timelock.time,baseBeg):ft_nearest(f1.timelock.time,plotEnd)));
            % Non-Faces
            d3=squeeze(f1.timelock.trial(stimView{1,3},cellind,ft_nearest(f1.timelock.time,baseBeg):ft_nearest(f1.timelock.time,plotEnd)));
            
            %             figure
            [d3smth, ~] = nandens(d3,smval,'gauss',fs,'convol');

            [d1smth, ~] = nandens(d1,smval,'gauss',fs,'convol');

            facewaveforms{fileloop} = d1smth;
            nonfacewaveforms{fileloop} = d3smth;
                     
%             figure;
%             hold on
%             h = area([0.123 0.194],[12 12]);
%             set(h,'FaceColor',[.75 .75 .75])
%             set(h,'EdgeColor','none')
%             
%             dofill(plotTim,d3,'b',1,smval,0,0,1,fs,1,1)
%             hold on;
%             dofill(plotTim,d1,'r',1,smval,0,0,1,fs,1,1)
%             ylabel('Spikes/s');
%             xlabel('Time From Stimulus Onset (s)');
%             titleLabel=[cellID,',PSTH Encoding Faces & Non-Faces'];
%             title(titleLabel);
%             n1=['Non-Faces, N = ',num2str(size(d3,1))];
%             n2=['Faces, N = ',num2str(size(d1,1))];
%             legend(n1,n2,'Location','Best')
%             xlim([baseBeg plotEnd])
%             line([0 0],get(gca,'YLim'),'LineWidth',1,'LineStyle','--','Color',[0 0 0]);
%             pause(1)
%             close
%             
            % 2nd Presentation:
            % Faces
            d2=squeeze(f2.timelock.trial(stimView{1,2},cellind,ft_nearest(f2.timelock.time,baseBeg):ft_nearest(f2.timelock.time,plotEnd)));
            % Non-Faces
            d4=squeeze(f2.timelock.trial(stimView{1,3},cellind,ft_nearest(f2.timelock.time,baseBeg):ft_nearest(f2.timelock.time,plotEnd)));
            
            %             figure
            [d4smth, ~] = nandens(d4,smval,'gauss',fs,'convol');
            %             plot(plotTim,d4smth,'b')
            %             hold on
            [d2smth, ~] = nandens(d2,smval,'gauss',fs,'convol');
            %             plot(plotTim,d2smth,'r')
            
            facewaveforms2{fileloop} = d2smth;
            nonfacewaveforms2{fileloop} = d4smth;
            
%                         figure;
%                         dofill(plotTim,d4,'b',1,smval,0,0,1,fs,1,1)
%                         hold on;
%                         dofill(plotTim,d2,'r',1,smval,0,0,1,fs,1,1)
%             ylabel('Spikes/s');
%             xlabel('Time From Stimulus Onset (s)');
%             titleLabel=[cellID,',PSTH Recognition Faces & Non-Faces'];
%             title(titleLabel);
%             n1=['Non-Faces, N = ',num2str(size(d4,1))];
%             n2=['Faces, N = ',num2str(size(d2,1))];
%             legend(n1,n2,'Location','Best')
%             xlim([baseBeg plotEnd])
%             line([0 0],get(gca,'YLim'),'LineWidth',1,'LineStyle','--','Color',[0 0 0]);
%             pause(1)
%             close
            
        end
    end
end
%%
clearvars -except nonfacewaveforms  nonfacewaveforms2 facewaveforms facewaveforms2 plotTim
%
% allwaveforms = [];
% for i = 1:length(nonfacewaveforms)
%     allwaveforms = [allwaveforms;nonfacewaveforms{i}-nanmin(nonfacewaveforms{i})];
% end
%
% for i = 1:length(nonfacewaveforms2)
%     allwaveforms = [allwaveforms;nonfacewaveforms2{i}-nanmin(nonfacewaveforms2{i})];
% end
%
% for i = 1:length(facewaveforms)
%     allwaveforms = [allwaveforms;facewaveforms{i}-nanmin(facewaveforms{i})];
% end
%
% for i = 1:length(facewaveforms2)
%     allwaveforms = [allwaveforms;facewaveforms2{i}-nanmin(facewaveforms2{i})];
% end
% for i = 1:size(allwaveforms,1);
%     allwaveforms(i,isnan(allwaveforms(i,:))) = nanmean(allwaveforms(i,:));
% end
allwaveforms = [];
for i = 1:size(facewaveforms,2)
    temp1 = facewaveforms{i};
    temp1(isnan(temp1)) = nanmean(temp1);
    temp2 = nonfacewaveforms{i};
    temp2(isnan(temp2)) = nanmean(temp2);
    allwaveforms = [allwaveforms;temp1-temp2];
end
for i = 1:size(facewaveforms2,2)
    temp1 = facewaveforms2{i};
    temp1(isnan(temp1)) = nanmean(temp1);
    temp2 = nonfacewaveforms2{i};
    temp2(isnan(temp2)) = nanmean(temp2);
    allwaveforms = [allwaveforms;temp1-temp2];
end

%%
[wcoeff,score,latent,tsquared,explained] = pca(allwaveforms,'VariableWeights','variance');
%%
figure
hold all
for i = 1:size(score,1)
    plot(score(i,1),score(i,2),'.')
end
%%
sil = zeros(1,10); %determines the number of clusters by comparing the ratio
%of intercluster and intracluster distances, faster mod of silhouette
for numclusts = 4:10
    T = kmeans(score(:,2:7),numclusts,'replicate',5);
    [silh] = InterVSIntraDist(score(:,2:7),T);
    sil(numclusts) = mean(silh);
end
sil(sil > 0.95*max(sil)) = 1;
numclusters = find(sil == max(sil));
idx = kmeans([score(:,1) score(:,2) score(:,3) score(:,4) score(:,5) score(:,6) score(:,7)],floor(mean(numclusters)),'replicate',5);

clr = 'rgbkmyrgb';
figure
hold on
for i = 1:length(idx)
    if ~isnan(idx(i))
        plot3(score(i,1),score(i,2),score(i,3),['.' clr(idx(i))])
    end
end
%%
count = zeros(1,max(idx));
avgwaveform = zeros(max(idx),size(allwaveforms,2));
for i = 1:length(idx)
    avgwaveform(idx(i),:) = avgwaveform(idx(i),:) + allwaveforms(i,:);
    count(idx(i)) = count(idx(i))+1;
end

figure
hold all
for i = 1:max(idx)
    plot(plotTim,avgwaveform(i,:)/count(i));
end
legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4')
%%
ind = find(idx ==6);
figure
hold all
for i = 1:length(ind)
    plot(plotTim,allwaveforms(ind(i),:))
end
%%
figure
hold all
for i = 38%:length(nonfacewaveforms);
    if ~isempty(facewaveforms2{i})
        plot(plotTim,facewaveforms2{i})
    end
    %     pause(0.5)
end
