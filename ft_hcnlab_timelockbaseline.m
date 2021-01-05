function data = ft_hcnlab_timelockbaseline(timeVec, data, baseline, baselinetype, dimord)

if ~exist('dimord','var') || isempty(dimord)
    dim_ord = 'chan_time';
end

baselineTimes = false(size(baseline,1),numel(timeVec));
for k = 1:size(baseline,1)
  baselineTimes(k,:) = (timeVec >= baseline(k,1) & timeVec <= baseline(k,2));
end

if length(size(data)) ~= length(strsplit(dim_ord,'_'))
  ft_error(['data matrix should have ' sum2str(length(strsplit(dim_ord,'_'))) ' dimensions '...
      ' for dimord = ' dimord]);
end

% compute mean of time/frequency quantity in the baseline interval,
% ignoring NaNs, and replicate this over time dimension
if length(size(data)) == 2
    if size(baselineTimes,1)==size(data,1)
        % do channel specific baseline
        meanVals = nan+zeros(size(data));
        stdVals = nan+zeros(size(data));
        for k = 1:size(baselineTimes,1)
            meanVals(k,:) = repmat(nanmean(data(k,baselineTimes(k,:)), 2), [1 size(data, 2)]);
            stdVals(k,:) = repmat(nanstd(data(k,baselineTimes(k,:)), 0, 2), [1 size(data, 2)]);
        end
    else
        meanVals = repmat(nanmean(data(:,baselineTimes), 2), [1 size(data, 2)]);
        stdVals = repmat(nanstd(data(:,baselineTimes), 0, 2), [1 size(data, 2)]);
    end
else % length(size(data)) == 3
    if size(baselineTimes,1)==size(data,2)
        % do channel specific baseline
        meanVals = nan+zeros(size(data));
        stdVals = nan+zeros(size(data));
        for k = 1:size(baselineTimes,1)
            meanVals(:,k,:) = repmat(nanmean(data(:,k,baselineTimes(k,:)), 3), [1 1 size(data, 3)]);
            stdVals(:,k,:) = repmat(nanstd(data(:,k,baselineTimes(k,:)), 0, 3), [1 1 size(data, 3)]);
        end
    else
        meanVals = repmat(nanmean(data(:,:,baselineTimes), 3), [1 1 size(data, 3)]);
        stdVals = repmat(nanstd(data(:,:,baselineTimes), 0, 3), [1 1 size(data, 3)]);
    end
end

if (strcmp(baselinetype, 'absolute'))
  data = data - meanVals;
elseif (strcmp(baselinetype, 'relative'))
  data = data ./ meanVals;
elseif (strcmp(baselinetype, 'relchange'))
  data = (data - meanVals) ./ meanVals;
elseif (strcmp(baselinetype, 'normchange')) || (strcmp(baselinetype, 'vssum'))
  data = (data - meanVals) ./ (data + meanVals);
elseif (strcmp(baselinetype, 'db'))
  data = 10*log10(data ./ meanVals);
elseif (strcmp(baselinetype,'zscore'))
    data=(data-meanVals)./stdVals;
elseif (strcmp(baselinetype, 'none'))
else
  ft_error('unsupported method for baseline normalization: %s', baselinetype);
end
