function [xm,ym,ys] = my_avg(x,y, bin_no, min_no)
if nargin < 4; min_no = 1; end

    bins = linspace(min(x(:)),max(x(:)),bin_no+1);
    
    xm = bins(1:end-1) + diff(bins(1:2))/2;
    ym = nan*zeros(1,length(bins)-1);
    ys = nan*zeros(1,length(bins)-1);
    for i = 1:length(bins)-1
        xids = find((x>=bins(i)).*(x<bins(i+1)));
        if length(xids) > min_no
        ym(i) = nanmean(y(xids));
        ys(i) = nanstd(y(xids));% / sqrt(length(xids));
        end
    end
end