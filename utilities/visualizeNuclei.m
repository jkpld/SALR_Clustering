function visualizeNuclei(markers)

if nargin < 1
    markers = cell(1,5);
end

imgs = cell(3,3);
truth = cell(3,3);
mk = cell(3,3);
% 67
idx = randperm(484,2);
[imgs{1,1}, truth{1,1}, mk{1,1}] = getTestObject('67',idx(1),markers{5});
[imgs{3,3}, truth{3,3}, mk{3,3}] = getTestObject('67',idx(2),markers{5});

idx = randperm(484,2);
[imgs{1,2}, truth{1,2}, mk{1,2}] = getTestObject('5',idx(1),markers{4});
[imgs{3,2}, truth{3,2}, mk{3,2}] = getTestObject('5',idx(2),markers{4});

idx = randperm(484,2);
[imgs{2,1}, truth{2,1}, mk{2,1}] = getTestObject('4',idx(1),markers{3});
[imgs{2,3}, truth{2,3}, mk{2,3}] = getTestObject('4',idx(2),markers{3});

idx = randperm(484,2);
[imgs{1,3}, truth{1,3}, mk{1,3}] = getTestObject('3',idx(1),markers{2});
[imgs{3,1}, truth{3,1}, mk{3,1}] = getTestObject('3',idx(2),markers{2});

idx = randperm(484,1);
[imgs{2,2}, truth{2,2}, mk{2,2}] = getTestObject('2',idx(1),markers{1});

maxI = cellfun(@(x) max(x(:)), imgs);
minI = cellfun(@(x) min(x(:)), imgs);

maxI = double(max(maxI(:)));
minI = double(max(minI(:)));

maxL = max(cellfun(@(x) max(size(x)), imgs(:)));


c_teal = [69, 178, 157]/255; % teal
c_darkblue = [51, 77, 92]/255; % dark blue
c_green = [161, 208, 68]/255; % green
c_yellow = [239, 201, 76]/255; % yellow
c_orange = [226, 122, 63]/255; % orange
c_red = [223, 73, 73]/255; % red
c_gray = [164, 164, 164]/255; % gray


fig = figure;
fig.Units = 'inch';
fig.Position(3:4) = [5,3]*0.75;
% axbg = axes('Position',[0 0 1 1],'Color',0.3*[1 1 1]);
% axbg.XAxis.Visible = 'off';
% axbg.YAxis.Visible = 'off';
for i = 3:-1:1
    for j = 3:-1:1
        pos = [(i-1)/3, (j-1)/3, 1/3, 1/3];
        ax(i,j) = axes('Position',pos);

        tmpIm = imgs{i,j};
        tmpTr = truth{i,j};
        tmpMk = mk{i,j};
        if size(tmpIm,2)<size(tmpIm,1)
            tmpIm = tmpIm';
            tmpTr = fliplr(tmpTr);
            tmpMk = fliplr(tmpMk);
        end

        tmpIm = uint16(double(intmax('uint16'))*(double(tmpIm)-minI)/(maxI-minI));
        imshow(tmpIm,'Parent',ax(i,j));

        line(tmpTr(:,2),tmpTr(:,1),'marker','o','MarkerSize',4,'MarkerFaceColor',c_red,'MarkerEdgeColor','k','LineStyle','none','Parent',ax(i,j))
        if ~isempty(tmpMk)
            line(tmpMk(:,2),tmpMk(:,1),'marker','o','MarkerSize',4,'MarkerFaceColor',c_green,'MarkerEdgeColor','k','LineStyle','none','Parent',ax(i,j))
        end
        if size(tmpIm,1)<size(tmpIm,2)
            ax(i,j).XLim = ax(i,j).XLim + [-1,1]*(maxL-ax(i,j).XLim(2))/2;
        end
    end
end



set(ax,'Visible','off')
% setTheme(fig,'dark')
fig.Color = 0.3*[1 1 1];

end



function [im, truth, markers] = getTestObject(n,num, markers)


load_img = @(n) imread(['exampleImages\testImage_image_LD' n 'P24.tif']);
truth = load_truth(n);
I = load_img(n);

objSize = round(size(I)/22);

[i,j] = ind2sub([22,22],num);
offset = [(i-1)*objSize(1), (j-1)*objSize(2)];

im = I( (1:objSize(1)) + offset(1), (1:objSize(2)) + offset(2) );

mask = im==0;
rowRemove = all(mask,2);
colRemove = all(mask,1);

offset = offset + [find(~rowRemove,1,'first'), find(~colRemove,1,'first')] - 1;
im(rowRemove,:) = [];
im(:,colRemove) = [];

truth = truth(truth(:,1)==num,2:3) - offset;

if ~isempty(markers)
    markers = markers(markers(:,3)==num,1:2) - offset;
end
end
