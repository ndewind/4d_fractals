% mandelbrot set scratch

res = [256,256,256,32];
% res = [16,16,16,16];
% center = [-0.2566,-0.7679];
center = [-0,-0];
% sz = [.00005,.00005];
sz = [2,2];
cr = linspace(center(1)-sz(1),center(1)+sz(1),res(1));
ci = linspace(center(2)-sz(2),center(2)+sz(2),res(2));
zSizeFactor = 1;
zr = linspace(-sz(1)*zSizeFactor,sz(1)*zSizeFactor,res(3));
zi = linspace(-sz(1)*zSizeFactor,sz(1)*zSizeFactor,res(4));
iterDepth = 2^4;
img = zeros(res);

for crindx = 1:numel(cr)
    fprintf('%05d/%05d\n',crindx,numel(cr))
    for ciindx = 1:numel(ci)
        for zrindx = 1:numel(zr)
            for ziindx = 1:numel(zi)
                zold = complex(zr(zrindx),zi(ziindx));
                c = complex(cr(crindx),ci(ciindx));
                %         zold = complex(zstartr,zstarti);
                niter = 0;
                iterStop = false;
                while niter < iterDepth & iterStop == false
                    niter = 1+niter;
                    znew = zold^2 + c;
                    img(ciindx,crindx,zrindx,ziindx) = niter;
                    if abs(znew) > 2
                        iterStop = true;
                    end
                    zold = znew;
                end
                
            end
        end
    end
end

load('4_cube_256')

customColorScale = repmat(linspace(0,1,iterDepth)',1,3);

% gamma correction
customColorScale = customColorScale.^(2.5/log2(iterDepth));

% colorDepth = 16;
% normedMinImgVal = min(img(:))/iterDepth;
% customColorScale = repmat(linspace(normedMinImgVal,1,colorDepth)',1,3);
% customColorScale = cat(1,customColorScale,ones(iterDepth - colorDepth,3));

customColorScale(end,:) = [0,0,0]; % make set pixels black for contrast

figure(1)
% for kzr = 1:numel(zr)
%     imagesc(cr,ci,img(:,:,kzr,kzr),[min(img(:)),max(img(:))])
%     set(gca,'dataaspectratio',[1,1,1])
%     colormap(customColorScale)
%     pause(0.0625)
% end
%

imgInSet = img == iterDepth;
imgSetList = zeros(sum(imgInSet(:)),4);
cpnt = 0;
for crindx = 1:numel(cr)
    fprintf('%05d/%05d\n',crindx,numel(cr))
    for ciindx = 1:numel(ci)
        for zrindx = 1:numel(zr)
            for ziindx = 1:numel(zi)
                if imgInSet(crindx,ciindx,zrindx,ziindx)
                    cpnt = cpnt + 1;
                    imgSetList(cpnt,:) = ...
                        [cr(crindx),ci(ciindx),zr(zrindx),zi(ziindx)];
                    % plot3(crindx,ciindx,zrindx,'.k');hold on;
                end
            end
        end
    end
end
%hold off;
% figure(1)
% for zrindx = 1:1:numel(zr)
%     sliceIndx = imgSetList(:,3) == zr(zrindx);% & ...
%     %rand(size(imgSetList,1),1) < 0.1;
%     w = scatter3(imgSetList(sliceIndx,1),...
%         imgSetList(sliceIndx,2),...
%         imgSetList(sliceIndx,4),...
%         2,...
%         '.k','MarkerEdgeAlpha',.1);
%     grid on
%     axis([-2,2,-2,2,-2,2])
%     set(gca,'dataaspectratio',[1,1,1])
%     set(gca,'projection','perspective')
%     pause(0.125)
% end

% setup video
outputVideo = VideoWriter(fullfile(pwd,'zi_slices.avi'));
outputVideo.FrameRate = 8;
open(outputVideo)

figure(2)
for kzi = 1:1:numel(zi)
    fv = isosurface(cr,ci,zr,img(:,:,:,kzi),iterDepth - .5);
    plot3(0,0,0,'.k')
    p1 = patch(fv,'FaceColor','black','EdgeColor','none');
    view(-37.5+180+(360*kzi/numel(zi)),30)
    daspect([1,1,1])
    axis([-2,2,-2,2,-2,2])
    %     axis tight
    camlight
    camlight(-80,-10)
    lighting gouraud
    grid on
    title(sprintf('%d',kzi))
    hold off;
    set(gca,'projection','perspective')
    pause(0.125)
    saveas(2,'temp.jpg')
    
    thisFrame = imread(fullfile(pwd,'temp.jpg'));
    writeVideo(outputVideo,thisFrame)
end
close(outputVideo)


% setup video
outputVideo = VideoWriter(fullfile(pwd,'zr_slices.avi'));
outputVideo.FrameRate = 8;
open(outputVideo)

figure(2)
for kzr = 1:1:numel(zr)
    fv = isosurface(cr,ci,zi,squeeze(img(:,:,kzr,:)),iterDepth - .5);
    plot3(0,0,0,'.k')
    p1 = patch(fv,'FaceColor','black','EdgeColor','none');
    view(3)
    daspect([1,1,1])
    axis([-2,2,-2,2,-2,2])
    camlight
    camlight(-80,-10)
    lighting gouraud
    grid on
    title(sprintf('%d',kzr))
    hold off;
    set(gca,'projection','perspective')
    pause(0.125)
    saveas(2,'temp.jpg')
    
    thisFrame = imread(fullfile(pwd,'temp.jpg'));
    writeVideo(outputVideo,thisFrame)
end
close(outputVideo)