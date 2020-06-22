% mandelbrot set scratch
% draws some videos from a pre calculated 4 cube
load('4_cube_64 center_n0.2566_n0.7679 size_00005.mat','img','cr','ci','zr','zi','iterDepth')
res = size(img);

% zi translation video
outputVideo = VideoWriter(fullfile(pwd,'zi_slices_zoom1.avi'));
outputVideo.FrameRate = 16;
open(outputVideo)

f = figure(2);
delete(f.Children)
figdim = [0,3,19,15];
set(f,'units','centimeter','position',figdim,'paperunits','centimeter','paperposition',figdim)
axis3D = axes('Position', [.0 .25 0.65 0.65]);
plot3(axis3D,0,0,0,'.w'); % just using this to set the axis to a 3D plot
axis_crxci = axes('Position',[.65,.72,.25,.2]);
axis_crxzr = axes('Position',[.65,.42,.25,.2]);
axis_cixzr = axes('Position',[.65,.1,.25,.2]);
for kzi = 1:1:numel(zi)
    set(f,'CurrentAxes',axis3D)
    fv = isosurface(cr,ci,zr,img(:,:,:,kzi),iterDepth - .5);
    cla % clear the old patch
    p1 = patch(axis3D,fv,'FaceColor','black','EdgeColor','none');
    p_crxci = patch(axis3D,[-2,2,2,-2],[-2,-2,2,2],repmat(ceil(res(3)/2),1,4),...
        'faceColor','red','EdgeColor','none','FaceAlpha',.3);
    p_crxzr = patch('XData',repmat(cr(ceil(res(2)/2)),1,4),...
        'YData',[-2,2,2,-2],...
        'ZData',[-2,-2,2,2],...
        'faceColor','green','EdgeColor','none','FaceAlpha',.3);
    p_cixzr = patch('XData',[-2,2,2,-2],...
        'YData',repmat(ci(ceil(res(1)/2)),1,4),...
        'ZData',[-2,-2,2,2],...
        'faceColor','blue','EdgeColor','none','FaceAlpha',.3);
    view(axis3D,-37.5+180+(360*kzi/numel(zi)),30)
%     daspect(axis3D,[1,1,1])
    axis(axis3D,[min(cr),max(cr),min(ci),max(ci),min(zr),max(zr)])
    camlight(axis3D)
    camlight(axis3D,-80,-10)
    lighting gouraud
    grid on
    title(sprintf('zi = %0.3f',zi(kzi)))
    xlabel('ci');ylabel('cr');zlabel('zr');
    set(gca,'projection','perspective')
    
    set(f,'CurrentAxes',axis_crxci)
    imagesc(cr,ci,img(:,:,ceil(res(3)/2),kzi)'==iterDepth)
    set(axis_crxci,'dataaspectratio',[1,1,1],'YDir','normal','box','off')
    colormap(axis_crxci,[1,1,1;1,0,0])
    xlabel('cr');ylabel('ci')
    
    set(f,'CurrentAxes',axis_crxzr)
    imagesc(cr,zr,squeeze(img(:,ceil(res(3)/2),:,kzi))'==iterDepth)
    set(axis_crxzr,'YDir','normal','box','off')
    colormap(axis_crxzr,[1,1,1;0,1,0])
    xlabel('cr');ylabel('zr')
    
    set(f,'CurrentAxes',axis_cixzr)
    imagesc(ci,zr,squeeze(img(ceil(res(3)/2),:,:,kzi))'==iterDepth)
    set(axis_cixzr,'YDir','normal','box','off')
    colormap(axis_cixzr,[1,1,1;0,0,1])
    xlabel('ci');ylabel('zr')
    
    pause(0.125)
    saveas(2,'temp.jpg')
    
    thisFrame = imread(fullfile(pwd,'temp.jpg'));
    writeVideo(outputVideo,thisFrame)
end
close(outputVideo)