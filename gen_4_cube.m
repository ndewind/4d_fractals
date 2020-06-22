% generate high res 4-cube


% res = [512,512,512,512];
res = [256,256,256,256];
% res = [64,64,64,64];

center = [-0.2566,-0.7679];
% center = [-0,-0];
sz = [.00005,.00005];
% sz = [2,2];
cr = linspace(center(1)-sz(1),center(1)+sz(1),res(1));
ci = linspace(center(2)-sz(2),center(2)+sz(2),res(2));
zSizeFactor = 200;
zr = linspace(-sz(1)*zSizeFactor,sz(1)*zSizeFactor,res(3));
zi = linspace(-sz(1)*zSizeFactor,sz(1)*zSizeFactor,res(4));
iterDepth = 2^8;
img = zeros(res,'uint16');
tic
for crindx = 1:numel(cr)
    fprintf('%05d/%05d\n',crindx,numel(cr))
    toc
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
                    img(crindx,ciindx,zrindx,ziindx) = niter;
                    if abs(znew) > 2
                        iterStop = true;
                    end
                    zold = znew;
                end
                
            end
        end
    end
end
%  
% for k = 1:numel(zi)
%     figure(1)
%     
%     subplot(1,2,1)
%     imagesc(ci,cr,img(:,:,ceil(res(3)/2),k),[0,iterDepth])
%     set(gca,'dataaspectratio',[1,1,1])
%     colormap(gray)
%     
%     subplot(1,2,2)
%     imagesc(ci,cr,img(:,:,k,ceil(res(3)/2)),[0,iterDepth])
%     set(gca,'dataaspectratio',[1,1,1])
%     colormap(gray)
%     
%     pause(1)
% end

save('4_cube_64 center_n0.2566_n0.7679 size_00005.mat','img','cr','ci','zr','zi','iterDepth','-v7.3')