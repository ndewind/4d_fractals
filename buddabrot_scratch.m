iterDepth = 2^8;

goodZ = complex(nan(iterDepth,1000),nan(iterDepth,1000));
for k = 1:1000
    thiscr = rand*4-2;
    thisci = rand*4-2;
    
    c = complex(thiscr,thisci);
        zold = complex(0,0);
        niter = 0;
        iterStop = false;
        while niter < iterDepth & iterStop == false
            niter = 1+niter;
            znew = zold^2 + c;
            if abs(znew) > 2
                iterStop = true;
            else
                goodZ(niter,k) = znew;
            end
            zold = znew;
        end
end
badZ = isnan(goodZ(end,:));
goodZ = goodZ(~badZ);
                plot(goodZ,'.k');
set(gca,'dataaspectratio',[1,1,1])