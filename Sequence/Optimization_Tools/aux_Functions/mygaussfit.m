%% Gaussian fit for a given distribution y with initial conditions x0

function [mu,sg,curve] = mygaussfit(t,y,x0,noplot)

 % non-linear least squares problem (lsqnonlin)
    opts = optimset('display','off');
    [xsol]=lsqnonlin(@(x) Ob_gauss(x,y,t),x0,[],[],opts);
    mu=xsol(1);
    sg=xsol(2);
    curve=Ob_gauss(xsol,0*t,t);
    if ~strcmpi(noplot,'noplot')
        figure, plot(t,y);
        hold on
        plot(t,curve,'--r')
        legend('measured','fitting'), title(['Mu: ' num2str(mu) ' Std: ' num2str(sg)])
    end
end
