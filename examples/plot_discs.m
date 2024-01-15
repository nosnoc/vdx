function plot_discs(h,pos,r)
    figure
    axis equal
    xlim([-50,50])
    ylim([-50,50])
    n = length(r);
    circles = {};
    tau = linspace(0, 2*pi)';
    vert = @(p,r) [r*cos(tau)+p(1),r*sin(tau)+p(2)];
    for ii=1:n
        p = pos((2*ii-1):(2*ii),1);
        v = vert(p,r(ii));
        circles{ii} = patch(v(:,1),v(:,2),'g');
    end
    pause(h(1));
    for jj=2:(length(h))
        for ii=1:n
            p = pos((2*ii-1):(2*ii),jj+1);
            circles{ii}.Vertices = vert(p,r(ii));
        end
        drawnow limitrate;
        %pause(h(jj));
        pause(0.1)
    end
end
