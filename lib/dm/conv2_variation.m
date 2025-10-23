

function dmg = conv2_variation(dmny, dmnx, xyzv, inf, stat)
  if ~exist('stat')
    stat.width.percent  = 0;        % vary inf width
    stat.position.pixel = 0;    % vary inf position 
  end

  if any(stat.width.percent) | any(stat.position.pixel)
    if stat.rseed == 0  % Static
        randvalue   = ones(3,size(xyzv,1));
    else
        rng("default")
        rng(stat.rseed)
        randvalue   = randn(3,size(xyzv,1));
    end
  end 


  dmg  = zeros(dmny, dmnx);             % initialize subsampled DM grid
  if ~rem(inf,2) % even 
    error('length(inf) should be odd');
  end
  hl    = (length(inf)-1)/2;


  for ii = 1:size(xyzv,1),
    if any(xyzv(ii,3))
        xinx    = xyzv(ii,1)+[-hl:hl]; 
        yinx    = xyzv(ii,2)+[-hl:hl]; 
        if any(stat.width.percent) | any(stat.position.pixel)
            inff            = inf_modify(inf, stat.width.percent*randvalue(1,ii), stat.position.pixel*randvalue(2:3,ii));
            dmg(yinx,xinx)  = dmg(yinx,xinx) + xyzv(ii,3)*inff;
            continue;
        end
        dmg(yinx,xinx)  = dmg(yinx,xinx) + xyzv(ii,3)*inf;
    end
  end
return

