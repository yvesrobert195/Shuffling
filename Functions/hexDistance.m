function d=hexDistance(startx,starty, destx,desty)
  if (startx == destx)
    d=abs(desty - starty);
  elseif (starty == desty)
    d=abs(destx - startx);
  else
    dx = abs(destx - startx);
    dy = abs(desty - starty);
    if starty < desty
      d=dx + dy - (ceil(dx/2.0));
    else
      d=dx + dy - (floor(dx/2.0));
    end
  end