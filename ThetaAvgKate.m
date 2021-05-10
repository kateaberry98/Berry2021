function Dir = ThetaAvgKate(THETA)
% Calculates average theta from nearest neighbour angles and assigns to a
% direction

xangleavg = mean(cos(THETA));
yangleavg = mean(sin(THETA));

if abs(yangleavg)>1e-14 || abs(xangleavg)>1e-14
    
    AvgT = atan(yangleavg/xangleavg)+pi.*(xangleavg<0);
    
    %pick direction
    if AvgT >-pi/4 && AvgT<=pi/4
        Dir = 'E';
    elseif AvgT >pi/4 && AvgT<=3*pi/4
        Dir = 'N';
    elseif AvgT >3*pi/4 && AvgT<=5*pi/4
        Dir = 'W';
    elseif AvgT >5*pi/4 || AvgT<=-pi/4
        Dir = 'S';
    else
        Dir = 'O'; %balanced vector means no angle / bias
    end
else
    Dir = 'O';
end
