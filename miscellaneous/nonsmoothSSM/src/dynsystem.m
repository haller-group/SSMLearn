function dy = dynsystem(y, p, mode,Fl,Fr,fvec)
%PIECEWISE   'hspo'-compatible encoding of vector field.

% Parameters
ffre  = p(1,:);
famp  = p(2,:);

% Variables
x  = y(1:end-2,:);
f1 = y(end-1,:);
f2 = y(end,:);
f  = f1.^2 + f2.^2;

dy = zeros(size(y,1),size(y,2));
switch mode
    case 'left'
        dy(1:end-2,:) = Fl(x) + fvec*famp.*f1;
        dy(end-1,:) = f1 - ffre.*f2 - f.*f1;
        dy(end,:)   = f2 + ffre.*f1 - f.*f2;
    case 'right'
        dy(1:end-2,:) = Fr(x) + fvec*famp.*f1;
        dy(end-1,:) = f1 - ffre.*f2 - f.*f1;
        dy(end,:)   = f2 + ffre.*f1 - f.*f2;
end

end
