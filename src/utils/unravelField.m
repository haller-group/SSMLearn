function [Vx, Vy, Vz] = unravelField(field, shape)
    % Return the U, V, W components of a flow field which have been stacked
    % on top of each other.
    NN = floor(size(field,1)/3);
    Vx = field(1:NN);
    Vy = field(NN+1:2*NN);
    Vz = field(2*NN+1:end);
    Vx = reshape(Vx, shape);
    Vy = reshape(Vy, shape);
    Vz = reshape(Vz, shape);
end
