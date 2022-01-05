function [data, res] = add_slot_IP(prob, data, command, varargin)
%DUFF_ADD_IP   Slot function: Add to bifurcation data.
%
% Store initial point on first trajectory segment to bifurcation data

res = {};
switch command
  case 'init'
    res   = 'X0';
  case 'data'
    chart = varargin{1};
    [data, uidx] = coco_get_func_data(prob, 'po.orb.coll', ...
      'data', 'uidx');
    res  = chart.x(uidx(data.coll_seg.maps.x0_idx));
end

end
