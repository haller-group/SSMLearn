function [data, res] = add_slot_Traj(prob, data, command, varargin)
% Store full trajectory segment

res = {};
switch command
  case 'init'
    res   = 'Traj';
  case 'data'
    chart = varargin{1};
    [data, uidx] = coco_get_func_data(prob, 'po.orb.coll', ...
      'data', 'uidx');
    res  = struct('StateVector',chart.x(uidx(data.coll_seg.maps.xbp_idx)), ... 
                   'NormalizedTime',data.coll_seg.mesh.tbp);
end

end
