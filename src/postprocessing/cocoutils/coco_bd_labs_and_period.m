function [labs,perls,normls] = coco_bd_labs_and_period(bd, pt_type)

if nargin<2
  pt_type = [];
end

if ischar(pt_type)
	if strcmpi(pt_type, 'all')
		pt_type = '';
  end
end

labs   = coco_bd_col(bd, 'LAB', false);
perls  = coco_bd_col(bd, 'po.period');
normls = coco_bd_col(bd, '||po.orb.x||_{L_2[0,T]}');

if isempty(pt_type)
  labs = [ labs{:} ];
else
  types = coco_bd_col (bd, 'TYPE');
  idx   = strcmp(pt_type, types);
  labs  = [ labs{idx} ];
  perls  = perls(find(idx==1));
  normls  = normls(find(idx==1));
end
