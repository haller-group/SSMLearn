function [data, res] = add_slot_hspo_FRC(prob, data, command, varargin)
% Slot function: Add to bifurcation data.

res = {};
switch command
    case 'init'
        res = { 'seg1.MAX|X|', 'seg2.MAX|X|', 'seg1.period' , ...
                'seg2.period', 'seg1.X0', 'seg2.X0', 'phase', 'switches' };
    case 'data'
        
        % Extract segment 1
        chart = varargin{1};
        [data, uidx] = coco_get_func_data(prob, 'hspo.orb.bvp.seg1.coll', ...
            'data', 'uidx');
        T1 = chart.x(uidx(data.coll_seg.maps.T_idx));
        t1 = transpose(data.coll_seg.mesh.tbp)*T1;
        x1 = chart.x(uidx(data.coll_seg.maps.xbp_idx));
        n = length(x1)/length(t1);
        x1 = reshape(x1,n,length(t1));
        % Extract segment 2
        [data, uidx] = coco_get_func_data(prob, 'hspo.orb.bvp.seg2.coll', ...
            'data', 'uidx');
        T2 = chart.x(uidx(data.coll_seg.maps.T_idx));
        t2 = transpose(data.coll_seg.mesh.tbp)*T2;
        x2 = chart.x(uidx(data.coll_seg.maps.xbp_idx));
        x2 = reshape(x2,n,length(t2));
        % Compute values
        amps_sol = [max(abs(x1),[],2) max(abs(x2),[],2)]; 
        x0_sol = [x1(:,2) x2(:,2)]; pers_sol = [T1 T2];
        
        
        x = [x1(1:end-2,:) x2(1:end-2,:)]; %x = [x1(1,:) x2(1,:)];
        t = [t1 t2+t1(end)];
        cosps = [x1(end-1,:) x2(end-1,:)]; sinps = [x1(end,:) x2(end,:)];
        fun_vals = x.*repmat((cosps-1i*sinps),n-2,1);
        z_temp = 0.5*sum( ( fun_vals(:,1:end-1)+fun_vals(:,2:end) ).*...
                                 repmat(diff(t),n-2,1), 2 ) / (T1+T2);
        phs_sol = atan2(imag(z_temp), real(z_temp));
        % Assume switching function equal to x(ndof+1) = 0;
        ndofs = (n-2)/2;
        q1_dot_seg_1_in = x1(ndofs+1,2:end-1);
        q1_dot_seg_2_in = x2(ndofs+1,2:end-1);
        if sum(sign(q1_dot_seg_1_in))>0
            q1_dot_seg_1_in_orig = q1_dot_seg_1_in;
            q1_dot_seg_1_in = q1_dot_seg_2_in;
            q1_dot_seg_2_in = q1_dot_seg_1_in_orig;
        end
        n_switches_sol = 2+(sum(q1_dot_seg_1_in>0)>0)+...
                           (sum(q1_dot_seg_2_in<0)>0);
        res = { amps_sol(:,1), amps_sol(:,2), pers_sol(1), pers_sol(2), ...
                x0_sol(:,1), x0_sol(:,2), phs_sol, n_switches_sol };
end

end
