function [t_perturbed,R_perturbed] = perturb(t,R,amplitude)
    
    t_perturbed = t;
	r = rodrigues(R);
    
    for i=1:3
        t_perturbed(i,1) = t_perturbed(i,1) + (rand-0.5)*2.0*amplitude;
        r(i,1) = r(i,1) + (rand-0.5)*2.0*amplitude;
    end
    
    R_perturbed = rodrigues(r);
end