function f = shd_error_opt(x,wl,aw,bbw,solz,a0,a1,alw_shade)

    [~,alphaw_shade,~,~,~]=model_alphaw(x,wl,aw,bbw,solz,a0,a1);    
 
    f=get_costfunction_error(alphaw_shade,alw_shade);
  
end

