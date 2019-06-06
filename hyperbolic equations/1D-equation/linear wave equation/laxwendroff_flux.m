function num_flux = laxwendroff_flux(model,u,v)
% LAXWENDROFF_FLUX
temp = 0.5*(u+v)-0.5*model.a*(model.dt/model.dx)*(v-u);
num_flux = model.a*temp;
end

