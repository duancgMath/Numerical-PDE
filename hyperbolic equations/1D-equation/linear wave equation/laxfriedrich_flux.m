function num_flux = laxfriedrich_flux(model,u,v)
% LAXFRIEDRICH_FLUX
num_flux = 0.5*model.a*(u+v)-0.5*(model.dx/model.dt)*(v-u);
end

