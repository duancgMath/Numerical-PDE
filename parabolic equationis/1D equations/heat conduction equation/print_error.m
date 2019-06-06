function print_error(model,dt_vec,dx_vec,err_vec,n_grid_vec)
% PRINT_ERROR

n = length(dx_vec);

% ----------------------------------
% print order
fprintf('\t Table of error (T=%4.2f)\n', model.t_total);
fprintf('-------------------------------------------------\n');
fprintf('  #dt \t\t #dx \t\t #Error \t #Order \n');
fprintf('-------------------------------------------------\n');
for i = 1:n
    fprintf('%6.3e \t', dt_vec(i));
    fprintf('%6.3e \t', dx_vec(i));
    fprintf('%6.3e \t', err_vec(i));
    if (i == 1)
        fprintf(' - \n');
    else
        fprintf('%6.3f \n', ...
            log(err_vec(i)/err_vec(i-1))/log(dx_vec(i)/dx_vec(i-1)));
    end
end
fprintf('-------------------------------------------------\n');


% ----------------------------------
% plot 
figure
loglog(n_grid_vec,err_vec,'.','markersize',15)
hold on
loglog(n_grid_vec,n_grid_vec.^(-2),'--') 
grid on

end

