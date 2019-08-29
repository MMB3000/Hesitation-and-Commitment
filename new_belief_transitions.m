
%%
n_recty = 14;                 %maximum time allowed (must be whole number)
alpha=0.25;

chunk_duration_msec=250;      %in miliseconds (includes rise and fall time)

t_num = n_recty/(chunk_duration_msec/1000);

[ gn, gpn, rn, rpn, pggn ] = new_belief( alpha, t_num );
%% Plot with regular scaling
figure
timestep = 5;
imagesc(gn{timestep},gpn{timestep},pggn{timestep}')
colorbar()
xticks(round(gn{timestep},2));
yticks(round(gpn{timestep},2));
xlabel('g');
ylabel('g''');
title('p(g''|g) n=5');
%% Plot adjusting the size of the squares to g and g'
figure
timestep = 5; %very computationally heavy for more than 20
X = gn{timestep};
Y = gpn{timestep};
spacing = round(diff(Y)/min(diff(Y)));
Y_spaced = Y(1):min(diff(Y)):Y(end);
data_spaced = repelem(pggn{timestep}',spacing([1 1:end]),1,1);
Y = Y_spaced;
spacing = round(diff(X)/min(diff(X)))';
X_spaced = X(1):min(diff(X)):X(end);
data_spaced_2 = repelem(data_spaced,1,spacing([1 1:end]));

imagesc(X_spaced,Y_spaced,data_spaced_2)
colorbar()
xlabel('g');
ylabel('g''');
title('p(g''|g) n=5');
xticks(round(gn{timestep},2));
yticks(round(gpn{timestep},2));

