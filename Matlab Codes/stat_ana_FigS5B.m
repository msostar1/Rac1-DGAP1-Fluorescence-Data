% A comparison of the measured circumferences of cells showcasing first-order oscillatory patterns,
% second-order oscillatory patterns, and first-order stationary patterns

varName = 'data_stats_mod_v_size';
path1 = 'Insert folder path'; 
load(strcat(path1, varName));  

siz1 = length(perim_one);
siz2 = length(perim_two);
siz3 = length(perim_polar);

g1 = ones(siz1, 1);
g2 = 2 * ones(siz2, 1);
g3 = 3 * ones(siz3, 1);
gg = vertcat(g1, g2, g3);

perim = vertcat(perim_one', perim_two', perim_polar');

[p_KW, table, stats] = kruskalwallis(perim, gg);
c = multcompare(stats, 'CType', 'tukey-kramer', 'display', 'off'); % last column displays p-values (first two columns display compared groups)