function plot_range_growthrate_patch(x,real_data_max,real_data_min,mycolors)
fill([x,fliplr(x)], [real_data_max,fliplr(real_data_min)],mycolors,'EdgeColor','none','FaceAlpha','0.1');

