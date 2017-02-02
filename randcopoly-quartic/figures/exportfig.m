function exportfig(filename)
addpath('figures/export_fig')

set(gcf, 'Color', 'white'); % white bckgr
export_fig( gcf, ...      % figure handle
    filename,... % name of output file without extension
    '-painters', ...      % renderer
    '-jpg', ...           % file format
    '-r72' );             % resolution in dpi