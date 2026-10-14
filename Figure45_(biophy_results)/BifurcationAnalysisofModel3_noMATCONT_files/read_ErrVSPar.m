clear
close all
matFiles = dir('CT=*_ErrVSPar_*.mat');
if isempty(matFiles)
    error('No files matching "CT=*_ErrVSPar_*.mat" found in %s', folderPath);
end

%% --- Parse (CT number, hill type, param) out of each filename ---
ctNumber = strings(numel(matFiles), 1);
hillType = strings(numel(matFiles), 1);
paramName = strings(numel(matFiles), 1);

pat = 'CT=([^_]+)_ErrVSPar_([A-Za-z0-9]+hill)_\{([^}]+)\}';

for i = 1:numel(matFiles)
    tok = regexp(matFiles(i).name, pat, 'tokens', 'once');
    if isempty(tok)
        warning('Skipping "%s": does not match CT=<num>_ErrVSPar_<type>hill_{<param>}_ pattern.', ...
            matFiles(i).name);
        continue;
    end
    ctNumber(i)  = string(tok{1});
    hillType(i)  = string(tok{2});
    paramName(i) = string(tok{3});
end

validMask = ctNumber ~= "";
matFiles  = matFiles(validMask);
ctNumber  = ctNumber(validMask);
hillType  = hillType(validMask);
paramName = paramName(validMask);

% Group key = combination of all three fields
groupKey = ctNumber + "||" + hillType + "||" + paramName;
uniqueKeys = unique(groupKey);

fprintf('Found %d unique (CT, hill type, param) group(s).\n', numel(uniqueKeys));

%%
for k = 1:numel(uniqueKeys)
    thisKey = uniqueKeys(k);
    idx = find(groupKey == thisKey);

    parts = split(thisKey, "||");
    thisCT = parts(1);
    thisHill = parts(2);
    thisParam = parts(3);
    label = sprintf('CT=%s_%s_{%s}', thisCT, thisHill, thisParam);

    all_err_combined = [];
    rangealpha_combined = [];

    for j = 1:numel(idx)
        fname = fullfile(pwd, matFiles(idx(j)).name);
        data = load(fname, 'all_err', 'rangealpha');

        if ~isfield(data, 'all_err') || ~isfield(data, 'rangealpha')
            warning('Skipping "%s": missing all_err or rangealpha.', matFiles(idx(j)).name);
            continue;
        end

        err_vec = data.all_err(:);
        alpha_vec = data.rangealpha(:);

        if numel(err_vec) ~= numel(alpha_vec)
            warning('Skipping "%s": all_err (%d elements) and rangealpha (%d elements) sizes do not match.', ...
                matFiles(idx(j)).name, numel(err_vec), numel(alpha_vec));
            continue;
        end

        all_err_combined = [all_err_combined; err_vec]; %#ok<AGROW>
        rangealpha_combined = [rangealpha_combined; alpha_vec]; %#ok<AGROW>
    end

    if isempty(all_err_combined)
        warning('No usable data found for %s, skipping.', label);
        continue;
    end

    % Sort by rangealpha for a cleaner plot (sortIdx applied to both
    % arrays identically, so the pairing between all_err and rangealpha
    % is preserved)
    [rangealpha_sorted, sortIdx] = sort(rangealpha_combined);
    all_err_sorted = all_err_combined(sortIdx);

    % Save combined data
    outMat = fullfile(pwd, sprintf('%s_combined.mat', label));
    save(outMat, 'all_err_sorted', 'rangealpha_sorted');

    % Plot
    fig = figure('Name', label);
    plot(rangealpha_sorted, all_err_sorted, 'LineWidth',3,'Color','k');
    hold on;

    % Overlay mean +/- std per unique rangealpha value
    % [uAlpha, ~, grp] = unique(rangealpha_sorted);
    % meanErr = accumarray(grp, all_err_sorted, [], @mean);
    % stdErr  = accumarray(grp, all_err_sorted, [], @std);
    % errorbar(uAlpha, meanErr, stdErr, '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'auto');

    xlabel(char(thisParam));
    ylabel('all_err');
    title(sprintf('%s', label), ...
        'Interpreter', 'none');
    legend({'all points'}, 'Location', 'best');
    set(gca,'Yscale','log')
    % grid on;
    hold off;

    outFig = fullfile(pwd, sprintf('%s_all_err_vs_rangealpha.fig', label));
    saveas(fig, outFig);
    % fprintf('%s: combined %d points from %d file(s). Saved plot to %s\n', ...
    %     label, numel(all_err_sorted), numel(idx), outFig);
legend('off')
end