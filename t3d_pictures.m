function out = t3d_pictures(varargin)
% T3D_PICTURES  The COMSOL pictures (and the final .mph with all plot groups) from a SAVED model, without re-solving.
% Use it when t3d_run stopped or was interrupted after the solves (the early <case>_3d_solved.mph exists), or to redraw
% the pictures with other settings (plain lines instead of tubes, no per-rail plots, ...).
%
%   t3d_pictures('tag', 'run30')                 reads out_run30\<case>_3d_solved.mph (or <case>_3d.mph), writes out_run30\comsol_images\
%   t3d_pictures('tag', 'run30', 'tubes', 0)     rails as plain lines (fast; ev6 needs this)
%   t3d_pictures('tag', 'run30', 'rail_plots', 0)
%   t3d_pictures('tag', 'run30', 'model', 'T3D') the model t3d_run left in the COMSOL server (after Ctrl+C), nothing loaded from disk
%
% OPTIONS   'tag' (output folder out_<tag>), 'mph' (file; default the solved / final .mph of that folder), 'tubes' ('auto' | 0 | 1),
%           'rail_plots' (1), 'topview' (1), 'tend_years' (colour range of the nucleation plot; default: last horizon), 'nper' (60),
%           'save_mph' (1 = save <case>_3d.mph with the plot groups at the end)
% Everything else (rails, horizons, output times) is read from data\<case>\ and from the csv files that t3d_run wrote
% into the folder (t3d_rail_result.csv, t3d_rail_smax.csv, t3d_nodes.csv).

import com.comsol.model.*
import com.comsol.model.util.*

here = fileparts(mfilename('fullpath'));  Y = 365.25 * 86400;
dd_ = dir(fullfile(here, 'data'));  dd_ = dd_([dd_.isdir] & ~ismember({dd_.name}, {'.', '..'}));
if numel(dd_) ~= 1, error('t3d:data', 'expected exactly one case folder under %s, found %d', fullfile(here, 'data'), numel(dd_)); end
o = struct('tag','run', 'mph','', 'model','', 'tubes','auto', 'rail_plots',1, 'topview',1, 'tend_years',[], 'nper',60, 'save_mph',1, 'data', fullfile(here, 'data', dd_(1).name));
for k = 1:2:numel(varargin), o.(varargin{k}) = varargin{k+1}; end
outDir = fullfile(here, ['out_' o.tag]);  imgDir = fullfile(outDir, 'comsol_images');  if ~exist(imgDir, 'dir'), mkdir(imgDir); end
fb = {};  t0 = tic;

R  = dlmread(fullfile(o.data, 'rails.csv'), ',', 1, 0);
S  = dlmread(fullfile(o.data, 'segments.csv'), ',', 1, 0);
PB = dlmread(fullfile(o.data, 'power_blocks.csv'), ',', 1, 0);
C  = jsondecode(fileread(fullfile(o.data, 'constants.json')));
K  = jsondecode(fileread(fullfile(o.data, 'stack.json')));
caseName = K.case;  dieNames = K.die_names;  if ~iscell(dieNames), dieNames = cellstr(dieNames); end
sc = C.em.sigma_crit;  nd = numel(dieNames);
% ---------------- what t3d_run solved: rails, output times, horizons ----------------
RR = dlmread(fullfile(outDir, 't3d_rail_result.csv'), ',', 1, 0);          % rail, die, isV, index, n_seg, t_nuc, smax...
[~, idx] = ismember(RR(:, 1), R(:, 1));
if any(idx == 0), error('t3d:data', 't3d_rail_result.csv names rails that are not in rails.csv'); end
nr = numel(idx);  Sr = R(idx, 5);  W = R(idx, 8);  die = R(idx, 2);  isV = R(idx, 3) == 1;  x0r = R(idx, 6);  y0r = R(idx, 7);
hdr = fgetl(fopen(fullfile(outDir, 't3d_rail_smax.csv')));  tl = str2double(strsplit(hdr, ','));  tl = tl(2:end);  nt = numel(tl);
hn = fgetl(fopen(fullfile(outDir, 't3d_nodes.csv')));  tok = regexp(hn, 'sigma_([0-9.]+)yr_Pa', 'tokens');
hz = cellfun(@(c) str2double(c{1}), tok);  nh = numel(hz);
ih = zeros(1, nh);  for h = 1:nh, [~, ih(h)] = min(abs(tl - hz(h) * Y)); end
if isempty(o.tend_years), o.tend_years = max(hz); end
fclose('all');
% ---------------- layers and rail levels (as in t3d_run, real thicknesses) ----------------
LJ = K.layers;  if isstruct(LJ), LJ = num2cell(LJ); end
L = struct('name', {}, 'kind', {}, 't', {}, 'sx', {}, 'sy', {}, 'die', {});
for i = 1:numel(LJ)
    j = LJ{i};  di = -1;  if strcmpi(j.kind, 'die'), di = j.die_index; end
    L(end+1) = struct('name', j.name, 'kind', lower(j.kind), 't', j.thickness, 'sx', j.size_x, 'sy', j.size_y, 'die', di); %#ok<AGROW>
end
nl = numel(L);  zb = [0, cumsum([L.t])];  dieL = zeros(1, nd);
for i = 1:nl, if L(i).die >= 0, dieL(L(i).die + 1) = i; end, end
ztop = zb(dieL + 1);  dz = 2e-6;  zH = ztop + dz;  zV = ztop + 2 * dz;
zr = zH(die + 1)';  zr(isV) = zV(die(isV) + 1)';
SX = K.die_size_x;  SY = K.die_size_y;  cx = 0.5 * SX;  cy = 0.5 * SY;
stackL = find(~ismember({L.kind}, {'spreader', 'sink'}));
aj = cell(nr, 1);
for q = 1:nr
    Ls = S(R(idx(q), 14) + (1:Sr(q)), 3);
    if isV(q), aj{q} = y0r(q) + [0; cumsum(Ls)]; else, aj{q} = x0r(q) + [0; cumsum(Ls)]; end
end
across = x0r .* isV + y0r .* ~isV;
% ---------------- the saved model ----------------
mphf = o.mph;
if ~isempty(o.model)                                           % a model that is still in the COMSOL server (t3d_run created 'T3D')
    model = ModelUtil.model(o.model);  mphf = ['server model ' o.model];  fprintf('using %s\n', mphf);
else
    if isempty(mphf)
        mphf = fullfile(outDir, sprintf('%s_3d_solved.mph', caseName));
        if ~exist(mphf, 'file'), mphf = fullfile(outDir, sprintf('%s_3d.mph', caseName)); end
    end
    if ~exist(mphf, 'file'), error('t3d:mph', 'no saved model: %s  (if the model is still in the server: t3d_pictures(..., ''model'', ''T3D''))', mphf); end
    fprintf('loading %s ...\n', mphf);  model = mphload(mphf);  fprintf('loaded (%.0f s)\n', toc(t0));
end
% a model saved with plot groups (the final .mph): start again from a clean result node
try, for t = cell(model.result.tags())', model.result.remove(t{1}); end, catch, end
try, for t = cell(model.result.export.tags())', model.result.export.remove(t{1}); end, catch, end
try, for t = cell(model.result.dataset.tags())', if ~isempty(regexp(t{1}, '^dset(Stack|Rails)', 'once')), model.result.dataset.remove(t{1}); end, end, catch, end
try, model.view.remove('vtop'); catch, end
hasT = false;  try, model.physics('ct');  hasT = true; catch, end
uname = 'u';  try, uname = char(model.physics('ce').field('dimensionless').field()); catch, end
if hasT, Tn = 'Tt';  try, Tn = char(model.physics('ct').field('dimensionless').field()); catch, end, else, Tn = 'Thsv'; end
tg = cell(model.result.dataset.tags());  dsU = '';  dsT = '';
for k = 1:numel(tg)
    try, info = mphsolinfo(model, 'dataset', tg{k});  if numel(info.solvals) == nt, dsU = tg{k}; elseif isempty(dsT), dsT = tg{k}; end, catch, end
end
if isempty(dsU), error('t3d:mph', 'no data set with %d output times in %s', nt, mphf); end
if ~hasT || isempty(dsT), dsT = dsU; end
% ---------------- selections by coordinates (as in t3d_run) ----------------
m = 1e-4;  railE = [];  railEq = cell(nr, 1);  dieE = cell(nd, 1);
for d = 1:nd
    i = dieL(d);
    dieE{d} = mphselectbox(model, 'geom1', [cx - 0.5 * L(i).sx - m, cx + 0.5 * L(i).sx + m; cy - 0.5 * L(i).sy - m, cy + 0.5 * L(i).sy + m; ztop(d) + 0.6 * dz, ztop(d) + 2.4 * dz], 'edge');
    railE = [railE, dieE{d}]; %#ok<AGROW>
end
railE = unique(railE);
if numel(railE) ~= sum(Sr), error('t3d:sel', 'found %d rail edges, expected %d', numel(railE), sum(Sr)); end
ua = sort(unique(across));  if numel(ua) > 1, wsel = 0.2 * min(diff(ua)); else, wsel = 20e-6; end
wsel = max(min(wsel, 20e-6), 1e-6);
for q = 1:nr
    a = aj{q};
    if isV(q), box = [x0r(q) - wsel, x0r(q) + wsel; a(1) - m, a(end) + m; zr(q) - 0.4 * dz, zr(q) + 0.4 * dz];
    else,      box = [a(1) - m, a(end) + m; y0r(q) - wsel, y0r(q) + wsel; zr(q) - 0.4 * dz, zr(q) + 0.4 * dz]; end
    railEq{q} = mphselectbox(model, 'geom1', box, 'edge');
    if numel(railEq{q}) ~= Sr(q), error('t3d:sel', 'rail %d: found %d edges, expected %d', R(idx(q), 1), numel(railEq{q}), Sr(q)); end
end
stackDoms = [];
for i = stackL
    e = 0.05 * min([L.t]);  bx = [cx - 0.5 * L(i).sx - m, cx + 0.5 * L(i).sx + m; cy - 0.5 * L(i).sy - m, cy + 0.5 * L(i).sy + m];
    v_ = mphselectbox(model, 'geom1', [bx; zb(i) - e, zb(i + 1) + e], 'domain');  stackDoms = [stackDoms, v_(:)']; %#ok<AGROW>
end
fprintf('%s: %d rails, %d output times, horizons %s; drawing ...\n', caseName, nr, nt, mat2str(hz));
% ---------------- plot groups and COMSOL-rendered images ----------------
pgs = {};
solT = '';  solU = '';
try, solU = char(model.result.dataset(dsU).getString('solution'));  solT = solU; catch, end
% data sets restricted to the die stack (temperature) and to the rails (stress): what the plots show, and what the
% automatic view frames.  The spreader and sink stay in the model; they only leave the pictures.
dsTd = dsT;  dsUr = dsU;
try
    d_ = model.result.dataset.create('dsetStack', 'Solution');  d_.set('solution', solT);
    d_.selection.geom('geom1', 3);  d_.selection.set(stackDoms);  dsTd = 'dsetStack';
    d_ = model.result.dataset.create('dsetRails', 'Solution');  d_.set('solution', solU);
    d_.selection.geom('geom1', 1);  d_.selection.set(railE);  dsUr = 'dsetRails';
    for d = 1:nd
        d_ = model.result.dataset.create(sprintf('dsetRails%d', d), 'Solution');  d_.set('solution', solU);
        d_.selection.geom('geom1', 1);  d_.selection.set(dieE{d});
    end
catch
    fb{end+1} = 'selection data sets not created: the COMSOL pictures show the whole package (zoom in the .mph)';
end
vtop = '';                                                    % an orthographic view from above for the per-die pictures
if o.topview
    try
        vw = model.view.create('vtop', 3);  Sd = max(SX, SY);  Dcam = 10 * Sd;  zt_ = zb(dieL(nd) + 1);
        vw.camera.set('projection', 'orthographic');
        vw.camera.set('position', [cx, cy, zt_ + Dcam]);  vw.camera.set('target', [cx, cy, zt_]);  vw.camera.set('up', [0 1 0]);
        vw.camera.set('rotationpoint', [cx, cy, zt_]);  vw.camera.set('viewoffset', [0 0]);
        vw.camera.set('zoomangle', 2 * atand(0.62 * Sd / Dcam));          % frames the die with a small margin
        vtop = 'vtop';
    catch
        fb{end+1} = 'top view not created: per-die pictures use the default 3-D view';
    end
end
try   % temperature of the whole package
    pg = model.result.create('pgT_all', 'PlotGroup3D');  pg.label('HotSpot temperature on the package (K; package faces show 0)');  pg.set('data', dsT);
    sf = pg.create('surf1', 'Surface');  sf.set('expr', Tn);  sf.set('colortable', 'Inferno');  pgs{end+1} = 'pgT_all';
catch, fb{end+1} = 'plot: temperature of the package not created'; end
try   % temperature of the die stack only
    pg = model.result.create('pgT_stack', 'PlotGroup3D');  pg.label('HotSpot temperature painted on the die stack (K)');  pg.set('data', dsTd);
    sf = pg.create('surf1', 'Surface');  sf.set('expr', Tn);  sf.set('colortable', 'Inferno');  pgs{end+1} = 'pgT_stack';
catch, fb{end+1} = 'plot: temperature of the die stack not created'; end
tube3 = sprintf('%.3g', 1.5e-3 * max(SX, SY));  tube2 = sprintf('%.3g', 2e-3 * max(SX, SY));     % tube radius of the rail plots (6 / 8 um on a 4-mm die)
% tubes turn every line element into a small cylinder: fine up to a few million elements, hopeless beyond (ev6: 5.6 M
% elements, one picture did not finish in two hours).  'tubes','auto' keeps them below the threshold, plain lines above.
if ischar(o.tubes), useTube = sum(Sr) * o.nper <= 2.5e6; else, useTube = logical(o.tubes); end
if ~useTube, fprintf('rail pictures with plain lines (%d line elements; tubes only up to 2.5e6)\n', sum(Sr) * o.nper); end
for h = 1:nh   % 3-D stress on the rails above grey dies, one picture per horizon
    try
        tgp = sprintf('pgS3d_%d', h);  pg = model.result.create(tgp, 'PlotGroup3D');  pg.label(sprintf('Rail stress after %g years, 3-D (Pa)', hz(h)));
        pg.set('data', dsUr);  pg.setIndex('looplevel', ih(h), 0);
        ln = pg.create('line1', 'Line');  ln.set('expr', uname);  ln.set('colortable', 'RainbowLight');
        if useTube, try, ln.set('linetype', 'tube');  ln.set('radiusexpr', tube3);  ln.set('tuberadiusscaleactive', true);  ln.set('tuberadiusscale', 1); catch, end
        else, try, ln.set('linetype', 'line');  ln.set('linewidth', 3); catch, end, end
        try, ln.set('rangecoloractive', true);  ln.set('rangecolormin', -sc);  ln.set('rangecolormax', sc); catch, end
        try
            sd = pg.create('surf1', 'Surface');  sd.set('data', dsTd);  sd.set('expr', '1');  sd.set('coloring', 'uniform');  sd.set('color', 'gray');
            sd.set('alpha', 0.25);  try, sd.set('solutionparams', 'parent'); catch, end
        catch
        end
        pgs{end+1} = tgp;
    catch, fb{end+1} = sprintf('plot: 3-D rail stress after %g years not created', hz(h)); end
end
for d = 1:nd   % per-die top view at every horizon (the 2-D stress map of that die's power grid)
    for h = 1:nh
        try
            tgp = sprintf('pgS2d_%s_%d', dieNames{d}, h);  pg = model.result.create(tgp, 'PlotGroup3D');
            pg.label(sprintf('%s rail stress after %g years, top view (Pa)', dieNames{d}, hz(h)));
            pg.set('data', sprintf('dsetRails%d', d));  pg.setIndex('looplevel', ih(h), 0);
            if ~isempty(vtop), pg.set('view', vtop); end
            ln = pg.create('line1', 'Line');  ln.set('expr', uname);  ln.set('colortable', 'RainbowLight');
            if useTube, try, ln.set('linetype', 'tube');  ln.set('radiusexpr', tube2);  ln.set('tuberadiusscaleactive', true);  ln.set('tuberadiusscale', 1); catch, end
        else, try, ln.set('linetype', 'line');  ln.set('linewidth', 3); catch, end, end
            try, ln.set('rangecoloractive', true);  ln.set('rangecolormin', -sc);  ln.set('rangecolormax', sc); catch, end
            pgs{end+1} = tgp;
        catch, fb{end+1} = sprintf('plot: top view of %s after %g years not created', dieNames{d}, hz(h)); end
    end
end
try   % nucleation time of every rail, 3-D
    pg = model.result.create('pgTn', 'PlotGroup3D');  pg.label('Nucleation time of every rail (years, log10)');  pg.set('data', dsUr);  pg.setIndex('looplevel', nt, 0);
    ln = pg.create('line1', 'Line');  ln.set('expr', 'log10(tnuc_yr)');  ln.set('colortable', 'Viridis');  ln.set('colortablerev', true);
    if useTube, try, ln.set('linetype', 'tube');  ln.set('radiusexpr', tube3);  ln.set('tuberadiusscaleactive', true);  ln.set('tuberadiusscale', 1); catch, end
        else, try, ln.set('linetype', 'line');  ln.set('linewidth', 3); catch, end, end
    try, ln.set('rangecoloractive', true);  ln.set('rangecolormin', log10(0.1));  ln.set('rangecolormax', log10(10 * o.tend_years)); catch, end
    pgs{end+1} = 'pgTn';
catch, fb{end+1} = 'plot: nucleation time not created'; end
nimg = 0;
for k = 1:numel(pgs)
    try
        im = model.result.export.create(sprintf('img%d', k), 'Image');
        im.set('sourceobject', pgs{k});  im.set('filename', fullfile(imgDir, [pgs{k} '.png']));
        try, im.set('size', 'manualweb');  im.set('unit', 'px');  im.set('width', '2400');  im.set('height', '1800'); catch, end
        im.run;  nimg = nimg + 1;
    catch
        fb{end+1} = ['COMSOL image of plot group ' pgs{k} ' not exported (open the .mph and export it from the plot window)'];
    end
end
if o.rail_plots   % one 1-D stress plot per rail, all horizons as curves
    railImg = fullfile(imgDir, 'rails');  if ~exist(railImg, 'dir'), mkdir(railImg); end
    nfail = 0;
    for q = 1:nr
        try
            tgp = sprintf('pgR%d', q);  pg = model.result.create(tgp, 'PlotGroup1D');
            pg.label(sprintf('rail %d (%s %s%d): stress along the rail', R(idx(q), 1), dieNames{die(q) + 1}, char('H' + isV(q) * ('V' - 'H')), R(idx(q), 4)));
            pg.set('data', dsU);  pg.set('looplevelinput', {'manual'});  pg.set('looplevel', {sprintf('%d ', ih)});
            lg = pg.create('lngr1', 'LineGraph');  lg.selection.set(railEq{q});  lg.set('expr', uname);
            lg.set('xdata', 'expr');  if isV(q), lg.set('xdataexpr', 'y'); else, lg.set('xdataexpr', 'x'); end
            try, lg.set('legend', true);  lg.set('legendmethod', 'automatic'); catch, end
            try, pg.set('xlabel', 'position along the rail (m)');  pg.set('ylabel', 'hydrostatic stress (Pa)');  pg.set('xlabelactive', true);  pg.set('ylabelactive', true); catch, end
            if q == 1, pgs{end+1} = tgp; end
            im = model.result.export.create(sprintf('imgR%d', q), 'Image');
            im.set('sourceobject', tgp);  im.set('filename', fullfile(railImg, sprintf('rail_%05d_%s_%s%d.png', R(idx(q), 1), dieNames{die(q) + 1}, char('H' + isV(q) * ('V' - 'H')), R(idx(q), 4))));
            try, im.set('size', 'manualweb');  im.set('unit', 'px');  im.set('width', '1600');  im.set('height', '900'); catch, end
            im.run;  nimg = nimg + 1;
            if q > 1, model.result.remove(tgp);  model.result.export.remove(sprintf('imgR%d', q)); end      % keep one in the .mph as a template
        catch
            nfail = nfail + 1;
        end
    end
    if nfail, fb{end+1} = sprintf('%d of %d per-rail plots could not be exported', nfail, nr); end
end
try, model.label(sprintf('%s stack: rail stress with HotSpot temperatures (no thermal solve)', caseName)); catch, end

sec_plots = toc(t0);
if o.save_mph
    try, mphsave(model, fullfile(outDir, sprintf('%s_3d.mph', caseName)));
    catch err, fb{end+1} = ['final .mph not saved (' regexprep(err.message, '\s+', ' ') ')']; end
end
meta = struct('case_name', caseName, 'source_mph', mphf, 'n_images', nimg, 'seconds', sec_plots, 'fallbacks', {fb});
fid = fopen(fullfile(outDir, 't3d_pictures_meta.json'), 'w');  fprintf(fid, '%s', jsonencode(meta));  fclose(fid);
fprintf('%d pictures in %s (%.0f s)\n', nimg, imgDir, sec_plots);
if ~isempty(fb), fprintf('API fall-backs taken:\n');  fprintf('   %s\n', fb{:}); end
out = meta;
end
