function varargout = fdra2c(varargin)
% WriteKvf(fn, c, allow_overwrite)
%   c is a struct containing string and numeric fields. Write it to the file
% fn. This fill will be used by the program fdra2c. The file is not
% necessarily portable to other machines.
%
% c = ReadKvf(fn)
%   Load a struct from the file fn.
  [varargout{1:nargout}] = feval(varargin{:});
end

function WriteKvf (varargin)
  kvf('Write', varargin{:});
end

function c = ReadKvf (varargin)
  c = kvf('Read', varargin{:});
end

function t = qload_t (cf, fld, stride, do_stride)
  if (nargin < 2) fld = 'v'; end
  if (nargin < 3) stride = 1; end
  if (nargin < 4) do_stride = numel(stride) == 1; end
  fn = sprintf('%s_%s', cf.save_filename, fld);
  if (exist([fn '.tsdf'], 'file')) fn = [fn '.tsdf']; else fn = [fn '.sdf']; end
  t = SaveStreamData('Read', fn, stride, 1, do_stride);
end

function s = qload (cf, stride, varargin)
  o = popt(varargin,...
           {{'do_stride' numel(stride) == 1};
            {'want_slip' 1}; {'want_theta' 0}; {'want_dlte' 0}; {'want_p' 0};
            {'want_t' 0}; {'s'}; {'rid' -1}});
  o.want_v = 1;
  if (isempty(o.s)) s = GetGeometry(cf, o.rid);
  else s = o.s; end
  s.v_creep = cf.v_creep;

  function Load (fld)
    if (eval(sprintf('o.want_%s', lower(fld))))
      eval(sprintf('s.t_%s = qload_t(cf, ''%s'', stride, o.do_stride);',...
                   fld, fld));
      s.(fld) = SaveStreamData('Read', sprintf('%s_%s.sdf', cf.save_filename,...
                                               fld),...
                               stride, [], o.do_stride);
      s.(fld)(1,:) = [];
    else
      eval(sprintf('s.t_%s = [];', fld));
      s.(fld) = [];
    end
  end
  
  try
    Load('v');
    s.t = s.t_v;
    Load('slip'); Load('theta'); Load('dlte'); Load('p'); Load('T');
  catch
    lasterr
  end
  s.stride = stride;
  s.fn = cf.save_filename;
end

function s = CombineQloads (ss)
  assert(iscell(ss) && ~isempty(ss));
  if (len(ss) == 1) s = ss{1}; return; end
  for (i = 2:len(ss)) assert(size(ss{1}.v, 1) == size(ss{i}.v, 1)); end

  s = ss{1};
  s = rmfield(s, 'fn');
  s = rmfield(s, 'stride');
  s.fns = {}; s.strides = {};
  for (i = 1:len(ss))
    s.fns{i} = ss{i}.fn;
    s.strides{i} = ss{i}.stride;
  end
  flds = {'t' 'v' 't_slip' 'slip' 't_theta' 'theta'};
  for (i = 1:len(flds)) s.(flds{i}) = []; end
  
  function Combine (tf, vf)
    for (i = 1:len(ss))
      s.(tf) = [s.(tf) ss{i}.(tf)];
      s.(vf) = [s.(vf) ss{i}.(vf)];
    end
    [s.(tf) p] = unique(s.(tf));
    s.(vf) = s.(vf)(:, p);
  end
  Combine('t', 'v');
  Combine('t_slip', 'slip');
  Combine('t_theta', 'theta');
end

function s = GetGeometry (cf, rid)
  if (nargin < 2) rid = -1; end
  s.rmesh_filename = cf.rmesh_filename;
  if (rid < 0) rid = cf.rmesh_filename; end
  s.rmi = rmi_Init(rid);
  s.fh = @(I) fdra2c('rmi_Image', s.rmi, I);
  s.rmmi = rmi_Init(rid, 'use_max', 0);
  s.fhm = @(I) fdra2c('rmi_Image', s.rmmi, I);
  if (ischar(rid)) ridi = dc3dm.mRead(rid); else ridi = rid; end
  s.rs = dc3dm.mRects(ridi);
  if (ischar(rid)) dc3dm.mClear(ridi); end
  d = dc3dm.mData(s.rs);
  s.dx = d.dx; s.dy = d.dy; s.xlim = d.xlim; s.ylim = d.ylim;
  s.cs = (s.rs(1:2,:) + 0.5*s.rs(3:4,:))';
  s.hmat_filename = [cf.hm_filename '_comp11.hmat'];
end

function MakeGifBig (cf, clim, giffn, is, varargin)
  function lim = PadLim(lim)
    if (numel(lim) < 2 || diff(lim) <= 0) lim = [0 1]; end
    lim = lim + 0.05*diff(lim)*[-1 1];
  end
  function Plot (sidx)
    clf;
    sp(611); plot(t, log10(mr), 'k-');
    set(gca, 'ylim', PadLim(log10([min(mr(mr > 0)) max(mr)])),...
             'xlim', PadLim(t([1 end])),...
             'xtick', [], 'ytick', []);
    axis off;
    sp(612); plot(ms(1:sidx) - ms(1), log10(mr(1:sidx)), 'k-');
    axis tight; axis off;
    A = s.fhm(o.transform(s.(o.fld)(:,sidx)));
    sp(3,1,2:3); imagesc(A); caxis(clim);
    colormap(jet(1024));
    axis equal; axis tight; axis xy; set(gca, 'xtick', [], 'ytick', []);
    xlim([611 1899]); ylim([100 424]);
    drawnow;
  end
  
  if (nargin < 4) is = []; end
  t = fdra2c('qload_t', cf);
  if (isempty(is)) is = 1:numel(t); end
  t = t(is)/(3600*24);
  o = popt(varargin, {{'delay' 0},
                      {'fld' 'v'},
                      {'transform' @log10}});
  nis = numel(is);
  fprintf('nis %d\n', nis);
  
  mr = zeros(size(t));
  ms = zeros(size(t));
  rid = dc3dm.mRead(cf.rmesh_filename);
  rs = dc3dm.mRects(rid);
  dc3dm.mClear(rid);
  w = rs(3,:).*rs(4,:);
  w = w/sum(w);
  
  clf;
  % Calibrate the colormap.
  s = fdra2c('qload', cf, is(1), 'do_stride', 0);
  sp(311); plot(1:10, -5:4, 'k-');
  sp(3,1,2:3);
  imagesc(clim(1) + (clim(2) - clim(1))*...
          rand(size(s.fh(o.transform(s.(o.fld)(:,1))))));
  warning('off', 'MATLAB:getframe:RequestedRectangleExceedsFigureBounds');
  f = getframe(gcf);
  [im map] = rgb2ind(f.cdata,1024);
  mem = 4*size(im,1)*size(im,2)*nis;
  %if (mem > 4e9) error(sprintf('Mem is %s: too big!', hbytes(mem))); end

  chunk = ceil(1e9/(4*size(s.(o.fld),1)));
  isi = 1;
  k = 1;
  first = true;
  im(:,:,1,nis) = 0;
  fprintf('%d: ', nis);
  while (isi < nis)
    up = min(isi + chunk - 1, nis);
    sis = is(isi:up);
    isi = isi + chunk;
    s = fdra2c('qload', cf, sis, 'do_stride', 0, 's', s);
    for (i = 1:numel(s.t))
      fprintf('%d ', k);
      ms(k) = sum(w(:).*s.slip(:,i));
      mr(k) = sum(w(:).*s.v(:,i));
      Plot(i);
      f = getframe(gcf);
      if (first)
	im = uint8(zeros(size(f.cdata,1), size(f.cdata,2)));
	first = false;
      end
      im(:,:,1,k) = rgb2ind(f.cdata, map);
      k = k + 1;
    end
  end
  fprintf(1,'\n');
  imwrite(im, map, giffn, 'DelayTime', o.delay, 'LoopCount', inf);
end

function bc = LoadBc (fn, nr)
  if (nargin < 2) nr = []; end
  fid = fopen(fn, 'r');
  bc = fread(fid, inf, 'double');
  if (isempty(nr))
    nc = 4;
  else
    assert(mod(numel(bc), nr) == 0);
    nc = numel(bc)/nr;
  end
  bc = reshape(bc, nr, nc);
  bc = sum(bc, 2);
  fclose(fid);
end

function [c rs] = rmi_Init (rid, varargin)
  o = popt(varargin, {{'use_max' 1}});
  rs = dc3dm.mRects(rid);
  d = dc3dm.mData(rs);
  if (o.use_max) dx = d.Dx; dy = d.Dy;
  else dx = d.dx; dy = d.dy; end
  c.x = (d.xlim(1) + 0.5*d.dx) : dx : d.xlim(2);
  c.y = (d.ylim(1) + 0.5*d.dy) : dy : d.ylim(2);
  [X Y] = meshgrid(c.x, c.y);
  c.id = dc3dm.mIds(rid, X(:), Y(:));
end

function A = rmi_Image (c, A)
  A = reshape(A(c.id), numel(c.y), numel(c.x));
end

% ------------------------------------------------------------------------------
% Image viewer.
%   s = fdra2c('qload', cf, ...);
%   fdra2c('vi_Start', s); fdra2c('vi_Start', s, 'fld', 'slip');

function vi_Start (s, varargin)
  o = popt(varargin, {{'fld' 'v'},
                      {'transform' @log10},
                      {'clim', []},
                      {'use_max' 1}});
  
  hd.ui.fig = figure();
  set(hd.ui.fig, 'name', sprintf('viewimg %s', s.fn), 'numbertitle', 'off');
  set(hd.ui.fig, 'KeyPressFcn', @vi_cb_fig_Kbd);
  hd.ui.clim.lo = uicontrol(...
    hd.ui.fig, 'style', 'edit', 'string', '',...
    'tooltip', 'Set caxis; leave blank for auto setting.',...
    'unit', 'normalized', 'position', [0.8 0 0.1 0.05],...
    'callback', @vi_cb_caxis);
  hd.ui.clim.hi = uicontrol(...
    hd.ui.fig, 'style', 'edit', 'string', '',...
    'tooltip', 'Set caxis; leave blank for auto setting.',...
    'unit', 'normalized', 'position', [0.9 0 0.1 0.05],...
    'callback', @vi_cb_caxis);
  hd.ui.msg = uicontrol(...
    hd.ui.fig, 'style', 'text', 'string', '',...
    'tooltip', 'Info', 'unit', 'normalized', 'position', [0 0.95 1 0.05],...
    'horizontalalignment', 'left');
  x = 0; dx = 0.05;
  hd.ui.left = uicontrol(...
    hd.ui.fig, 'style', 'pushbutton', 'string', '<',...
    'tooltip', 'Flip backward.',...
    'unit', 'normalized', 'position', [x 0 dx 0.05],...
    'callback', @(ho,ed)vi_cb_pb_Flip(ho,ed,-1));
  x = x + dx;
  hd.ui.setidx = uicontrol(...
    hd.ui.fig, 'style', 'edit', 'string', '1',...
    'tooltip', 'Set image index.',...
    'unit', 'normalized', 'position', [x 0 2*dx 0.05],...
    'callback', @vi_cb_ed_Flip);
  x = x + 2*dx;
  hd.ui.setstride = uicontrol(...
    hd.ui.fig, 'style', 'edit', 'string', '1',...
    'tooltip', 'Set stride.',...
    'unit', 'normalized', 'position', [x 0 dx 0.05],...
    'callback', @vi_cb_ed_Stride);
  x = x + dx;
  hd.ui.right = uicontrol(...
    hd.ui.fig, 'style', 'pushbutton', 'string', '>',...
    'tooltip', 'Flip forward.',...
    'unit', 'normalized', 'position', [x 0 dx 0.05],...
    'callback', @(ho,ed)vi_cb_pb_Flip(ho,ed,1));
  x = x + dx;
  hd.ui.play = uicontrol(...
    hd.ui.fig, 'style', 'pushbutton', 'string', 'Play',...
    'tooltip', 'Play as movie.', ...
    'unit', 'normalized', 'position', [x 0 2*dx 0.05],...
    'callback', @vi_cb_pb_Play);
  x = x + 2*dx;
  hd.ui.stop = uicontrol(...
    hd.ui.fig, 'style', 'pushbutton', 'string', 'Stop',...
    'tooltip', 'Stop movie.',...
    'unit', 'normalized', 'position', [x 0 2*dx 0.05],...
    'callback', @vi_cb_pb_Stop, 'buttondownfcn', @vi_cb_pb_Stop);
  x = x + 2*dx;
  hd.ui.zoom = uicontrol(...
    hd.ui.fig, 'style', 'pushbutton', 'string', 'Zoom',...
    'tooltip', 'Toggle zoom.',...
    'unit', 'normalized', 'position', [x 0 2*dx 0.05],...
    'callback', @vi_cb_pb_Zoom);
  x = x + 2*dx;
  hd.ui.indices = uicontrol(...
    hd.ui.fig, 'style', 'pushbutton', 'string', 'Indices',...
    'tooltip', 'Print indices visible in top plot.',...
    'unit', 'normalized', 'position', [x 0 2*dx 0.05],...
    'callback', @vi_cb_pb_PrintIndices);
  x = x + 2*dx;
  hd.ui.freehand = uicontrol(...
    hd.ui.fig, 'style', 'pushbutton', 'string', 'Freehand',...
    'tooltip', 'Create a freehand shape. Record data in global g_fh.',...
    'unit', 'normalized', 'position', [x 0 2*dx 0.05],...
    'callback', @vi_cb_pb_Freehand);
  global g_fh; g_fh = [];
  x = x + 2*dx;

  hd.s = s;
  hd.t.i = 1;
  hd.t.stride = 1;
  hd.t.transform = o.transform;
  if (isempty(o.clim))
    hd.t.ca = [min(o.transform(s.(o.fld)(:))) max(o.transform(s.(o.fld)(:)))];
  else
    hd.t.ca = o.clim;
  end
  if (strcmp(o.fld, 'v') && hd.t.ca(1) < -15) hd.t.ca(1) = -15; end
  hd.t.oca = hd.t.ca;
  hd.t.fld = o.fld;
  hd = vi_InitImage(s, hd, o);
  vi_ShowImage(s, hd, hd.t.i);
  
  set(hd.ui.clim.lo, 'string', sprintf('%3.1f', hd.t.ca(1)));
  set(hd.ui.clim.hi, 'string', sprintf('%3.1f', hd.t.ca(2)));
  set(hd.ui.fig, 'handlevisibility', 'callback', 'toolbar', 'figure');
  guidata(hd.ui.fig, hd);
end

function vi_cb_fig_Kbd (ho, ed)
  switch (ed.Character)
    case {'g' 'j' 'p'}
      delta = -1;
    case {'f' 'h' 'n'}
      delta = 1;
    otherwise
      return;
  end
  vi_Flip(ho, delta);
end

function vi_Flip (ho, delta)
  hd = guidata(ho);
  hd.t.i = mod(hd.t.i - 1 + delta*hd.t.stride, numel(hd.s.t)) + 1;
  vi_ShowImage(hd.s, hd, hd.t.i);
  try
    guidata(hd.ui.fig, hd);
  catch
    % Don't do anything. I just want to silence the error that results from
    % closing the window while it's still flipping images.
  end
end

function vi_cb_pb_Flip (ho, ed, delta)
  vi_Flip(ho, delta);
end

function vi_cb_ed_Flip (ho, ed)
  hd = guidata(ho);
  n = str2num(get(hd.ui.setidx, 'string'));
  if (numel(n) == 1 && n >= 1 && n <= numel(hd.s.t))
    hd.t.i = n;
    guidata(ho, hd);
  end
  vi_ShowImage(hd.s, hd, hd.t.i);
end

function vi_cb_ed_Stride (ho, ed)
  hd = guidata(ho);
  n = str2num(get(hd.ui.setstride, 'string'));
  if (numel(n) == 1 && round(n) == n && n >= 1 && n <= numel(hd.s.t))
    hd.t.stride = n;
    guidata(ho, hd);
  else
    set(hd.ui.setstride, 'string', num2str(hd.t.stride));
  end
end

function vi_cb_pb_Play (ho, ed)
  hd = guidata(ho);
  hd.stop = 0;
  guidata(ho, hd);
  while (1)
    % Stop button pushed?
    try
      hd = guidata(ho);
    catch
      % Almost certainly here because the window was closed.
      return;
    end
    if (hd.stop) return; end
    % Nope, so continue.
    vi_Flip(ho, 1);
    drawnow;
    pause(0.1);
  end
end

function vi_cb_pb_Stop (ho, ed)
  hd = guidata(ho);
  hd.stop = 1;
  guidata(ho, hd);
end

function vi_cb_pb_Zoom (ho, ed)
  zoom;
end

function vi_cb_pb_PrintIndices (ho, ed)
  hd = guidata(ho);
  subplot(311);
  xl = get(gca, 'xlim');
  yl = get(gca, 'ylim');
  mask = log10(hd.t.mr) >= yl(1) & log10(hd.t.mr) <= yl(2);
  mr_lo = find(hd.t.t >= xl(1) & mask, 1);
  mr_hi = find(hd.t.t <= xl(2) & mask, 1, 'last');
  str = get(hd.ui.msg, 'string');
  is = strfind(str, '|');
  if (numel(is) == 2) str = str(1:is(2)-2); end
  set(hd.ui.msg, 'string', sprintf(...
    '%s | %d to %d',...
    str, mr_lo, mr_hi));
  % On the command line for convenience:
  fprintf('(%d:%d)\n', mr_lo, mr_hi);
end

function vi_cb_pb_Freehand (ho, ed)
  subplot(3,1,2:3);
  h = impoly('closed', false);
  p = h.getPosition();
  hd = guidata(ho);
  global g_fh;
  if (isempty(g_fh)) g_fh.ps = {}; end
  if (numel(g_fh.ps) < hd.t.i || isempty(g_fh.ps{hd.t.i}))
    g_fh.ps{hd.t.i} = {};
  end
  g_fh.ps{hd.t.i}{end+1} = p;
end

function vi_cb_caxis (ho, ed)
  hd = guidata(ho);
  if (isempty(get(hd.ui.clim.lo, 'string')) ||...
      isempty(get(hd.ui.clim.hi, 'string')))
    cl = hd.t.oca;
  else
    cl = [str2num(get(hd.ui.clim.lo, 'string'))
          str2num(get(hd.ui.clim.hi, 'string'))];
  end
  if (length(cl) == 2 && cl(1) < cl(2))
    subplot(3,1,2:3); caxis(cl);
    hd.t.ca = cl;
    guidata(hd.ui.fig, hd);
  else
    cl = hd.t.ca;
  end
  set(hd.ui.clim.lo, 'string', sprintf('%3.1f', cl(1)));
  set(hd.ui.clim.hi, 'string', sprintf('%3.1f', cl(2)));
end

function hd = vi_InitImage (s, hd, o)
  hd.t.rmi = fdra2c('rmi_Init', s.rmesh_filename, 'use_max', o.use_max);
  hd.t.t = s.t/(3600*24);
  rs = dc3dm.mRects(s.rmesh_filename);
  w = rs(3,:).*rs(4,:);
  hd.t.mr = sum(repmat(w(:), 1, size(s.v, 2)) .* s.v)/sum(w);
end

function vi_ShowImage (s, hd, idx)
  subplot(3,1,2:3);
  xl = xlim(); yl = ylim();
  I = fdra2c('rmi_Image', hd.t.rmi, s.(hd.t.fld)(:, idx));
  imagesc(hd.t.rmi.x*1e-3, hd.t.rmi.y*1e-3, hd.t.transform(I));
  caxis(hd.t.ca);
  tb = 'NA'; ta = 'NA';
  if (idx > 1) tb = ToTimeStr(s.t(idx) - s.t(idx-1)); end
  if (idx < numel(s.t)) ta = ToTimeStr(s.t(idx+1) - s.t(idx)); end
  axis xy; axis equal; axis tight;
  if (any(xl ~= [0 1])) xlim(xl); ylim(yl); end
  colorbar;
  colormap(jet(1024));
  
  subplot(311);
  plot(hd.t.t, log10(hd.t.mr), 'b.-', hd.t.t(idx), log10(hd.t.mr(idx)), 'ro');
  xlabel('Time [day]'); ylabel('log_{10} v [m/s]');
  axis tight;
    
  set(hd.ui.setidx, 'string', num2str(hd.t.i));

  set(hd.ui.msg, 'string', sprintf(...
    '%d of %d | %s after, %s before',...
    idx, numel(hd.s.t), tb, ta));
  drawnow;
end

function t = ToTimeStr (t)
  if (t <= 0.1) t = sprintf('%1.1e [s]', t);
  elseif (t <= 1000) t = sprintf('%1.1f [s]', t);
  elseif (t <= 24*3600) t = sprintf('%1.1f [hr]', t/3600);
  elseif (t <= 24*3600*365.25) t = sprintf('%1.1f [day]', t/(24*3600));
  else sprintf('%1.1f [yr]', t/(24*3600*365.25)); end
end

function u = TimeUnit (ustr)
  switch (lower(ustr(1)))
    case 's'
      u.fac = 1;
      u.str = 'sec';
    case 'h'
      u.fac = 1/3600;
      u.str = 'hr';
    case 'd'
      u.fac = 1/(24*3600);
      u.str = 'day';
    case 'y'
      u.fac = 1/(24*3600*365.25);
      u.str = 'yr';
  end
end

function Imagesc (s, v)
% Use as
%   img = @(s,v) fdr2c('Imagesc',s,v);
%   imv = @(s) fdr2c('Imagesc',s,log10(s.v));
  [v x] = Interp(s, v);
  imagesc(1:len(s.t), x*1e-3, v);
end

% ------------------------------------------------------------------------------
% Analysis.

function is = TrimTs (t, min_stride, max_stride, min_dt, max_dt)
  is = zeros(numel(t), 1);
  n = 1;
  is(1) = 1;
  for (i = 2:numel(t)-1)
    include = (i - is(n) == max_stride || t(i+1) - t(is(n)) >= max_dt) &&...
              t(i) - t(is(n)) >= min_dt && i - is(n) >= min_stride;
    if (include) n = n + 1; is(n) = i; end
  end
  n = n + 1; is(n) = numel(t);
  is = is(1:n);
end

function [um rm] = SizeFactor (rmesh_filename)
  rs = dc3dm.mRects(rmesh_filename);
  nc = rs(3,:)/min(rs(3,:));
  um = sum(nc.^2);
  rm = numel(nc);
  pr('um: %d rm: %d fac %1.1f\n', um, rm, um/rm);
end

function [hsb hsbma] = Hstar (d_c, s_normal, a, b)
  mu_hat = 3e10/(1 - 0.25);
  hsb = mu_hat*d_c ./ (s_normal.*b);
  hsbma = pi*mu_hat*d_c ./ (4*s_normal.*(b - a)); hsbma(a >= b) = nan;
end

function lcodx = GetTrueLcodx (cf)
  hstar = Hstar(cf.d_c, cf.s_normal, cf.a, cf.b);
  rs = dc3dm.mRects(cf.rmesh_filename);
  lcodx = hstar(:)./rs(3,:)';
end

function c = ixy_Compute (cf, is, varargin)
  o = popt(varargin, {{'c'}});
  if (nargin < 2) is = []; end
  if (isempty(is) || numel(is) == 1)
    t = qload_t(cf);
    if (isempty(is))
      is = 1:numel(t);
    else
      is = 1:is:numel(t);
    end
  end

  if (isempty(o.c))
    c = ixy_Init(cf, varargin{:});
  else
    c = o.c;
  end
  c.is = is;
  c.Ix = []; c.Iy = []; c.t = [];
  
  nis = numel(is);
  chunk = ceil(1e8/(4*c.nrs));
  pr('%d: ', nis);
  k = 0;
  while (1)
    sis = k + (1:chunk);
    sis(sis > numel(is)) = [];
    
    iss = is(sis);
    if (~isempty(o.c))
      is_new = setdiff(iss, o.c.is);
      [is_old p_old] = intersect(o.c.is, iss);
    else
      is_new = iss;
    end

    if (~isempty(is_new))
      [Ix Iy t] = ixy_IntegrateImages(c, is_new, varargin{:});
    else
      Ix = []; Iy = []; t = [];
    end
    
    if (~isempty(o.c))
      isa = [is_old is_new];
      [iss p] = sort(isa);
      Ix = [o.c.Ix(:,p_old), Ix]; Ix = Ix(:,p);
      Iy = [o.c.Iy(:,p_old), Iy]; Iy = Iy(:,p);
      t = [o.c.t(p_old), t]; t = t(p);
    end
    c.Ix = [c.Ix Ix];
    c.Iy = [c.Iy Iy];
    c.t = [c.t t];
    
    k = k + chunk;
    pr('%d ', sis(end));
    if (sis(end) == numel(is)) break; end
  end
  pr('\n');
end

function [h Ix Iy] = ixy_ShowChunk (c, is, varargin)
  if (nargin < 2 || isempty(is)) is = 1:numel(c.t); end
  o = popt(varargin, {{'n' numel(is)}, {'unit' 'day'}});
  unit = TimeUnit(o.unit);
  tlin = linspace(c.t(is(1)), c.t(is(end)), o.n);
  Ix = interp1(c.t(is), c.Ix(:, is).', tlin).';
  Iy = interp1(c.t(is), c.Iy(:, is).', tlin).';
  iml = @(y, I) imagesc((tlin - tlin(1))*unit.fac, y, log10(I));
  h(1) = sp(211); iml(c.x*1e-3, Ix); cb; axis xy;
  ylabel('Along strike [km]');
  h(2) = sp(212); iml(c.y*1e-3, Iy); cb; axis xy;
  ylabel('Along dip [km]'); xlabel(sprintf('Time [%s]', unit.str));
  h = linkprop(h, 'xlim');
end

function h = ixy_ShowRaw (c, is, varargin)
  o = popt(varargin, {{'transform' @log10}});
  if (isempty(o.transform)) o.transform = @(x) x; end
  if (nargin < 2 || isempty(is)) is = 1:numel(c.t); end
  img = @(x, I) imagesc(1:numel(is), x*1e-3, o.transform(I));
  h(1) = sp(211); img(c.x, c.Ix); cb; axis xy;
  ylabel('along strike [km]');
  h(2) = sp(212); img(c.y, c.Iy); cb; axis xy;
  ylabel('along dip [km]'); xlabel('time step');
  h = linkprop(h, 'xlim');
end

function c = ixy_Init (cf, varargin)
  function v = ixy_default_fn (I, X, Y)
    v = mean(I);
  end
  c = popt(varargin, {{'xlim'}, {'ylim'}, {'nx'}, {'ny'}, ...
                      {'fn' @ixy_default_fn}});
  c.cf = cf;

  rs = dc3dm.mRects(cf.rmesh_filename);
  c.nrs = size(rs, 2);
  d = dc3dm.mData(rs);
  
  if (isempty(c.xlim)) c.xlim = d.xlim; end
  if (isempty(c.ylim)) c.ylim = d.ylim; end
  if (isempty(c.nx)) c.nx = round(diff(c.xlim)/d.dx); end
  if (isempty(c.ny)) c.ny = round(diff(c.ylim)/d.dy); end

  c.x = CC(linspace(c.xlim(1), c.xlim(2), c.nx + 1));
  c.y = CC(linspace(c.ylim(1), c.ylim(2), c.ny + 1));
  [c.X c.Y] = meshgrid(c.x, c.y);
  c.id = dc3dm.mIds(cf.rmesh_filename, c.X, c.Y);
end

function [Ix Iy t] = ixy_IntegrateImages (c, is, varargin)
  o = popt(varargin, {{'fld' 'v'}});
  s = fdra2c('qload', c.cf, is, 'do_stride', 0);
  t = s.t;
  ni = numel(s.t);
  Ix = zeros(numel(c.x), ni);
  Iy = zeros(numel(c.y), ni);
  for (i = 1:ni)
    I = s.(o.fld)(c.id, i);
    I = reshape(I, size(c.id));
    Ix(:,i) = c.fn(I, c.X, c.Y).';
    Iy(:,i) = c.fn(I.', c.X.', c.Y.').';
  end
end

function ixy_AutoShow (c, varargin)
% This is one way to use c from ixy_Compute, but I'm not at all happy with
% it. Going to try something more manual next.
  o = popt(varargin, {{'short' 1}, {'n'}});
  ts = GroupTimes(c.t, o.short);
  if (isempty(o.n)) o.n = ceil(numel(c.t)/numel(ts)); end
  Ix = []; Iy = [];
  cuts = 0;
  k = 0;
  for (i = 1:numel(ts))
    is = k + (1:numel(ts{i}));
    if (numel(ts{i}) == 1)
      Ix = [Ix c.Ix(:, is)];
      Iy = [Iy c.Iy(:, is)];
    else
      tlin = linspace(ts{i}(1), ts{i}(end), o.n);
      if (1) % Hack for when we stored t as float.
        dt = double(diff(ts{i}));
        dt(dt == 0) = 10*eps*double(ts{i}(1));
        ts{i} = double(ts{i}(1)) + [0 cumsum(dt)];
      end
      Ix = [Ix interp1(ts{i}, c.Ix(:, is).', tlin).'];
      Iy = [Iy interp1(ts{i}, c.Iy(:, is).', tlin).'];
    end
    k = k + numel(ts{i});
    cuts(end+1) = k;
  end
  iml = @(y, I) imagesc(1:size(Ix,2), y, log10(I));
  sp(211); iml(c.x*1e-3, Ix); cb;
  sp(212); iml(c.y*1e-3, Iy); cb;
end

function ts = GroupTimes (t, short)
% Group times into short-dt and long-dt sections.
  assert(numel(t) >= 2);
  ts = {t(1:2)};
  if (t(2) - t(1) <= short) mode = 1; else mode = 2; end
  for (i = 3:numel(t))
    dt = t(i) - t(i-1);
    if (mode == 1)
      if (dt <= short)
        ts{end}(end+1) = t(i);
      else
        ts{end+1} = [ts{end}(end) t(i)];
        ts{end-1}(end) = [];
        if (isempty(ts{end-1})) ts = {ts{1:end-2} ts{end}}; end
        mode = 2;
      end
    else
      ts{end}(end+1) = t(i);
      if (dt <= short && i < numel(t))
        if (t(i+1) - t(i) <= short)
          mode = 1;
          ts{end+1} = [];
        end
      end
    end
  end
  % Check.
  tc = [];
  for (i = 1:numel(ts)) tc = [tc; ts{i}(:)]; end
  assert(sum(tc - t(:)) == 0);
end

function SubsampleSdfs (fn, is)
% Subsample all .tsdf and .sdf files associated with fn*.sdf to just the
% indices 'is'. If 'is' is a scalar, then it is a stride. NB. The original
% file is overwritten. That's the point of this routine.
  dr = fileparts(fn);
  d = dir([fn '*.sdf']);
  tmp_fn = [fn '.tmp'];
  
  reply = input(sprintf('Overwrite %s*.?sdf? Y/N [N]:', fn), 's');
  if (isempty(reply)) reply = 'n'; end
  if (lower(reply) ~= 'y')
    fprintf('Ok, returning.\n');
    return;
  end
  
  for (i = 1:numel(d))
    [~, name, suffix] = fileparts(d(i).name);
    fn = [dr filesep name];
    fprintf('%s\n', fn);
    ExtractSomeCols(fn, tmp_fn, is);
    movefile([tmp_fn '.sdf'], [fn '.sdf']);
    movefile([tmp_fn '.tsdf'], [fn '.tsdf']);
  end
end

function ExtractSomeCols (ofn, nfn, is)
  assert(~exist([nfn '.sdf'], 'file'));
  if (numel(is) == 1)
    t = SaveStreamData('Read', [ofn '.tsdf'], 1, 1, 1);
    is = 1:is:numel(t);
    if (is(end) < numel(t)) is(end+1) = numel(t); end
  end
  CopySdf([ofn '.sdf'], [nfn '.sdf'], is, 1);
  CopySdf([ofn '.tsdf'], [nfn '.tsdf'], is, 2);
end

function A = CopySdf (src, dst, cidxs, realp)
  d = SaveStreamData('Read', src, 1, [], 0);
  chunk = ceil(1e9/(6*numel(d)));
  ssd = SaveStreamData('Init', dst, 'realp', realp);
  while (1)
    if (numel(cidxs) > chunk)
      is = cidxs(1:chunk);
      cidxs = cidxs(chunk+1:end);
    elseif (isempty(cidxs))
      break;
    else
      is = cidxs;
      cidxs = [];
    end
    d = SaveStreamData('Read', src, is, [], 0);
    ssd = SaveStreamData('Write', ssd, d, 1);
  end
  SaveStreamData('Finalize', ssd);
end

function s_Map (cf, is, fn)
% Run fn on qload struct for chunks of indices is from fdra kvf cf.
%   fn has signature
%     fn('init', s, n), where s is for is(1) and n = numel(is);
%     fn('run', s), where s contains a chunk of is;
%     fn('done', s), where s is the final chunk that was already processed.
  assert(~isempty(is));
  chunk = ceil(1e9 / (8*numel(cf.a)));
  progress = 0; nis = numel(is);
  s = fdra2c('qload', cf, is(1), 'do_stride', 0);
  fn('init', s)
  while (1)
    if (numel(is) > chunk)
      cis = is(1:chunk);
      is = is(chunk+1:end);
    elseif (isempty(is))
      break;
    else
      cis = is;
      is = [];
    end
    s = fdra2c('qload', cf, cis, 'do_stride', 0, 's', s);
    fn('run', s);
    p = round(100*(nis - numel(is))/nis);
    if (p > progress) progress = p; pr('%d ', progress); end
  end
  fn('done', s);
  pr('\n');
end

function c = TraceSlip (cf, varargin)
  persistent q;
  function Fn (stat, s, varargin)
    switch (stat)
      case 'init'
        rs = dc3dm.mRects(s.rmesh_filename);
        [x y] = dc3dm.mCC(rs);
        uy = unique(y);
        y_want = uy(q.o.yis);
        q.is = [];
        for (i = 1:numel(y_want))
          yis = find(y == y_want(i));
          yis = yis(:).';
          q.is = [q.is yis(1:40:end)];
        end
        q.y = y(q.is);
        q.x = x(q.is);
        pr('#is %d\n', numel(q.is));
        q.slip = zeros(numel(q.is), numel(q.t));
        q.k = 0;
        
      case 'run'
        cis = q.k + (1:numel(s.t));
        q.slip(:, cis) = s.slip(q.is, :);
        q.k = q.k + numel(s.t);
        
      case 'done'
    end
  end
  
  q.o = popt(varargin, {{'yis' []}});
  if (isempty(q.o.yis)) error('No y indices specified.'); end
  q.t = fdra2c('qload_t', cf, 'slip');
  s_Map(cf, 1:numel(q.t), @Fn);
  c = struct('t', q.t, 'slip', q.slip, 'x', q.x, 'y', q.y, 'is', q.is);
  clear persistent q;
end

function c = GetDenseAvgV (cf, varargin)
  persistent q;
  function Fn (stat, s, varargin)
    switch (stat)
      case 'init'
        q.mr = zeros(size(q.t));
        q.k = 0;
      
      case 'run'
        cis = q.k + (1:numel(s.t));
        q.mr(cis) = sum(repmat(q.rs(3,:)', 1, size(s.v, 2)).*s.v) /...
            sum(q.rs(3,:));
        q.k = q.k + numel(s.t);
        
      case 'done'
    end
  end
  
  q.rs = dc3dm.mRects(cf.rmesh_filename);
  q.d = dc3dm.mData(q.rs);
  q.o = popt(varargin, {{'xlim' q.d.xlim}; {'ylim' q.d.ylim}});
  
  q.t = fdra2c('qload_t', cf, 'slip');
  s_Map(cf, 1:numel(q.t), @Fn);
  
  c = struct('t', q.t, 'mr', q.mr);
  clear persistent q;
end

function c = FindMaxV (cf, varargin)
  persistent q;
  function Fn (stat, s, varargin)
    switch (stat)
      case 'init'
        q.idx = [];
        q.v_max = [];
        
      case 'run'
        [v_max idx] = max(s.v);
        q.idx = [q.idx idx];
        q.v_max = [q.v_max v_max];
        
      case 'done'
    end
  end

  q = struct;
  rs = dc3dm.mRects(cf.rmesh_filename);
  [q.x q.y] = dc3dm.mCC(rs);
  q.md = dc3dm.mData(rs);
  
  q.t = fdra2c('qload_t', cf, 'v');
  s_Map(cf, 1:numel(q.t), @Fn);
  
  c = q;
  clear persistent q;
end

function ts_Plot (c, varargin)
  o = popt(varargin,...
           {{'nx', 5}});
  uy = unique(c.y);
  is = [];
  for (i = 1:numel(uy))
    iis = find(c.y(:).' == uy(i));
    stride = ceil(numel(iis) / o.nx);
    is = [is iis(1:stride:end)];
  end
  plot(s2y(c.t)*365.25, c.slip(is,:));
end

% ------------------------------------------------------------------------------
% Util

function c = CC (v)
% Cell-centered from node-centered.
  c = 0.5*(v(1:end-1) + v(2:end));
end
