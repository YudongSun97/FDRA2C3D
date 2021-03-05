function varargout = wrk(varargin)
  [varargout{1:nargout}] = feval(varargin{:});
end

function Addpath()
  addpath matlab;
end

function x = randspace(xlo, xhi, n, alpha)
  if (nargin < 4) alpha = 1; end
  x = cumsum(alpha + rand(1, n));
  x = (x - x(1))/(x(n) - x(1));
  x = xlo + (xhi - xlo)*x;
end

function c = WriteMeshRectKvf()
  d = 500;
  c.x = 10 + randspace(0, 0.9*d, 900);
  c.y = 20 + randspace(0, d, 1000);
  [X Y] = meshgrid(c.x, c.y);
  c.f = cos(5*X/d).^2 + sin(4*Y/d).^2;
  c.min_len = 1e-3;
  c.max_len = d;
  %c.f = c.f*c.max_len/max(c.f(:));
  c.f = c.f/max(c.f(:)); c.f = c.min_len + (c.max_len - c.min_len)*c.f;
  c.save_filename = 'test';
  kvf('Write', 'testmr.kvf', c, true);
end

function rm = ReadRectMesh(fn)
  rm.fn = fn;
  fid = fopen([fn '.rect'], 'r');
  if (fid == -1) error(sprintf('Can''t read %s\n', fn)); end
  rm.domain = fread(fid, 4, 'double');
  rm.ropts = fread(fid, 2, 'double');
  n = fread(fid, 1, 'int64');
  rm.rs = reshape(fread(fid, 4*n, 'double'), 4, n);
  fclose(fid);
end

function DrawRectMesh(rm, c)
  clf;
  fn = @(x)x;%@log10;
  
  sp(221); imagesc(c.x, c.y, fn(c.f)); cb; axis xy;
  caxis(fn(rm.ropts)); ca = caxis();
  
  rmesh('read', [rm.fn '.ser']);
  [X Y] = meshgrid(c.x, c.y);
  t = tic(); is = rmesh('getids', X, Y); et = toc(t);
  pr('getrects et: %1.3e\n', et);
  
  sp(222);
  img = reshape(max(rm.rs(3:4,is(:))), numel(c.y), numel(c.x));
  imagesc(fn(img)); cb; axis xy; caxis(ca);
  
  sp(223); imagesc(c.f < img & c.f >= rm.ropts(1));
  axis xy; title('Did we miss anything?');

  sp(224); imagesc(img - c.f);
  axis xy; cb;

  rmesh('free');
end

function scratch()
  c=wrk('WriteMeshRectKvf');tic;system('./bin/MakeRectUnstructMesh testmr.kvf');toc,rm=wrk('ReadRectMesh',c.save_filename);wrk('DrawRectMesh',rm,c);pr('nrect = %d\n',size(rm.rs,2));
end
