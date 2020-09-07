function varargout = SaveStreamData(fn,varargin)
% Functions to save a stream of data to a file. The data should be
% column-oriented, so that one is building a wide matrix by appending
% columns.
% 
% ssd = SaveStreamData('Init',filename,[varargin])
%   Save a streaming matrix having nr rows to a designated file.
%     varargin (parameter-value pairs):
%     'nr': Number of rows in data. If not provided, it will be deduced in
%   Write. Note that this means the data matrix passed to write must be
%   sized correctly.
%     'bfr_ncalls': Store in memory for ncalls to Write before writing to
%   disk.
%     'realp': 1 (single) or 2 (double, default). Precision at which to store
%   data.
%
% ssd = SaveStreamData('Write',ssd,data,[force_write])
%   data is a matrix having ncols columns. If force_write is provided and is
% true, then force a write despite buffering.
%
% SaveStreamData('Finalize',ssd)
%   Write any buffered data and close the ssd struct.
%
% data = SaveStreamData('Read',filename,[cidxs],[ridxs],[isStride])
%   Read in the contents of a SSD file. If the optional argument cidxs is
%   provided, then data contains the columns listed in cidxs in ascending
%   order. cidxs can contain column numbers that exceed the number of rows in
%   the file. [ridxs] does the equivalent for rows. Set cidxs=[] if you want
%   to specifiy ridxs but not cidxs. If isStride=1 is given as the fourth
%   argument, cidxs is interpreted as a stride and should be a scalar. Set
%   ridxs=[] if you want to specifiy isStride but not ridxs.
%
% See the function 'Test' below for an example.

% AMB ambrad@cs.stanford.edu
  
% This fdra2c version does not use 'iee-le'. Until I get machine-independent
% C++ I/O, I can't do it here.

  [varargout{1:nargout}] = feval(fn,varargin{:});
end

function ssd = Init(fn,varargin)
  [nr bfr_ncalls realp] = process_options(...
      varargin,'nr',[],'bfr_ncalls',0,'realp',2);
  if (realp ~= 1 && realp ~= 2) error('realp is 1 or 2.'); end
  fid = fopen(fn,'r');
  ssd.fid = fopen(fn,'wb');
  ssd.nr = nr;
  ssd.realp = realp;
  switch (realp)
   case 1, ssd.realp_str = 'single';
   case 2, ssd.realp_str = 'double';
  end
  if (~isempty(nr)) WriteHeader(ssd.fid,nr,realp); end
  ssd.bfr.ncalls = bfr_ncalls;
  ssd.bfr.i = 1;
  ssd.bfr.data = [];
end

function WriteHeader(fid,nr,realp)
  fwrite(fid,nr,'int32');
  fwrite(fid,realp,'int32');
end

function ssd = Write(ssd,data,force)
  if (isempty(ssd.nr))
    ssd.nr = size(data,1);
    if (ftell(ssd.fid) == 0)
      WriteHeader(ssd.fid,ssd.nr,ssd.realp);
    end
  end
  if (ssd.bfr.ncalls)
    if (nargin < 3) force = 0; end
    write = force;
    ssd.bfr.data = [ssd.bfr.data data];
    if (write)
      ssd.bfr.i = 1;
    else
      if (ssd.bfr.i == ssd.bfr.ncalls)
	write = 1;
	ssd.bfr.i = 1;
      else
	ssd.bfr.i = ssd.bfr.i + 1;
      end
    end
    if (write)
      fwrite(ssd.fid,ssd.bfr.data(:),ssd.realp_str);
      ssd.bfr.data = [];
    end
  else
    fwrite(ssd.fid,data(:),ssd.realp_str);
  end
end

function Finalize(ssd)
  if (ssd.bfr.ncalls) ssd = Write(ssd,[],1); end
  fclose(ssd.fid);
end

function nr = NbrRows(fn)
  fid = Fopen(fn);
  nr = fread(fid,1,'int32');
  fclose(fid);
end

function fid = Fopen(fn)
  fid = fopen(fn,'rb');
  if (fid == -1)
    fid = fopen([fn '.dat'],'rb');
  end
  if (fid == -1) error(sprintf('Can''t read %s',fn)); end
end

function A = Read(fn,cidxs,ridxs,doStride)
  if (nargin < 4) doStride = 0; end
  if (nargin < 3) ridxs = []; end
  if (nargin < 2) cidxs = []; end
  doStride = doStride || isempty(cidxs);
  fid = Fopen(fn);
  A = ReadFid(fid,cidxs,ridxs,doStride);
  fclose(fid);
end

function A = ReadFid(fid,cidxs,ridxs,doStride)
  if (isempty(cidxs) && isempty(ridxs))
    % Read the whole file
    [nr realp_str] = ReadHeader(fid);
    A  = fread(fid,inf,realp_str);
    nc = length(A)/nr;
    A = reshape(A,nr,nc);
  elseif (doStride)
    % Read cols with a stride and possibly selected rows
    [nr realp_str nbytes real_class] = ReadHeader(fid);
    if (isempty(ridxs)) ridxs = 1:nr; end
    stride = cidxs;
    if (length(stride) > 1) error('stride is a scalar'); end
    if (isempty(stride)) stride = 1; end
    chunk = min(1000,ceil(5e7/length(ridxs)));
    A = zeros(length(ridxs),chunk,real_class);
    cnt = 1;
    while (true)
      data = fread(fid,nr,realp_str);
      % Check for EOF
      if (length(data) ~= nr)
	% Resize A if necessary
	A = A(:,1:cnt-1);
	break;
      end
      % Alloc more for A?
      if (size(A,2) == cnt)
	A = [A A];
      end
      % Save this column
      A(:,cnt) = data(ridxs);
      cnt = cnt + 1;
      if (stride > 1)
	% Skip ahead
	stat = fseek(fid,nbytes*(stride-1)*nr,'cof');
	if (stat < 0)
	  A = A(:,1:cnt-1);
	  break;
	end
      end
    end
  else
    % Read selected cols and possibly rows
    [nr realp_str nbytes real_class] = ReadHeader(fid);
    if (isempty(ridxs)) snr = nr; else snr = length(ridxs); end
    ni = length(cidxs);
    A = zeros(snr,ni,real_class);
    d = 0;
    for(i = 1:ni)
      % Nbr of cols between next desired one and the one we're reading now
      d = cidxs(i) - d - 1;
      if (d > 0)
	% Skip ahead.
	stat = fseek(fid,nbytes*d*nr,'cof');
	if (stat < 0)
	  A = A(:,1:i-1);
	  break;
	end 
      end
      if (ni == 1)
	% Halve the amount of memory used in this special case
	A = fread(fid,nr,realp_str);
	if (~isempty(ridxs)) A = A(ridxs); end
      else
        if (isempty(ridxs))
          try
            A(:,i) = fread(fid,nr,realp_str);
          catch
            A = A(:,1:i-1);
            break;
          end
        else
          data = fread(fid,nr,realp_str);
          if (length(data) ~= nr)
            A = A(:,1:i-1);
            break;
          end
          A(:,i) = data(ridxs);
        end
      end
      d = cidxs(i);
    end
  end
end

function [nr realp_str nbytes real_class] = ReadHeader(fid)
  nr = fread(fid,1,'int32');
  realp = fread(fid,1,'int32');
  % For backwards compatibility. Can fail if the four bytes happen to look
  % exactly like an int32 1 or 2.
  if (realp ~= 1 && realp ~= 2)
    % We have an old file. Move back to just after the nr entry.
    fseek(fid,4,'bof');
    realp_str = 'double=>double';
    nbytes = 8;
    real_class = 'double';
  else
    switch (realp)
     case 1
      realp_str = 'single=>single';
      nbytes = 4;
      real_class = 'single';
     case 2
      realp_str = 'double=>double';
      nbytes = 8;
      real_class = 'double';
    end
  end
end

function Test
  fn = 'ssdtest.dat';
  ssd = SaveStreamData('Init',fn);
  TestSsd(ssd,fn);
  ssd = SaveStreamData('Init',fn,'realp',1);
  TestSsd(ssd,fn);
end

function TestSsd(ssd,fn)
  A0 = (1:3)'*(1:4);
  for (i = 1:size(A0,2))
    ssd = SaveStreamData('Write',ssd,A0(:,i));
  end
  SaveStreamData('Finalize',ssd);
  system(['ls -ltrh ' fn]);
  A1 = SaveStreamData('Read',fn);
  A2 = SaveStreamData('Read',fn,2:2:100);
  A3 = SaveStreamData('Read',fn,2,[],true);
  A0,A1,A2,A3
end
