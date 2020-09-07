function varargout = ode(fn,varargin)
% ODE. Explicit 1-2 and 2-3 (same as Matlab's) pairs. ode('Pair23') is meant
% to replace Matlab's ode23 function when a large problem is solved. It does
% not save more in memory than necessary. A user-supplied function is given
% the solution at each time step to save by the method desired.
%
% o = ode('Options',...): Same as 'odeset' except that 'OutputFcn' has the
% additional argument 'dt':
%     stop = OutputFcn(t,dt,y,msg,...).
%
% ode('Pair??',odeFcn,tspan,y0,opts,...): '??' is either '12' or '23'. odeFcn
% has the signature
%     [yd sigerr yset] = odeFcn(t,dt,y,...).
%   sigerr is nominally set to 0; set it to a non-0 value to signal an error
% and request a smaller time step.
%   yset resets the value of y to yset(is), where is = find(~isnan(yset)).
%   dt is passed to odeFcn, unlike in ode23, so that the user has access to
% the exact time difference between stages.
  
% AMB ambrad@cs.stanford.edu

  [varargout{1:nargout}] = feval(fn,varargin{:});
    
function o = Options(varargin)
% Can use this or just call odeset.
  [o.RelTol o.AbsTol o.OutputFcn o.InitialStep] =...
      process_options(varargin,'RelTol',1e-3,'AbsTol',1e-6,...
		      'OutputFcn',[],'InitialStep',1e-3);
  
function Pair12(odeFcn,tspan,y0,opts,varargin)
% 2nd-order method: explicit trapezoidal rule, also called Heun's method. Error
% control is modelled on Matlab's ode*.

  pow = 1/2;
  threshold = opts.AbsTol(:)/opts.RelTol;
  tdir = sign(diff(tspan));
  n = length(y0);
  
  hof = ~isempty(opts.OutputFcn);
  if(hof) of = opts.OutputFcn; end
  
  have_yset = nargout(odeFcn) >= 3;

  f = zeros(n,2);
  t = tspan(1);
  absh = opts.InitialStep;
  [f(:,1) serr] = feval(odeFcn,t,0,y0,varargin{:});
  y = y0;
  
  if(hof)
    feval(of,tspan,0,y,'init',varargin{:});
  end
  
  while(t < tspan(end)) % Loop to integrate
    hmin = 16*eps(t);
    nofailed = true;
    while(1) % Loop for one step
      h = tdir*absh;
      if(t + h > tspan(end))
	h = tspan(end) - t;
	absh = abs(h);
      end
      tnew = t + h;
      h = tnew - t;
    
      % Explicit Euler
      ye = y + h*f(:,1);
      [f(:,2) serr] = feval(odeFcn,tnew,h,ye,varargin{:});
      if(serr)
	nofailed = false;
	err = inf;
      else
	% Explicit trapezoidal rule
	yt = y + f*([1; 1]*h/2);
	% Relative difference between Euler and trapezoidal rule
	err = absh/2*norm((f*[1; -1])./...
			  max(max(abs(y),abs(yt)),threshold),inf);
      end

      if(err > opts.RelTol)
	if(absh <= hmin)
	  if(hof)
	    feval(of,t,h,y,'tolfail',varargin{:});
	  end
	  warning(sprintf('ode: Failure at t=%e. Unable to meet tol.\n',t));
	  return;
	end
	
	if(nofailed)
	  nofailed = false;
	  absh = max(hmin,absh*max(0.5,0.8*(opts.RelTol/err)^pow));
	else
	  absh = max(hmin,0.5*absh);
	end
      else
	if(nofailed)
	  fac = 0.8*(opts.RelTol/err)^pow;
	  absh = min(5,fac)*absh;
	end
	
	if(have_yset)
	  [f(:,1) serr odeFcnSet] = feval(odeFcn,tnew,h,yt,varargin{:});
	else
	  [f(:,1) serr] = feval(odeFcn,tnew,h,yt,varargin{:});
	end
	if(serr)
	  nofailed = false;
	  absh = max(hmin,0.5*absh);
	  continue;
	else
	  t = tnew;
	  if(have_yset)
	    yset_idxs = find(~isnan(odeFcnSet));
	    yt(yset_idxs) = odeFcnSet(yset_idxs);
	  end
	  y = yt;
	  if(hof)
	    stop = feval(of,t,h,y,'',varargin{:});
	    if(stop) return; end
	  end
	  break;
	end
      end
    end
  end

function Pair23(odeFcn,tspan,y0,opts,varargin)
% Clone of Matlab's ode23. Cite
%   P. Bogacki, L.F. Shampine, "A 3(2) Pair of Runge-Kutta Formulas",
%   Appl. Math Lett. 2(4), 321-325, 1989.

  threshold = opts.AbsTol(:)/opts.RelTol;
  tdir = sign(diff(tspan));
  n = length(y0);
  
  hof = ~isempty(opts.OutputFcn);
  if(hof) of = opts.OutputFcn; end
  
  have_yset = nargout(odeFcn) >= 3;

  f = zeros(n,4);
  t = tspan(1);
  absh = opts.InitialStep;
  [f(:,1) serr] = feval(odeFcn,t,0,y0,varargin{:});
  y = y0;
  
  if(hof)
    feval(of,tspan,0,y,'init',varargin{:});
  end
  
  % Taken from ode23.m
  pow = 1/3;
  A = [1/2, 3/4, 1];
  B = [1/2         0               2/9
       0           3/4             1/3
       0           0               4/9
       0           0               0  ]; 
  E = [-5/72; 1/12; 1/9; -1/8];
  
  while(t < tspan(end)) % Loop to integrate
    hmin = 16*eps(t);
    nofailed = true;
    while(1) % Loop for one step
      h = tdir*absh;
      if(t + h > tspan(end))
	h = tspan(end) - t;
	absh = abs(h);
      end
      tnew = t + h;
      h = tnew - t;
    
      hB = h*B;
      [f(:,2) serr] = feval(odeFcn,t+h*A(1),h*A(1),y+f*hB(:,1),varargin{:});
      if(~serr)
	[f(:,3) serr] = feval(odeFcn,t+h*A(2),h*A(2),y+f*hB(:,2),varargin{:});
	if(~serr)
	  yt = y + f*hB(:,3);
	  if(have_yset)
	    [f(:,4) serr odeFcnSet] = feval(odeFcn,tnew,h,yt,varargin{:});
	  else
	    [f(:,4) serr] = feval(odeFcn,tnew,h,yt,varargin{:});
	  end
	  if(~serr)
	    if(have_yset)
	      yset_idxs = find(~isnan(odeFcnSet));
	      yt(yset_idxs) = odeFcnSet(yset_idxs);
	    end
	    % h*f*E is the difference between the 3rd- and 2nd-order methods
	    err = absh*norm((f*E)./max(max(abs(y),abs(yt)),threshold),inf);
	  end
	end
      end
      if(serr)
	err = inf;
	nofailed = false;
      end
	    
      if(err > opts.RelTol)
	if(absh <= hmin)
	  if(hof)
	    feval(of,t,h,y,'tolfail',varargin{:});
	  end
	  warning(sprintf('ode: Failure at t=%e. Unable to meet tol.\n',t));
	  return;
	end
	
	if(nofailed)
	  nofailed = false;
	  absh = max(hmin,absh*max(0.5,0.8*(opts.RelTol/err)^pow));
	else
	  absh = max(hmin,0.5*absh);
	end
      else
	break;
      end
    end

    if(nofailed)
      fac = 0.8*(opts.RelTol/err)^pow;
      absh = min(5,fac)*absh;
    end
    
    t = tnew;
    y = yt;
    f(:,1) = f(:,4);
    if(hof)
      stop = feval(of,t,h,y,'',varargin{:});
      if(stop) return; end
    end
  end
