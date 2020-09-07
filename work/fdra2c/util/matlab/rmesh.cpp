#include <stdio.h>
#include <iostream>
#include <string>
#include "mex.h"
#include "util/matlab/MexUtil.hpp"
#include "util/include/RectMeshUnstruct.hpp"
using namespace std;
using namespace util;

static bool GetStringv(const mxArray* ms, vector<char>& s)
{
  int strlen = mxGetNumberOfElements(ms) + 1;
  s.resize(strlen);
  if (mxGetString(ms, &s[0], strlen) != 0) return false;
  return true;
}

static bool GetString(const mxArray* ms, string& s)
{
  vector<char> vs;
  if (!GetStringv(ms, vs)) return false;
  s = string(&vs[0]);
  return true;
}

// Clean up command.
static bool GetCommand(const mxArray* mcmd, string& fn)
{
  vector<char> vfn;
  if (!GetStringv(mcmd, vfn)) return false;
  for (int i = 0; i < vfn.size() - 1; i++) vfn[i] = tolower(vfn[i]);
  fn = string(&vfn[0]);
  return true;
}

static rmesh::RectMeshUnstruct* _rm = NULL;

void mexFunction(int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs)
{
  string fn;
  if (nrhs < 1 || !GetCommand(prhs[0], fn))
    mexErrMsgTxt("[...] = rmesh('cmd',...)");
  if (fn == "read") {
    if (nlhs != 0 || nrhs != 2) mexErrMsgTxt("rmesh('read', filename)");
    if (_rm) {
      delete _rm;
      _rm = NULL;
    }
    string filename;
    if (!GetString(prhs[1], filename))
      mexErrMsgTxt("Arg 2 should be a string.");
    try {
      _rm = new rmesh::RectMeshUnstruct(filename);
    } catch (const FileException& e) {
      mexErrMsgTxt(e.GetMsg().c_str());
    }
  } else if (fn == "free") {
    if (nlhs != 0 || nrhs != 1) mexErrMsgTxt("rmesh('free')");
    if (_rm) {
      delete _rm;
      _rm = NULL;
    }
  } else if (fn == "getrects") {
    if (!_rm) mexErrMsgTxt("No current RectMeshUnstruct.");
    if (nlhs != 1 || nrhs != 3 ||
        mxGetNumberOfElements(prhs[1]) != mxGetNumberOfElements(prhs[2]))
      mexErrMsgTxt("ids = rmesh('getrects', x, y) with numel(x) == numel(y)");
    int n = mxGetNumberOfElements(prhs[1]);
    double* px = mxGetPr(prhs[1]);
    double* py = mxGetPr(prhs[2]);
    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[1]), mxGetN(prhs[1]), mxREAL);
    double* pi = mxGetPr(plhs[0]);
    for (int i = 0; i < n; i++)
      pi[i] = (int) _rm->GetRectId(px[i], py[i]) + 1; // change to base-1 idx
  } else {
    mexErrMsgTxt((fn + "is not a command.").c_str());
  }
}
