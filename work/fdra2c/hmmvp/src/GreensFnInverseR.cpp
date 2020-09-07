class InverseRGreensFn : public ImplGreensFn {
public:
  virtual void Init(const KeyValueFile* kvf) throw (Exception);
  virtual Hd* ComputeHd() { return NewHd(_x); }
  virtual bool Call(const vector<uint>& rs, const vector<uint>& cs,
                    double* B);

private:
  Matd _x;
  uint _order;
  double _delta;

  double Eval(uint i, uint j);
};

inline double InverseRGreensFn::Eval(uint i, uint j)
{
  // Distance between points i and j.
  double r = 0.0;
  for (uint k = 1; k <= 3; k++) {
    double d = _x(k, i) - _x(k, j);
    r += d*d;
  }
  r = sqrt(r + _delta);

  // If there is no softening, then I could return inf, but I prefer 0.
  if (r == 0.0 && _delta == 0.0) return 0.0;

  if (_order > 0)
    return 1.0 / pow(r, _order);
  else
    return log(r);
}

void InverseRGreensFn::Init(const KeyValueFile* kvf) throw (Exception)
{
  double d;
  const Matd* m;

  if (!kvf->GetMatd("X", m)) throw Exception("Missing X.");
  _x = *m;
  if (_x.Size(1) != 3) throw Exception("X must be 3xN.");

  if (kvf->GetDouble("order", d)) _order = (uint) d;

  _delta = 0.0;
  kvf->GetDouble("delta", _delta);
  if (_delta < 0.0) throw Exception("delta must be >= 0.");
}

bool InverseRGreensFn::
Call(const vector<uint>& rs, const vector<uint>& cs, double* B)
{
  for (uint k = 0, ic = 0; ic < cs.size(); ic++)
    for (uint ir = 0; ir < rs.size(); ir++, k++)
      B[k] = Eval(rs[ir], cs[ic]);
  return true;
}
