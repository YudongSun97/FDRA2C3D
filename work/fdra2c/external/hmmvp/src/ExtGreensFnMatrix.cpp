/* hmmvp: Software to form and apply Hierarchical Matrices
 *   Version 1.3
 *   Andrew M. Bradley
 *   ambrad@cs.stanford.edu
 *   CDFM Group, Geophysics, Stanford
 *   https://pangea.stanford.edu/research/CDFM/software
 * hmmvp is licensed as follows:
 *   Open Source Initiative OSI - Eclipse Public License 1.0
 *   http://www.opensource.org/licenses/eclipse-1.
 *
 * This function is written by Camilla Cattania. It reads precalculated Green's function matrix.
*/

class ExtGreensFn : public ImplGreensFn {
public:
  virtual void Init(const KeyValueFile* kvf) throw (Exception);
  virtual Hd* ComputeHd (double eta) { return NewHd(_xs, _xr, NULL, eta); }
  virtual bool Call(const CompressBlockInfo& cbi, const vector<UInt>& rs,
                    const vector<UInt>& cs, double* B) const;

private:
  Matd _G;
  Matd _xs, _xr;
  UInt _ns, _nr;
};

void ExtGreensFn::Init (const KeyValueFile* kvf) throw (Exception) {
  double d;
  const Matd* m;

  if (!kvf->GetMatd("xs", m)) throw Exception("Missing xs.");
  _xs = *m;

  if (!kvf->GetMatd("xr", m)) throw Exception("Missing xr.");
  _xr = *m;

  if (!kvf->GetMatd("G", m)) throw Exception("Missing G.");
  _G = *m;

  if (kvf->GetDouble("Ns", d)) _ns = (UInt) d;
  if (kvf->GetDouble("Nr", d)) _nr = (UInt) d;

  if (_G.Size(1) != _nr) throw Exception("X must be Nr x Ns.");
  if (_G.Size(2) != _ns) throw Exception("X must be Nr x Ns.");

}

bool ExtGreensFn::
Call (const CompressBlockInfo& cbi, const vector<UInt>& rs,
      const vector<UInt>& cs, double* B) const {
  for (UInt k = 0, ic = 0; ic < cs.size(); ic++)
    for (UInt ir = 0; ir < rs.size(); ir++, k++)
      B[k] = _G(rs[ir], cs[ic]);
  return true;
}
