namespace hmmvp {
  namespace es {

    inline double GetTractionComp
    (const LameParms& lp, const dc3::Elem& es, const double* disl,
     const dc3::Elem& eo, size_t component)
    {
      double u[3], du[9], s[6];
      Dc3d(lp, es, disl, eo.Center(), u, du);
      DuToS(lp, du, s);
      double B;
      switch (component) {
      case 0:
        ProjectStress(s, eo.Normal(), eo.AlongStrike(), NULL, &B, NULL, NULL);
        break;
      case 1:
        ProjectStress(s, eo.Normal(), eo.AlongDip(),    NULL, &B, NULL, NULL);
        break;
      case 2:
        ProjectStress(s, eo.Normal(), NULL,             NULL, &B, NULL, &B  );
        break;
      }
      return B;
    }
    
  }
}
