#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
#if defined(__cplusplus)
extern "C" {
#endif

extern void _cacumm_reg(void);
extern void _cagk_reg(void);
extern void _cal2_reg(void);
extern void _can2_reg(void);
extern void _capmp_reg(void);
extern void _cat_reg(void);
extern void _distr_reg(void);
extern void _exp2nmdar_reg(void);
extern void _ExpAmpar_reg(void);
extern void _fc_reg(void);
extern void _h_reg(void);
extern void _KahpM95_reg(void);
extern void _kaprox_reg(void);
extern void _kcasimple_reg(void);
extern void _kd_reg(void);
extern void _kdrca1_reg(void);
extern void _kir_reg(void);
extern void _km_reg(void);
extern void _leak_chan_reg(void);
extern void _na3n_reg(void);
extern void _naxn_reg(void);
extern void _syn1_reg(void);

void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"cacumm.mod\"");
    fprintf(stderr, " \"cagk.mod\"");
    fprintf(stderr, " \"cal2.mod\"");
    fprintf(stderr, " \"can2.mod\"");
    fprintf(stderr, " \"capmp.mod\"");
    fprintf(stderr, " \"cat.mod\"");
    fprintf(stderr, " \"distr.mod\"");
    fprintf(stderr, " \"exp2nmdar.mod\"");
    fprintf(stderr, " \"ExpAmpar.mod\"");
    fprintf(stderr, " \"fc.mod\"");
    fprintf(stderr, " \"h.mod\"");
    fprintf(stderr, " \"KahpM95.mod\"");
    fprintf(stderr, " \"kaprox.mod\"");
    fprintf(stderr, " \"kcasimple.mod\"");
    fprintf(stderr, " \"kd.mod\"");
    fprintf(stderr, " \"kdrca1.mod\"");
    fprintf(stderr, " \"kir.mod\"");
    fprintf(stderr, " \"km.mod\"");
    fprintf(stderr, " \"leak_chan.mod\"");
    fprintf(stderr, " \"na3n.mod\"");
    fprintf(stderr, " \"naxn.mod\"");
    fprintf(stderr, " \"syn1.mod\"");
    fprintf(stderr, "\n");
  }
  _cacumm_reg();
  _cagk_reg();
  _cal2_reg();
  _can2_reg();
  _capmp_reg();
  _cat_reg();
  _distr_reg();
  _exp2nmdar_reg();
  _ExpAmpar_reg();
  _fc_reg();
  _h_reg();
  _KahpM95_reg();
  _kaprox_reg();
  _kcasimple_reg();
  _kd_reg();
  _kdrca1_reg();
  _kir_reg();
  _km_reg();
  _leak_chan_reg();
  _na3n_reg();
  _naxn_reg();
  _syn1_reg();
}

#if defined(__cplusplus)
}
#endif
