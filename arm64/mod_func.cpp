#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
#if defined(__cplusplus)
extern "C" {
#endif

extern void _CaPQ_reg(void);
extern void _Kv31_reg(void);
extern void _ampa_reg(void);
extern void _cad_reg(void);
extern void _cadiv_reg(void);
extern void _can2_reg(void);
extern void _car_reg(void);
extern void _gabaa_reg(void);
extern void _nmda_reg(void);

void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"CaPQ.mod\"");
    fprintf(stderr, " \"Kv31.mod\"");
    fprintf(stderr, " \"ampa.mod\"");
    fprintf(stderr, " \"cad.mod\"");
    fprintf(stderr, " \"cadiv.mod\"");
    fprintf(stderr, " \"can2.mod\"");
    fprintf(stderr, " \"car.mod\"");
    fprintf(stderr, " \"gabaa.mod\"");
    fprintf(stderr, " \"nmda.mod\"");
    fprintf(stderr, "\n");
  }
  _CaPQ_reg();
  _Kv31_reg();
  _ampa_reg();
  _cad_reg();
  _cadiv_reg();
  _can2_reg();
  _car_reg();
  _gabaa_reg();
  _nmda_reg();
}

#if defined(__cplusplus)
}
#endif
