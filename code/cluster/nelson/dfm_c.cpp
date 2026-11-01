// C interface to the Nelson (2000) dead fuel moisture model as implemented in NFDRS4 (firelab/NFDRS4, public domain).
// run_dfm: hourly series -> stick mean weighted moisture (g/g) after each hour, for a stick of the given time lag class.
#include "deadfuelmoisture.h"
extern "C" int run_dfm(int lag, int n, const float* T, const float* RH, const float* SW, const float* R, float* out) {
    DeadFuelMoisture* d = lag == 1 ? DeadFuelMoisture::createDeadFuelMoisture1("s") : DeadFuelMoisture::createDeadFuelMoisture10("s");
    d->initializeEnvironment(20., 0.4, 0., 0., 20., 0.4, 0.2, 0.0218);
    int bad = 0;
    for (int i = 0; i < n; i++) {
        double t = T[i], h = RH[i], s = SW[i], r = R[i];
        if (!(t == t) || !(h == h)) { out[i] = out[i > 0 ? i - 1 : 0]; bad++; continue; }
        if (h < 0.01) h = 0.01; if (h > 1.0) h = 1.0; if (!(s == s) || s < 0) s = 0; if (s > 1400) s = 1400; if (!(r == r) || r < 0) r = 0;
        d->update(1.0, t, h, s, r, 0.0218, true);
        out[i] = (float)d->meanWtdMoisture();
    }
    delete d; return bad;
}
