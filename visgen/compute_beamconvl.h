#ifndef COMPUTE_BEAMCONV_H
#define COMPUTE_BEAMCONV_H

int
compute_beamconvl (/* beams, st1, st2, ionoscreen, freq, do_ion, beamconvl) */
struct stnbeam *stbeams,
int st1_index,
int st2_index,
struct iono_screen *ionoscreen,
double freq,
int do_ion,
struct conv_fn *beamconvl
);

void compute_beamconvl_set_debug(int level);

#endif

