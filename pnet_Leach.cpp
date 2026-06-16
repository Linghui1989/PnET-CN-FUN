#include "pnet_model.h"

void pnet_model::Leach(clim_struct* clim, int rstep, share_struct* share)
{
	//
	//PnET-CN leaching routine
	//

	// ---- seasonal f_retain lookup ----
	double f_retain;
	int doy = clim->doy[rstep];
	if (doy <= 90 || doy >= 335) f_retain = 0.60;  // Winter
	else if (doy <= 181)                      f_retain = 0.33;  // Spring
	else if (doy <= 273)                      f_retain = 0.49;  // Summer
	else                                      f_retain = 0.14;  // Fall

	f_retain = std::min(0.96, std::max(0.10, f_retain));

	share->NDrain = share->FracDrain * share->NO3 * (1 - f_retain);
	share->NO3 -= share->NDrain;
	share->NDrainYr = share->NDrainYr + share->NDrain;
	share->NDrainMgL = (share->NDrain * 1000) / (share->Drainage * 10+1.0E-6);  //to convert to mg/l

}
