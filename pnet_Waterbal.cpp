#include "pnet_model.h"

void pnet_model::Waterbal(site_struct* site, veg_struct* veg, clim_struct* clim, int rstep, share_struct* share)
{
	//
	// PnET ecosystem water balance routine.
	//


	int wday;
	double	CanopyGrossPsnMG, WUE, PotTransd, Rain, Snow, Evap, precrem, SnowFrac, Tavew, SnowMelt, waterin, waterind;
	double	Trans, TotSoilMoistEff, Transd, CanopyGrossPsnAct;
	double	dTemp;

	//Potential transpiration
	CanopyGrossPsnMG = share->CanopyGrossPsn * 1000.0 * (44.0 / 12.0);  // mg CO2/m2.day
	WUE = (veg->WUEConst / share->VPD) * share->DWUE;  // g CO2/kg water
	//DelAmax added to the following eqn to prevent increased transpiration at high CO2
//	PotTransd = (CanopyGrossPsnMG / share->DelAmax / WUE) / 10000.0;   //cm water/day
	PotTransd = (CanopyGrossPsnMG /share->DVPD/ share->DelAmax / WUE) / 10000.0;

	// Ts-only phase parameters (from your calibrated/selected Ts-only phase model)
	const double alpha = 0.6273;
	const double beta = -1.24;

	//snow fraction 
	// linear predictor
	double eta = alpha + beta * share->Tw;
	if (eta >= 50.0)  SnowFrac = 1.0;
	else if (eta <= -50.0) SnowFrac = 0.0;
	else SnowFrac = 1.0 / (1.0 + std::exp(-eta));

	Rain = clim->prec[rstep] * (1 - SnowFrac);
	Snow = clim->prec[rstep] * SnowFrac;
	Evap = Rain * veg->PrecIntFrac; // all intercepted rain evaporated
	precrem = Rain - Evap;


	//precip modify with urbanization
	if (share->In > 0) {
		Evap = Evap * (1 - share->ISA);
		Snow = Snow * (1 - share->ISA);
		precrem = precrem * (1 + share->ISA);
	}

	site->SnowPack = site->SnowPack + Snow;

	//========================================
	//Snowmelt
	if (site->SnowPack > 0)
	{

		if (clim->timestep == 0)  //monthly version
		{
			Tavew = share->Tave;
			if (share->Tave < 1.0)
			{
				Tavew = 1.0;
			}

			SnowMelt = 0.15 * Tavew * share->dayspan;
		}

		if (clim->timestep == 1)  //daily version
		{
			// --- initialize albedo if needed ---
			if (!std::isfinite(site->Albedo)) site->Albedo = 0.85;
			if (Snow >= 0.2) site->Albedo = 0.85;

			// --- dynamic albedo update ---
			double tau = (share->Tave > 0.0 ? 7.0 : 25.0);
			tau = tau / (1.0 + 3.0 * std::max(Rain, 0.0));
			if (tau < 1.0) tau = 1.0;
			site->Albedo = 0.55 + (site->Albedo - 0.55) * std::exp(-1.0 / tau);

			// ---- CC parameters (match your R runner / final knobs) ----
			const double b_temp = 0.200;   // cm/°C/day
			const double b_rad = 0.002;   // cm/MJ/m²
			const double k_cc = 0.001;   // cold content rate
			const double cc_cap = 5.0;     // cold content upper limit

			// --- CC safety guards ---
			if (!std::isfinite(share->ColdContent) ||
				share->ColdContent < 0.0)
				share->ColdContent = 0.0;
			if (site->SnowPack <= 1e-12)
				share->ColdContent = 0.0;

			// --- CC accumulation ---
			// matches R: cc <- cc + k_cc * (-Tave)  if Ice>0 && Tave<0
			if (site->SnowPack > 1e-12 && share->Tave < 0.0) {
				share->ColdContent += (-share->Tave) * k_cc;
				if (share->ColdContent > cc_cap)
					share->ColdContent = cc_cap;
			}

			// --- melt potential (standard degree-day + radiation) ---
			// matches R: melt_pot = b_temp*max(Tday,0) + b_rad*R_net
			double Tday_eff = (share->Tday > 0.0 ? share->Tday : 0.0);
			double sm_pot = b_temp * Tday_eff;

			double Rs = clim->par[rstep] * 0.04235;
			double R_net = std::max(Rs * (1.0 - site->Albedo), 0.0);
			sm_pot += b_rad * R_net;

			if (sm_pot < 0.0) sm_pot = 0.0;

			// --- CC whole-pack gating ---
			// matches R final logic (not partial r_pay)
			double sm_actual;
			if (sm_pot > 0.0 && share->ColdContent > 0.0) {
				if (sm_pot >= share->ColdContent) {
					sm_actual = sm_pot - share->ColdContent;
					share->ColdContent = 0.0;
				}
				else {
					sm_actual = 0.0;
					share->ColdContent -= sm_pot;
				}
			}
			else {
				sm_actual = sm_pot;
			}
			if (sm_actual < 0.0) sm_actual = 0.0;

			// ---------- total melt production ----------
			SnowMelt = sm_actual;

			if (SnowMelt < 0.0) SnowMelt = 0.0;
			if (SnowMelt > site->SnowPack) SnowMelt = site->SnowPack;

		}

	}
	else
	{
		SnowMelt = 0;
	}


	// (1) solid snow decreases
	double snowpack_ice_before = site->SnowPack;
	site->SnowPack -= SnowMelt;
	if (site->SnowPack < 0.0) site->SnowPack = 0.0;
	// rain goes directly to soil
	double rain_to_soil = precrem;   // default: no snow -> goes to soil
	

	// liquid input this timestep (cm/period)
	waterin = SnowMelt + rain_to_soil;          // keep name: waterin = total liquid input

	// daily-mean diagnostic (cm/day), keep name waterind
	waterind =  waterin / share->dayspan;

	//Transpiration
	Trans = 0;

	TotSoilMoistEff = 0;
	for (wday = 1; wday <= share->dayspan; wday++)
	{
			share->Water = share->Water + waterind;

			if (PotTransd > 0.0) {
				if (share->Water >= PotTransd / veg->f)
				{
					Transd = PotTransd;
				}
				else
				{
					Transd = share->Water * veg->f;
				}
				share->Water = share->Water - Transd;
				Trans += Transd;
			}
			else
			{
				Transd = 0;
			}

			if (share->Water < site->WHC)
			{
				dTemp = share->Water;
			}
			else
			{
				dTemp = site->WHC;
			}
			TotSoilMoistEff = TotSoilMoistEff + pow((dTemp / site->WHC), (1.0 + veg->SoilMoistFact));

			const double Dmax_BC = 0.022;  // cm/day
			const double perc_thr_frac = 0.92;   // 92% WHC
			const double freeze_SWE_thr = 6.0;   // cm
			const double freeze_T_thr = 1.0;    // °C

			bool soil_unfrozen = (site->SnowPack >= freeze_SWE_thr) || (share->Tave >= freeze_T_thr);

			if (soil_unfrozen && share->Water > perc_thr_frac * site->WHC) {
				double avail = share->Water - perc_thr_frac * site->WHC;
				double range = site->WHC - perc_thr_frac * site->WHC;
				double perc_BC = Dmax_BC * (avail / range);  // b_BC=1 
				if (perc_BC > avail) perc_BC = avail;
				if (perc_BC < 0.0)   perc_BC = 0.0;
				share->Water -= perc_BC;
			}

	}// end of daily loop


		share->MeanSoilMoistEff = TotSoilMoistEff / share->dayspan;  // bug correct, zzx
		if (share->MeanSoilMoistEff > 1)share->MeanSoilMoistEff = 1.0;


		//Water stress
		if (PotTransd > 0.0 && share->dayspan > 0)
			share->DWater = Trans / (PotTransd * share->dayspan);
		else
			share->DWater = 1.0;

		/////////////Manually turn water stress off  //
		if (site->WaterStress == 0) share->DWater = 1;

		if (share->DWater > 1.0) share->DWater = 1.0;
		//	share->DWater = 1.0;

		share->Dwatertot = share->Dwatertot + (share->DWater * share->dayspan);
		share->DwaterIx = share->DwaterIx + share->dayspan;



	//Calculate actural ozone effect and NetPsn with drought stress

	if ((share->O3EffectOnPSN == 1) && (share->CanopyGrossPsn > 0))
	{
		if (share->WUEO3Eff == 0)
		{
			//no O3 effect on WUE (assumes no stomatal impairment)
			share->CanopyDO3 = share->CanopyDO3Pot + ((1 - share->CanopyDO3Pot) * (1 - share->DWater));  // drought offset to O3 effect
		}
		else
		{
			//reduce the degree to which drought offsets O3 (assumes stomatal impairment in proportion to effects on psn)
			share->CanopyDO3 = share->CanopyDO3Pot + (1 - share->CanopyDO3Pot) * (1 - share->DWater / share->CanopyDO3Pot);

		}
		share->DroughtO3Frac = share->CanopyDO3Pot / share->CanopyDO3;  // relative drought offset
	}
	else
	{
		share->CanopyDO3 = 1;
		share->DroughtO3Frac = 1;
	}

	if (share->CanopyDO3 > 1)  share->CanopyDO3 = 1;

	if (share->CanopyGrossPsn > 0)
	{
		share->CanopyDO3Tot = share->CanopyDO3Tot + share->CanopyDO3 * share->dayspan;
		share->CanopyDO3TotDay = share->CanopyDO3TotDay + share->dayspan;

	}


	CanopyGrossPsnAct = share->CanopyGrossPsn * share->DWater;   // ' g C/m2.day
	share->CanopyGrossPsnActMo = CanopyGrossPsnAct * share->dayspan;  // g C/m2.month
	share->GrsPsnMo = share->CanopyGrossPsnActMo * share->CanopyDO3;
	share->NetPsnMo = share->CanopyNetPsn * share->DWater * share->dayspan * share->CanopyDO3; //Linghui 05182020

	share->TotCanopyGrossPsnMax = share->TotCanopyGrossPsnMax + share->CanopyGrossPsn * share->dayspan;  // psn without water stress

	/*
	if (share->Water > site->WHC)
	{
		share->Drainage = share->Water - site->WHC;
		share->Water = site->WHC;
	}
	else
	{
		share->Drainage = 0;
	}*/

	// runoff/outflow determined AFTER ET by WHC
	double Outflow;
	if (share->Water > site->WHC)
	{
		Outflow = share->Water - site->WHC;   // cm/period
		share->Water = site->WHC;
	}
	else
	{
		Outflow = 0.0;
	}

	share->Drainage = Outflow;

	//======================================================
	// NOW allocate snowmelt N AFTER runoff is known
	//======================================================
	double ratio = (snowpack_ice_before > 1e-12) ? SnowMelt / snowpack_ice_before : 0.0;
	ratio = std::min(1.0, std::max(0.0, ratio));

	// snow N lost with outflow vs retained in soil
	double relNO3 = share->SnowNO3 * ratio;
	double relNH4 = share->SnowNH4 * ratio;
	share->NO3 += relNO3;
	share->NH4 += relNH4;
	share->SnowNO3 -= relNO3;
	share->SnowNH4 -= relNH4;

	share->FracDrain = share->Drainage / (share->Water + share->Drainage);
	share->TotTrans = share->TotTrans + Trans;
	share->TotWater = share->TotWater + share->Water;
	share->TotPsn = share->TotPsn + share->NetPsnMo;
	share->TotDrain = share->TotDrain + share->Drainage;
	share->TotPrec = share->TotPrec + clim->prec[rstep];
	share->TotEvap = share->TotEvap + Evap;
	share->TotGrossPsn = share->TotGrossPsn + share->GrsPsnMo;
	share->ET = Trans + Evap;

}
