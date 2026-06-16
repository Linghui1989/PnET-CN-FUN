#include "pnet_model.h"

void pnet_model::Phenology(veg_struct* veg, clim_struct* clim, int rstep, share_struct* share, int GrowthPhase)
{
	//
	// Phenlology calculations for the PnET ecosystem model.
	//

	double OldFolMass, GDDFolEff, delGDDFolEff, FolMassNew, GDD, CDD;
	double t0, t0_chill, T_base, T_opt, T_min, T_max, F_crit1, F_crit2, C_req, Rc, Tc;
	int k;
	double spring_mod_select, fall_mod_select;
	double fall_t0, fall_T_base, fall_F_crit1, fall_F_crit2;

	if (GrowthPhase == 1)
	{
		switch (clim->timestep)
		{
		case 0:   // monthly
			share->dayspan = (double)getdays(clim->year[rstep], clim->doy[rstep]);
			break;
		case 1:   // daily
			share->dayspan = 1;
			break;
		default:
			share->dayspan = (double)getdays(clim->year[rstep], clim->doy[rstep]);  //default for monthly
			break;
		}


		GDD = share->Tave * share->dayspan;

		share->TaveYr += GDD / 365.0;
		share->PARYr += clim->par[rstep] * share->dayspan / 365.0;

		if (GDD < 0 || clim->doy[rstep] < 60)  GDD = 0;   // need modification for tropical region

		if (share->GDDTot < veg->GDDFolStart && share->GDDTot + GDD >= veg->GDDFolStart) { share->FolGDD = 1; }

		share->GDDTot += GDD;

		//add soil freeze algorithm
		if (clim->doy[rstep] < veg->GDDFolStart || clim->doy[rstep] > veg->SenescStart) {

			// Update freezing counter
			share->freezeDoy += (share->Tave <= 0);
			share->freezeDoy -= (share->Tave > 0);

			// Bound freezeDoy between 0 and 5
			if (share->freezeDoy < 0) share->freezeDoy = 0;
			if (share->freezeDoy > 5) share->freezeDoy = 5;

			// Binary freeze flag
			share->freeze = (share->freezeDoy == 5); //If the freezing-day counter equals 5, set freeze to 1; otherwise set it to 0.

		}
		else {
			// Growing season: no freezing
			share->freezeDoy = 0;
			share->freeze = 0;
		}


		
		spring_mod_select = 0;
		

		if (spring_mod_select == 0) {

			F_crit1 = veg->GDDFolStart; //GDDFolStart 332
			F_crit2 = veg->GDDFolEnd; //GDDFolEnd

			if (clim->doy[rstep] < 60)  Tc = 0, share->Tsum = 0;
			Tc = share->Tave;

		}
		else if (spring_mod_select == 1)
		{
			// SEQUENTIAL M1 (SM1) MODEL;
			// Temperature forcing, photoperiod, and chilling requirements
			// parameters for SM1 model

			t0 = 92;//par[1]
			t0_chill = 74;//par[2]
			T_base = 4;//par[3]
			T_opt = 9.96;//par[4]
			T_min = -1.34;//par[5]
			T_max = 14.86;//par[6]
			F_crit1 = veg->GDDFolStart;//par[7]
			F_crit2 = veg->GDDFolEnd;//par[7.5]
			C_req = 4.36;//par[8]

			//chilling triangular temperature response
			if (share->Tave <= T_min) Rc = 0;
			else if (share->Tave > T_min && share->Tave < T_opt) Rc = (share->Tave - T_min) / (T_opt - T_min);
			else if (share->Tave > T_opt && share->Tave < T_max) Rc = (share->Tave - T_max) / (T_opt - T_max);
			else if (share->Tave >= T_max) Rc = 0;

			if (clim->doy[rstep] <= t0_chill)  Rc = 0;
			share->Sc = Rc;

			if (share->Sc > C_req) k = 0;
			else int k = 1;

			//forcing
			Tc = share->Tave - T_base;
			if (Tc < 0)  Tc = 0;
			Tc = pow(((share->DayLength / 3600) / 24), k) * Tc;
			if (clim->doy[rstep] <= t0_chill)  Tc = 0, share->Tsum = 0, share->CDDTot = 0;
		}

		// BUDBURST
		share->Tsum += Tc * share->dayspan;

		if ((share->Tsum >= F_crit1) && (clim->doy[rstep] < veg->SenescStart))
		{
			OldFolMass = share->FolMass;
			GDDFolEff = (share->Tsum - F_crit1) / (F_crit2 - F_crit1);

			if (GDDFolEff < 0)GDDFolEff = 0;
			if (GDDFolEff > 1)GDDFolEff = 1;

			delGDDFolEff = GDDFolEff - share->OldGDDFolEff;
			share->FolMass = share->FolMass + (share->BudC * delGDDFolEff) / veg->CFracBiomass;
			share->FolProdCMo = (share->FolMass - OldFolMass) * veg->CFracBiomass;
			share->FolGRespMo = share->FolProdCMo * veg->GRespFrac;
			share->OldGDDFolEff = GDDFolEff;
		}
		else
		{
			share->FolProdCMo = 0;
			share->FolGRespMo = 0;
		}
	}
	else
	{
		share->FolLitM = 0;

		fall_mod_select = 1;

		if (fall_mod_select == 0) {

			if ((share->PosCBalMass < share->FolMass) && (clim->doy[rstep] > veg->SenescStart))
			{


				if ((share->PosCBalMass) > (veg->FolMassMin))
				{
					FolMassNew = share->PosCBalMass;
				}
				else
				{
					FolMassNew = veg->FolMassMin;
				}


				if (FolMassNew == 0)
				{
					share->LAI = 0;
				}
				else if (FolMassNew < share->FolMass)
				{
					share->LAI = share->LAI * (FolMassNew / share->FolMass);
				}
				if (FolMassNew < share->FolMass)
				{
					share->FolLitM = share->FolMass - FolMassNew;
				}
				share->FolMass = FolMassNew;
			}
		}
		else if (fall_mod_select == 1)
		{
			

			fall_t0 = veg->SenescStart; //226
			fall_T_base = 25.08; //25.08; //
			fall_F_crit1 = -384.80; // threshold for start of leaf senescence
			fall_F_crit2 = -784.22; // threshold for end of leaf senescence

			CDD = (share->Tave - fall_T_base) * share->dayspan;

			if (CDD > 0) CDD = 0;

			if (clim->doy[rstep] <= fall_t0)  CDD = 0;

			share->CDDTot += CDD;

			share->FolLitM = 0;

			if ((share->CDDTot <= fall_F_crit1) && (clim->doy[rstep] > fall_t0))
			{
				if (share->CDDTot > fall_F_crit2)
				{
					FolMassNew = share->FolMass * (fall_F_crit2 - share->CDDTot) / (fall_F_crit2 - fall_F_crit1);
				}
				else if (share->CDDTot < fall_F_crit2)
				{
					FolMassNew = veg->FolMassMin;
				}
				if (FolMassNew == 0)
				{
					share->LAI = 0;
				}
				else if (FolMassNew < share->FolMass)
				{
					share->LAI = share->LAI * (FolMassNew / share->FolMass);
				}
				if (FolMassNew < share->FolMass)
				{
					share->FolLitM = share->FolMass - FolMassNew;
				}
				share->FolMass = FolMassNew;
			}
		}

		
	}


}