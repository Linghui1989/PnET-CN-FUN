#include "pnet_model.h"

void pnet_model::Decomp(site_struct* site, veg_struct* veg, clim_struct* clim, int rstep, share_struct* share)
{
	//
	//PnET-CN decomposition routine
	//

	double tEffSoil, TMult, WMult, KhoAct, DHO, GrossNMin, SoilPctN, NReten, GrossNImmob, NetNMin, NetNitr, RootNSinkStr, PlantNUptake, NH4Up, NO3Up;
	double NH4, NO3;
	double FN2O, FNO, FN2, FNOd, FN2Od, NitLostGas, Kn, Fdeni;  // N gas emission
	double fNRatioNit = 1.0; 		// factor to adjust daily net Nitrification.
	double R_NOx_N2O;  // ratio of NO to N2O through nitrification
	double Wnit; 		// water stress
	double den_T, Q10, den_W, Wcrit, R_N2_N2O, den_CO2, Ccrit, Cmax;
	int iday;


	FN2O = 0;
	FNO = 0;
	FN2 = 0;
	FN2Od = 0;
	Kn = 0.01; //Kn=0.02, Kn is the fraction of nitrified N lost as N2O/NOx  (or0.05) (Davidson,2000;parton,2001)



	//Temperature effect on all soil processes
	for (iday = 1; iday <= share->dayspan; iday++)
	{
		//Add atmospheric N deposition
		if (share->SnowStatus == 1) {

			share->SnowNO3 += clim->NO3dep[rstep] / share->dayspan;
			share->SnowNH4 += clim->NH4dep[rstep] / share->dayspan;
		}
		else
		{
			share->NO3 += clim->NO3dep[rstep] / share->dayspan;
			share->NH4 += clim->NH4dep[rstep] / share->dayspan;

		}

		share->NdepTot = share->NdepTot + clim->NO3dep[rstep] / share->dayspan + clim->NH4dep[rstep] / share->dayspan;  //ZZX


		if (share->Tave > 1)
		{
			tEffSoil = share->Tave;
		}
		else
		{
			tEffSoil = 1;
		}

		TMult = (exp(0.1 * (share->Tave - 7.1)) * 0.68) * 1;

		WMult = share->MeanSoilMoistEff;
		if (site->WaterStress == 0) WMult = 1;

				//Add litter to Humus pool
		share->HOM = share->HOM + share->TotalLitterM / share->dayspan;
		share->HON = share->HON + share->TotalLitterN / share->dayspan;

		//Humus dynamics

		KhoAct = veg->Kho / 365.0; //* share->dayspan;
		DHO = share->HOM * (1 - exp(-KhoAct * TMult * WMult));
		GrossNMin = DHO * (share->HON / share->HOM);

		const double a_ndep = -1.814;
		const double b_ndep = 17.401;
		const double Ndep_ref = 0.517;  // 2003-2005 anchor

		double f_corr = 1.0;  // default: no correction

		if (clim->year[rstep] > 2000) {
			double Ndep_corr = (share->NdepTot_prev > 0.0)
				? share->NdepTot_prev
				: Ndep_ref;
			f_corr = (a_ndep + b_ndep * Ndep_corr) /
				(a_ndep + b_ndep * Ndep_ref);
			f_corr = std::max(0.05, std::min(2.0, f_corr));
		}

		share->EnvMaxFol = f_corr;

		GrossNMin *= f_corr;

		share->SoilDecResp = DHO * veg->CFracBiomass;
		share->SoilDecRespYr = share->SoilDecRespYr + share->SoilDecResp;
		share->GrossNMinYr = share->GrossNMinYr + GrossNMin;
		share->HON -= GrossNMin;
		share->HOM -= DHO;

		share->NetCBal = share->NetCBal - share->SoilDecResp; // updating NetCBal

		// Immobilization and net mineralization
		if (modeltype == 4)
		{
			// PnET-FUN-CORPSE uses DecompFUN_CORPSE() for belowground N cycling.
			// If this legacy routine is reached, do not apply the fixed NReten
			// immobilization formulation on top of the FUN-CORPSE bookkeeping.
			GrossNImmob = 0.0;
			NetNMin = GrossNMin;
		}
		else
		{
			SoilPctN = (share->HON / share->HOM) * 100;
			NReten = (veg->NImmobA + veg->NImmobB * SoilPctN) / 100;
			if (NReten > 1) { NReten = 1; }
			if (NReten < 0) { NReten = 0; }

			GrossNImmob = NReten * GrossNMin;

			double MicrobialCN = 8.0;
			double CImmob = GrossNImmob * MicrobialCN;

			share->HOM += CImmob;
			share->HON += GrossNImmob;
			NetNMin = GrossNMin - GrossNImmob;
		}

		share->GrossNImmobYr = share->GrossNImmobYr + GrossNImmob;
		share->NH4 = share->NH4 + NetNMin;


		// nitrification========================================================
		double nit_season_factor;
		int doy = clim->doy[rstep];

		if (doy <= 90 || doy >= 335)
			nit_season_factor = 0.462;        // winter
		else if (doy <= 181)
			nit_season_factor = 0.856;      // Spring
		else if (doy <= 273)
			nit_season_factor = 1.000;      // Summer
		else
			nit_season_factor = 0.819;      // Fall


		fNRatioNit = share->dayspan/30.4;  // daily step; Linghui Meng 20201221
		NetNitr = share->NH4 * share->NRatioNit * fNRatioNit * nit_season_factor;

		share->WFPS = share->Water / site->WHC * 0.75;
		if (share->WFPS > 1.0)share->WFPS = 1.0;

		Wnit = 1.0;
		NitLostGas = Kn * NetNitr * Wnit;
		NitLostGas = 0;

		R_NOx_N2O = pow(10.0, -3.79 * share->WFPS + 2.73);
		FN2O = NitLostGas / (1.0 + R_NOx_N2O);
		FNO = NitLostGas - FN2O;

		share->NO3 = share->NO3 + NetNitr;
		share->NH4 = share->NH4 - NetNitr - NitLostGas;


		// soil denitrification
		Q10 = 2.0;
		share->Tsoil = share->Tave;
		if (0.0 < share->Tsoil)  den_T = pow(Q10, 0.1 * (share->Tsoil - 20.0));
		else    den_T = 0.01;

		Wcrit = 0.75;
		den_W = 1.0 - 1.0 / (1.0 + exp(12.5 * (share->WFPS - Wcrit)));

		Ccrit = 1.8;
		Cmax = 3;
		den_CO2 = 1.0 - 1.0 / (1.0 + exp(6.8 * (share->SoilDecResp - Ccrit) / Cmax));
		den_CO2 = 1.0;

		site->Kden = 0.03;

		if (den_W > den_CO2) den_W = den_CO2;
		Fdeni = share->NO3 * site->Kden * den_T * den_W;

		R_N2_N2O = 1.44 * exp(1.33 * share->WFPS);

		double RnA, RnB;

		share->RnMax = 10.0;
		share->RnX1 = 0.7;
		share->RnY1 = 1.5;

		share->RnX2 = 1.0;
		share->RnY2 = share->RnMax;

		RnB = log(share->RnY2 / share->RnY1) / (share->RnX2 - share->RnX1);
		RnA = share->RnY1 * exp(-RnB * share->RnX1);

		R_N2_N2O = RnA * exp(RnB * share->WFPS);

		FN2Od = Fdeni * 1 / (1 + R_NOx_N2O + R_N2_N2O);
		FN2Od = 0.0;

		FNOd = R_NOx_N2O * FN2Od;
		FN2 = FN2Od * R_N2_N2O;

		share->FluxN2OYrDe += FN2Od;
		share->FluxNOYrDe += FNOd;

		share->FluxN2OYr += FN2O + FN2Od;
		share->FluxNOYr += FNO + FNOd;
		share->FluxN2Yr += FN2;

		share->NO3 -= FN2Od + FN2 + FNOd;

		//Plant Uptake

		RootNSinkStr = share->RootNSinkEff * TMult;
		RootNSinkStr = RootNSinkStr*(share->dayspan / 30.4)*1.5;

		if (site->SnowPack > 0) {
			RootNSinkStr = 0;
		}
		else if (share->GDDTot>50 && share->FolMass == 0)
		{
			RootNSinkStr *= 0.3;
		}

		if (share->RootMass < 0.9 * share->OldRoot) {

			share->UptakeEff = share->RootMass / share->OldRoot;
			
		}
		else
		{
			share->UptakeEff = 1;
			share->OldRoot = 0;
		}


		PlantNUptake = (share->NH4 + share->NO3) * RootNSinkStr *share->UptakeEff;
		if ((PlantNUptake + share->PlantN) > veg->MaxNStore)
		{
			PlantNUptake = veg->MaxNStore - share->PlantN;
			RootNSinkStr = PlantNUptake / (share->NO3 + share->NH4);
		}
		if (PlantNUptake < 0)
		{
			PlantNUptake = 0;
			RootNSinkStr = 0;
		}
		share->PlantN +=  PlantNUptake;
		share->PlantNUptakeYr += PlantNUptake;

		NH4Up = share->NH4 * RootNSinkStr;
		share->NH4 = share->NH4 - NH4Up;
		NO3Up = share->NO3 * RootNSinkStr;
		share->NO3 = share->NO3 - NO3Up;


		share->NetNMinYr = share->NetNMinYr + NetNMin;
		share->NetNitrYr = share->NetNitrYr + NetNitr;

	}
}
