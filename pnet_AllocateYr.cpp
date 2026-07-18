#include "pnet_model.h"

void pnet_model::AllocateYr(site_struct* site, veg_struct* veg, clim_struct* clim, int rstep, share_struct* share, int CN_Mode)
{
	//
	//Annual C allocation for the PnET ecosystem model.
	//

	int i;
	double EnvMaxFol, SppMaxFol, FolRegen, BiomLossFrac,RootLossFrac, TotalC, nr;
	double PotLightEff;  // potential light for growth
	double RemoveFrac;
	double BudC, BudN, BudNmin; //param use to calculate potental plant pool
	double folnconnew;
	double RootCloss, RootNloss; //C and N loss during root disturbance
	double FolMassMax;

	//Check for a disturbance year
	FolRegen = 100;
	BiomLossFrac = 0;
	RootLossFrac = 0;

	for (i = 1; i <= site->distyrs; i++)
	{
		if (clim->year[rstep] == site->distyear[i])
		{
			BiomLossFrac = site->distintensity[i];
			RemoveFrac = site->distremove[i];    // remove out of field
			RootLossFrac = site->distRootloss[i];
			if (modeltype != 4)
			{
				share->HOM *= (1 - site->distsoilloss[i]);
				share->HON *= (1 - site->distsoilloss[i]);
			}
			share->FolMass *= (1 - BiomLossFrac);
			if (veg->FolMassMax < site->folregen[i])veg->FolMassMax = site->folregen[i];
			break;
		}
	}

	if (RootLossFrac>0)
	{
		share->OldRoot = share->RootMass;
		RootCloss = share->RootMass* RootLossFrac;
		RootNloss = share->RootMassN * RootLossFrac;

		share->RootMass -= RootCloss;
		share->RootMassN -= RootNloss;

		share->DeadWoodM += RootCloss * 0.9;
		share->DeadWoodN += RootNloss * 0.9;

		share->RootMRespYr += RootCloss * 0.1;
		share->SnowNO3 += RootNloss * 0.1;
	}

	//Agriculture  Linghui Meng 05062020 move agriculture to year end calculation
	if (clim->year[rstep] >= site->agstart && clim->year[rstep] < site->agstop)
	{

		share->WoodMass = share->WoodMass * (1 - site->agrem * (1 - site->agrem));
		share->WoodMassN = share->WoodMassN * (1 - site->agrem * (1 - site->agrem));

		share->DeadWoodM = share->DeadWoodM * (1 - site->agrem * (1 - site->agrem)); //Linghui 03252020
		share->DeadWoodN = share->DeadWoodN * (1 - site->agrem * (1 - site->agrem)); //Linghui 03252020

		share->PlantC *= (1 - site->agrem * (1 - site->agrem));
		share->PlantN *= (1 - site->agrem * (1 - site->agrem));

		if (modeltype != 4)
		{
			share->HON += 1.5; //Linghui Meng 11262020 assume N addition
		}
	}
	
	//if the decidous tree foliar mass greater than 0 at the end of year, make it to 0 and relocate C and N
	//transfer foliar to litter
	if ((veg->FolReten == 1)&&(share->FolMass>0))
	{
		double FolNLoss, Retrans, FolLitN;
		share->FolLitM = share->FolMass;
		FolNLoss = share->FolLitM * (share->FolNCon / 100);
		Retrans = FolNLoss * veg->FolNRetrans;
		share->PlantN += Retrans;
		FolLitN = FolNLoss - Retrans;

		share->TotalLitterM += share->FolLitM;
		share->TotalLitterN += FolLitN;

		if (modeltype != 4)
		{
			share->HOM = share->HOM + share->TotalLitterM / share->dayspan;
			share->HON = share->HON + share->TotalLitterN / share->dayspan;
			// delete transferred litter from litter pool in legacy soil mode
			share->TotalLitterM -= share->TotalLitterM / share->dayspan;
			share->TotalLitterN -= share->TotalLitterN / share->dayspan;
		}

		share->FolMass = 0;
		//in case plantN greater than maximum N storage
		if (share->PlantN > veg->MaxNStore)
		{
			if (modeltype == 4)
			{
				share->FUNOverflowToNH4Yr += (share->PlantN - veg->MaxNStore);
			}
			share->NH4 += (share->PlantN - veg->MaxNStore);  // ZZX revised
			share->PlantN = veg->MaxNStore;
		}
	}

	share->NPPFolYr = share->FolProdCYr / veg->CFracBiomass;
	share->NPPWoodYr = share->WoodProdCYr / veg->CFracBiomass;
	share->NPPRootYr = share->RootProdCYr / veg->CFracBiomass;

//	share->PlantC -= share->WoodGRespYr + share->WoodMRespYr; Linghui 01212026

	share->PlantC += share->WoodC;
	share->WoodC = 0;
	


	if (share->DwaterIx > 0)
	{
		share->AvgDWater = share->Dwatertot / share->DwaterIx;
		//	share->AvgDWater = 1.0;
	}
	else
	{
		share->AvgDWater = 1.0;
	}


	if (share->PosCBalMassIx > 0)
	{
		share->avgPCBM = (share->PosCBalMassTot / share->PosCBalMassIx);
	}
	else
	{
		share->avgPCBM = share->FolMass;
	}


	share->CanopyDO3Avg = share->CanopyDO3Tot / share->CanopyDO3TotDay;


	if (share->LightEffCBalIx > 0)
	{
		PotLightEff = (share->LightEffCBalTot / share->LightEffCBalIx);
	}
	else
	{
		PotLightEff = 1.0;
	}

	share->LightEffMin = PotLightEff;

	EnvMaxFol = (share->AvgDWater * share->avgPCBM) * (1.0 + (veg->FolRelGrowMax * share->LightEffMin));
	SppMaxFol = share->avgPCBM * (1.0 + (veg->FolRelGrowMax * share->LightEffMin));

//	share->EnvMaxFol = EnvMaxFol;
//	share->SppMaxFol = SppMaxFol;

	if (EnvMaxFol < SppMaxFol)
	{
		FolMassMax = EnvMaxFol;
	}
	else
	{
		FolMassMax = SppMaxFol;
	}

	//add a limitation to foliar growth
	if (FolMassMax > veg->FolMassMax) {
		FolMassMax = veg->FolMassMax;
	}

//	veg->FolMassMin = (veg->FolMassMax - veg->FolMassMax * (1.0 / veg->FolReten));
	veg->FolMassMin = veg->FolMassMax * (1.0 - 1.0 / veg->FolReten);
	share->BudC = ((veg->FolMassMax - share->FolMass) * veg->CFracBiomass);

	share->BudC = fmin(share->BudC, share->PlantC); //if BudC greater than plantC, let BudC equal to plantC

	if (share->BudC < 0)
	{
		share->BudC = 0;
		if (veg->FolReten > 1) share->BudC = share->FolMass * 1 / (veg->FolReten - 1) * veg->CFracBiomass * 0.5;  // evergreen has half budc
	}

	//if PnET-CN model, need to calculate BudN and relocate BudC base on available N
	if (CN_Mode == 1)
	{

		share->NRatio = 1 + (share->PlantN / veg->MaxNStore) * veg->FolNConRange;

		if (share->NRatio < 1)
		{
			share->NRatio = 1;
		}

		if (share->NRatio > (1 + veg->FolNConRange))
		{
			share->NRatio = 1 + veg->FolNConRange;
		}

		if (share->NRatio == 1)
		{
			share->NRatioNit = 0;
		}
		else
		{

			nr = share->NRatio - 1 - (veg->FolNConRange / 3);
			if (nr < 0)nr = 0;

			share->NRatioNit = (nr / (0.6667 * veg->FolNConRange)) * (nr / (0.6667 * veg->FolNConRange));
			if (share->NRatioNit > 1) share->NRatioNit = 1;

		}

		if (veg->FolReten == 1) {
			double Nconc;

			Nconc = veg->FolNCon * share->NRatio;

			if (Nconc > veg->FolNCon + veg->FolNConRange) {

				Nconc = veg->FolNCon + veg->FolNConRange;
			}

			BudN = (share->BudC / veg->CFracBiomass) * Nconc / 100; //theritical N need 

			if (share->PlantN < BudN) {

				BudNmin = (share->BudC / veg->CFracBiomass) * (veg->FolNCon - veg->FolNConRange) / 100 * share->NRatio;

				if (share->PlantN < BudNmin) {

					share->BudC = share->BudC * (share->PlantN / BudNmin);
					BudN = share->PlantN;
				}
				else
				{
					BudN = share->PlantN;
				}

			}

			share->BudN = BudN;

			share->PlantN -= share->BudN;

			// save current foliar N for output
			share->FolNConOld = share->FolNCon;

			folnconnew = (share->FolMass * (share->FolNCon / 100) + share->BudN) / (share->FolMass + (share->BudC / veg->CFracBiomass)) * 100;

			share->FolNCon = folnconnew;
		}

	}

	if (share->BudC > share->PlantC) {
		share->BudC = share->PlantC;
	}

	share->PlantC -= share->BudC;
	share->WoodC = (1 - veg->PlantCReserveFrac) * share->PlantC;
//	share->WoodC = 0.9 * share->BudC;
	share->PlantC -= share->WoodC;


	//NEP calculation for PnET-II
	share->NEP = share->TotPsn - share->WoodMRespYr - share->WoodGRespYr - share->FolGRespYr - share->SoilRespYr;
	// save current foliar N for output
	share->FolNConOld = veg->FolNCon;

	//calculate mortality after 1990
	int yr;
	double mortality;
	yr = clim->year[rstep];
	if (yr > 1991) {
		mortality = (5e-30 * std::exp(0.0337 * yr))/100;
		if (mortality > veg->WoodTurnover) {
			veg->WoodTurnover = mortality;
		}
	}

	// PnET-CN Only -----------------------------------------------------------------
	if (CN_Mode == 1)
	{
		if (modeltype == 4)
		{
			// In PnET-FUN, HOM/HON are legacy mirrors of FUN-CORPSE pools only.
			SyncLegacySoilPoolsFromFUN(share);
		}

		share->RootNSinkEff = sqrt(1 - (share->PlantN / veg->MaxNStore));   // this is based on monthly rate

		//Annual total variables for PnET-CN
		share->NEP = share->TotPsn - share->SoilDecRespYr - share->WoodDecRespYr - share->WoodMRespYr
			- share->WoodGRespYr - share->FolGRespYr - share->RootMRespYr - share->RootGRespYr;
		share->FolN = (share->FolMass * share->FolNCon / 100.0);
		share->FolC = share->FolMass * veg->CFracBiomass;
		share->TotalN = share->FolN + share->WoodMassN + share->RootMassN + share->HON
			+ share->NH4 + share->NO3 + share->BudN + share->DeadWoodN + share->PlantN;
		share->TotalM = (share->BudC / veg->CFracBiomass) + share->FolMass + (share->WoodMass + share->WoodC / veg->CFracBiomass)
			+ share->RootMass + share->DeadWoodM + share->HOM + (share->PlantC / veg->CFracBiomass) + (share->RootC / veg->CFracBiomass);

		share->NdepTot_prev = share->NdepTot;
	}
	
	//calculate ISA in each year Linghui Meng 20220807
	if (modelmode == 7 && clim->year[rstep] >= site->yrurstart) {share->In = 1;}

	//calculate ISA for next year
	if (share->In>0&&clim->year[rstep] < site->yrurstart2) {

		share->ISA = share->ISA + (share->ISA1 / 120); //calculated ISA in during stage 1
	}
	else
	{
		share->ISA = share->ISA + (share->ISA2 / 30);

	}

}
