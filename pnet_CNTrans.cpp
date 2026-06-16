#include "pnet_model.h"

void pnet_model::CNTrans(site_struct* site, veg_struct* veg, clim_struct* clim, int rstep, share_struct* share)
{
	//
	//Carbon and nitrogen translocation routine
	//


	int i;
	double BiomLossFrac, RemoveFrac, RootTurnover, RootLitM, RootLitN, WoodLitM, WoodLitN, WoodMassLoss, WoodTransM;
	double WoodTransN, FolNLoss, Retrans, FolLitN;
	double FolNRetrans;
	//	int dist=1;
	int doyEnd;

	//Check for a disturbance year
	BiomLossFrac = 0;
	RemoveFrac = 0;


	doyEnd = 349;  // default value for monthly

	if (clim->timestep == 1)
	{
		doyEnd = 365;
		//		if (is_leapyear(clim->year[rstep])) doyEnd =366;

	}


	//	if (clim->doy[rstep] > 335)  // need to modify for daily version
	if (clim->doy[rstep] == doyEnd && (site->distyrs > 0))
		//	if (site->distyrs > 0)
	{
		for (i = 1; i <= site->distyrs; i++)
		{
			if (clim->year[rstep] == site->distyear[i])
			{
				BiomLossFrac = site->distintensity[i];  // mortality
				RemoveFrac = site->distremove[i];    // remove out of field
				share->HOM = share->HOM * (1 - site->distsoilloss[i]);
				share->HON = share->HON * (1 - site->distsoilloss[i]);
				break;
			}
		}
	}




	RootTurnover = veg->RootTurnoverA + (veg->RootTurnoverB * share->NetNMinLastYr) + (veg->RootTurnoverC * pow(share->NetNMinLastYr, 2)); //Matlab to C
	if (RootTurnover > 2.5)
	{
		RootTurnover = 2.5;
	}
	if (RootTurnover < 0.1)
	{
		RootTurnover = 0.1;
	}
	RootTurnover = RootTurnover * (share->dayspan / 365.0); //Matlab to C
	if (BiomLossFrac > RootTurnover)
	{
		RootTurnover = BiomLossFrac;
	}
	RootLitM = share->RootMass * RootTurnover;
	RootLitN = RootLitM * (share->RootMassN / share->RootMass);
	//	RootLitN = RootLitM * veg->RLPctN;
	//	share->PlantN = share->PlantN + RootLitN * (share->NRatio-1); //Linghui Meng 20210124
	/*
		if (share->PlantN > veg->MaxNStore) {
			share->NO3 = share->NO3 + (share->PlantN - veg->MaxNStore);
			share->PlantN = veg->MaxNStore;
		}
	*/
	share->RootMass = share->RootMass - RootLitM;
	//	share->RootMassN = share->RootMassN - RootLitN * share->NRatio;
	share->RootMassN = share->RootMassN - RootLitN;

	if (BiomLossFrac > 0)
	{
		WoodLitM = share->WoodMass * BiomLossFrac * (1 - RemoveFrac);
		WoodLitN = share->WoodMassN * BiomLossFrac * (1 - RemoveFrac);
		share->WoodMass = share->WoodMass * (1 - BiomLossFrac);
		share->WoodMassN = share->WoodMassN * (1 - BiomLossFrac);
		share->PlantC *= (1 - BiomLossFrac);
		share->PlantN *= (1 - BiomLossFrac);


	}
	else
	{

		WoodLitM = share->WoodMass * veg->WoodTurnover * (share->dayspan / 365.0); //Matlab to C
//		WoodLitN = WoodLitM * veg->WLPctN;
		WoodLitN = share->WoodMassN * veg->WoodTurnover * (share->dayspan / 365.0); //Matlab to C
		share->WoodMass = share->WoodMass - WoodLitM;
		//		share->PlantN = share->PlantN + WoodLitN * (share->NRatio-1); //Linghui Meng 20210124
		share->WoodMassN = share->WoodMassN - WoodLitN;
		/*
				if (share->PlantN > veg->MaxNStore) {
					share->NO3 = share->NO3 + (share->PlantN - veg->MaxNStore);
					share->PlantN = veg->MaxNStore;
				}
				share->WoodMassN = share->WoodMassN - WoodLitN*share->NRatio; //Linghui Meng 20210124

		*/
	}
	share->DeadWoodM = share->DeadWoodM + WoodLitM;
	share->DeadWoodN = share->DeadWoodN + WoodLitN;
	WoodMassLoss = share->DeadWoodM * veg->WoodLitLossRate * (share->dayspan / 365.0); //Matlab to C
//	if (WoodMassLoss >= 500 * (share->dayspan / 365.0)) { WoodMassLoss = 500 * (share->dayspan / 365.0); } //Linghui Meng 11262020
	WoodTransM = WoodMassLoss * (1 - veg->WoodLitCLoss);  // wood littWoodDecResper
	share->WoodDecResp = (WoodMassLoss - WoodTransM) * veg->CFracBiomass; // loss as CO2
	share->WoodDecRespYr = share->WoodDecRespYr + share->WoodDecResp;
	//	WoodTransN = WoodTransM * veg->WLPctN;
	WoodTransN = (WoodMassLoss / share->DeadWoodM) * share->DeadWoodN;
	//	share->NO3 = share->NO3 + (WoodMassLoss - WoodTransM)*veg->WLPctN; //Linghui MENG 20210124
	share->DeadWoodM = share->DeadWoodM - WoodMassLoss;
	//	share->DeadWoodN = share->DeadWoodN - WoodMassLoss*veg->WLPctN; //Linghui MENG 20210124
	share->DeadWoodN = share->DeadWoodN - WoodTransN;
	share->NetCBal = share->NetCBal - share->WoodDecResp;  // updating NetCBal


	//======================================================
	// add dynamic FolNRetrans algorithm
	FolNRetrans = (69.73 - 34.61 * (1.1 * share->PlantN / veg->MaxNStore)) / 100;
	FolNLoss = share->FolLitM * (share->FolNCon / 100);		
	Retrans = FolNLoss * FolNRetrans;
	share->PlantN += Retrans;
	FolLitN = FolNLoss - Retrans;
	veg->FolNRetrans = FolNRetrans;


	//======================================================
	if (BiomLossFrac > 0)
	{
		share->FolLitM = share->FolLitM + (share->FolMass * BiomLossFrac);
		FolLitN = FolLitN + (share->FolMass * BiomLossFrac * (share->FolNCon / 100));
		share->FolMass = share->FolMass * (1 - BiomLossFrac);
		share->PlantC = share->PlantC * (1 - BiomLossFrac);
		share->PlantN += (veg->MaxNStore - share->PlantN) * BiomLossFrac;  // increase PlantN?
	}

	share->TotalLitterM = share->FolLitM + RootLitM + WoodTransM;
	share->TotalLitterN = FolLitN + RootLitN + WoodTransN;


	share->TotalLitterMYr = share->TotalLitterMYr + share->TotalLitterM;
	share->TotalLitterNYr = share->TotalLitterNYr + share->TotalLitterN;
}

