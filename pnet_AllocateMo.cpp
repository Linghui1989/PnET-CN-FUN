#include "pnet_model.h"

void pnet_model::AllocateMo(veg_struct* veg, share_struct* share, int rstep, int CN_Mode)
{
	//
	// C allocation for the PnET ecosystem model.
	//


	double WoodMRespMo, GDDWoodEff, delGDDWoodEff, WoodProdCMo, WoodGRespMo, TMult, RootCAdd, RootAllocCMo, RootProdCMo, RootMRespMo, RootGRespMo;
	double rpotential, wpotential;
	double folnconnew;



	share->PlantC += share->NetPsnMo - share->FolGRespMo;
//	WoodMRespMo = share->CanopyGrossPsnActMo * veg->WoodMRespA;
	WoodMRespMo = share->GrsPsnMo * veg->WoodMRespA;
	share->WoodMRespYr = share->WoodMRespYr + WoodMRespMo;
	share->FolProdCYr = share->FolProdCYr + share->FolProdCMo;
	share->FolGRespYr = share->FolGRespYr + share->FolGRespMo;


	if (share->FolGDD == 1 && veg->FolReten !=1) {


//		share->BudN = (share->BudC / veg->CFracBiomass) * veg->FLPctN * (1 / (1 - veg->FolNRetrans)) * share->NRatio;//ZZX

		share->BudN = (share->BudC / veg->CFracBiomass) * share->FolNCon / 100; //Linghui 12132025

		if (share->BudN > share->PlantN)
		{
			if (share->PlantN < 0)
			{
				share->BudC *= 0.1;
				share->BudN *= 0.1;
			}
			else
			{
				share->BudC *= (share->PlantN / share->BudN);
				share->BudN *= (share->PlantN / share->BudN);
			}
		}

		folnconnew = (share->FolMass * (veg->FolNCon / 100) + share->BudN) / (share->FolMass + (share->BudC / veg->CFracBiomass)) * 100;
		share->FolNCon = folnconnew;

		share->PlantN = share->PlantN - share->BudN;
		share->FolGDD = 0;
	}


	//calculate root carbon addition
	TMult = (exp(0.1 * (share->Tave - 7.1)) * 0.68) * 1.0;
	RootCAdd = veg->RootAllocA * (share->dayspan / 365.0) + veg->RootAllocB * share->FolProdCMo;
	share->RootC += RootCAdd;
	
	//if plantc is not enough, stop allocate carbon
	if (share->PlantC < 0) {
		RootCAdd = 0;
	}
	share->PlantC -=  RootCAdd;
	//calculate root allocation
	RootAllocCMo = (share->dayspan / 365.0) * TMult;   // modified from matlab version
	if (RootAllocCMo > 1.0) RootAllocCMo = 1.0;
	RootAllocCMo = RootAllocCMo * share->RootC;
	RootProdCMo = RootAllocCMo / (1.0 + veg->RootMRespFrac + veg->GRespFrac);



	//calculate wood production
	if (share->GDDTot >= veg->GDDWoodStart)
	{
		GDDWoodEff = (share->GDDTot - veg->GDDWoodStart) / (veg->GDDWoodEnd - veg->GDDWoodStart);
		if (GDDWoodEff > 1.0)GDDWoodEff = 1;
		if (GDDWoodEff < 0)GDDWoodEff = 0;

		delGDDWoodEff = GDDWoodEff - share->OldGDDWoodEff;
		WoodProdCMo = share->WoodC * delGDDWoodEff;
		share->OldGDDWoodEff = GDDWoodEff;

	}
	else
	{
		WoodProdCMo = 0;
		WoodGRespMo = 0;
	}

	

	if (CN_Mode == 1)
	{
		double PotentialN, WoodN, RootN;

		WoodN = WoodProdCMo / veg->CFracBiomass * veg->RLPctN * share->NRatio;
		RootN = RootProdCMo / veg->CFracBiomass * veg->RLPctN * share->NRatio;
		PotentialN = WoodN + RootN;

		/*
		if (PotentialN > share->PlantN) {

			WoodProdCMo *= (share->PlantN / PotentialN);
			RootProdCMo *= (share->PlantN / PotentialN);
			PotentialN = share->PlantN;
		}
		*/
		//this is a function that detect if there is enough N to allocate
		if (share->PlantN < PotentialN) {
			//if the N is not enough to allocate root, coarse root and wood is 0, all N allocate to root
			if (share->PlantN <= RootN) {

				RootN = share->PlantN;
				WoodN = 0;

			}
			else
			{

				WoodN = share->PlantN - RootN;
			}

		}



		RootProdCMo = RootN * veg->CFracBiomass / (veg->RLPctN * share->NRatio);
		WoodProdCMo = WoodN * veg->CFracBiomass / (veg->WLPctN * share->NRatio);



	}


	share->RootMass = share->RootMass + (RootProdCMo / veg->CFracBiomass);
	share->RootMassN = share->RootMassN + ((RootProdCMo / veg->CFracBiomass) * veg->RLPctN * share->NRatio);
	share->PlantN = share->PlantN - ((RootProdCMo / veg->CFracBiomass) * veg->RLPctN * share->NRatio);

	share->WoodMass = share->WoodMass + (WoodProdCMo / veg->CFracBiomass);
	share->WoodMassN = share->WoodMassN + ((WoodProdCMo / veg->CFracBiomass) * veg->WLPctN * share->NRatio);
	share->PlantN = share->PlantN - ((WoodProdCMo / veg->CFracBiomass) * veg->WLPctN * share->NRatio);




	//substract cost from C pool
	share->RootC -= RootAllocCMo;
	share->WoodC -= WoodProdCMo;
	RootMRespMo = RootProdCMo * veg->RootMRespFrac;
	RootGRespMo = RootProdCMo * veg->GRespFrac;
	WoodGRespMo = WoodProdCMo * veg->GRespFrac;

	share->WoodProdCYr = share->WoodProdCYr + WoodProdCMo;
	share->WoodGRespYr = share->WoodGRespYr + WoodGRespMo;

	//substract resp from plant C
	//if wood resp greater than plant c, extra c extract from wood c pool
	if (share->PlantC > WoodGRespMo + WoodMRespMo) {

		share->PlantC -= WoodGRespMo + WoodMRespMo;
		
	}
	else
	{
		share->WoodC -= (WoodGRespMo + WoodMRespMo - share->PlantC);
		share->PlantC = 0;
	}
	


	share->WoodProdCMo = WoodProdCMo;
//	share->PlantC = share->PlantC - WoodMRespMo - WoodGRespMo;

	share->RootProdCYr = share->RootProdCYr + RootProdCMo;
	share->RootMRespYr = share->RootMRespYr + RootMRespMo;
	share->RootGRespYr = share->RootGRespYr + RootGRespMo;


		//PnET-CN Only -----------------------------------------------------------------
	if (CN_Mode == 1)
	{
		share->NetCBal = share->NetCBal + share->NetPsnMo - WoodMRespMo - WoodGRespMo - share->FolGRespMo - RootMRespMo - RootGRespMo;
		// needs -share->SoilDecResp - share->WoodDecResp, and will be updated in the respective routine.

	}
	else
	{
		//------------------------------------------------------------------------------
		share->NetCBal = share->NetCBal + share->NetPsnMo - share->SoilRespMo - WoodMRespMo - WoodGRespMo - share->FolGRespMo;
	}

}

