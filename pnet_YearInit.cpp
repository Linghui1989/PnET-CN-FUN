#include "pnet_model.h"

void pnet_model::YearInit(share_struct* share)
{
	//
	// Annual initialization for the PnET ecosystem model.
	//

	share->GDDTot = 0;
	share->WoodMRespYr = 0;
	share->SoilRespYr = 0;
	share->TotTrans = 0;
	share->TotPsn = 0;
	share->TotGrossPsn = 0;
	share->Drainage = 0;//Linghui 0527
	share->TotDrain = 0;
	share->TotPrec = 0;
	share->TotEvap = 0;
	share->TotWater =0;
	share->FolProdCYr = 0;
	share->WoodProdCYr = 0;
	share->RootProdCYr = 0;
	share->RootMRespYr = 0;
	share->FolGRespYr = 0;
	share->WoodGRespYr = 0;
	share->RootGRespYr = 0;
	share->OldGDDFolEff = 0;
	share->OldGDDWoodEff = 0;
	share->PosCBalMassTot = 0;
	share->PosCBalMassIx = 0;
	share->Dwatertot = 0;
	share->DwaterIx = 0;
	share->NDrainYr=0;
	share->NetCBal = 0;
	share->GrossNMinYr=0;
	share->PlantNUptakeYr=0;
	share->GrossNImmobYr=0;
	share->TotalLitterMYr=0;
	share->TotalLitterNYr=0;
	share->NetNitrYr=0;
	share->LightEffMin=1;
	share->SoilDecRespYr=0;
	share->WoodDecRespYr=0;
	share->NetNMinLastYr = share->NetNMinYr;
	share->NetNMinYr=0;
	share->NdepTot=0.0; //ZZX
//	share->NdepTot_prev = 0.0;
	share->LightEffCBalTot=0;
	share->LightEffCBalIx=0;

	share->TaveYr=0;
	share->PARYr=0;

	share->FluxN2OYrDe=0;
	share->FluxNOYrDe=0;
	share->FluxN2OYr=0;
	share->FluxNOYr=0;
	share->FluxN2Yr=0;

	share->TotCanopyGrossPsnMax=0;
	share-> CanopyDO3Tot = 0;			// Total O3 effect on photosynthesis for the whole canopy for the whole growing season
	share-> CanopyDO3TotDay =0;			// Total days of O3 effect on photosynthesis for the whole canopy for the whole growing season
	share->EnvMaxFol = 0;
	share->SppMaxFol = 0;
	for (int i = 0; i<51; i++)share->O3Effect[i]=0.0;

	share->TotO3Dose=0;

	share->CDDTot = 0;
	share->RootExudateCYr = 0;
	share->MycorrhizalCYr = 0;
	share->MycorrhizalNToPlantYr = 0;
	share->FUNFixNYr = 0;
	share->BGGrossNMobilizationYr = 0;
	share->BGGrossNImmobilizationYr = 0;
	share->BGNetNMineralizationYr = 0;
	share->BGOrgGrossNMineralizationYr = 0;
	share->BGRhizGrossNMineralizationYr = 0;
	share->BGBulkGrossNMineralizationYr = 0;
	share->BGOrgGrossNImmobilizationYr = 0;
	share->BGRhizGrossNImmobilizationYr = 0;
	share->BGBulkGrossNImmobilizationYr = 0;
	share->BGOrgNetNMineralizationYr = 0;
	share->BGRhizNetNMineralizationYr = 0;
	share->BGBulkNetNMineralizationYr = 0;
	share->BGMicrobialTurnoverNYr = 0;
	share->BGProtectionNYr = 0;
	share->BGDeprotectionNYr = 0;
	share->PlantNUptakeLast = 0;
	share->FUNPotentialRootProdC = 0;
	share->FUNPotentialWoodProdC = 0;
	share->FUNNdemandYr = 0;
	share->FUNNdemandGapYr = 0;
	share->FUNPotentialRootProdCYr = 0;
	share->FUNPotentialWoodProdCYr = 0;
	share->FUNNLimitRootProdCYr = 0;
	share->FUNNLimitWoodProdCYr = 0;
	share->FUNPlantNOverflowYr = 0;
	share->FUNOverflowToNH4Yr = 0;
	share->FUNOverflowToNO3Yr = 0;

	share->FolN = share->FolMass * share->FolNCon / 100.0;
	share->AnnualStartSoilOrgN = share->HON;
	share->AnnualStartMineralN = share->NH4 + share->NO3;
	share->AnnualStartPlantStoreN = share->PlantN + share->BudN;
	share->AnnualStartVegStructN = share->FolN + share->WoodMassN + share->RootMassN;
	share->AnnualStartDeadWoodN = share->DeadWoodN;
	share->AnnualStartTotalN = share->AnnualStartSoilOrgN + share->AnnualStartMineralN
		+ share->AnnualStartPlantStoreN + share->AnnualStartVegStructN + share->AnnualStartDeadWoodN;

}
