#include "pnet_model.h"

void pnet_model::initvars(site_struct* site, veg_struct* veg, share_struct* share)
{
	//
	//Initialize internal shared variable structures for the PnET ecosystem model
	//


	int i;

	//Build "share" structure
	share->Tave = 0;
	share->Tday = 0;
	share->Tnight = 0;
	share->DayLength = 0;
	share->NightLength = 0;
	share->VPD = 0;
	share->dayspan = 1;  // initial for daily
	share->GDDTot = 0;
	share->OldGDDFolEff = 0;
	share->FolLitM = 0;
	share->PosCBalMass = 0;
	share->PosCBalMassTot = 0;
	share->PosCBalMassIx = 0;
	share->LAI = 0;
	share->DVPD = 0;
	share->DayResp = 0;
	share->NightResp = 0;
	share->CanopyNetPsn = 0;
	share->CanopyGrossPsn = 0;
	share->Dwatertot = 0;
	share->DwaterIx = 0;
	share->GrsPsnMo = 0;
	share->NetPsnMo = 0;
	share->FolGRespMo = 0;
	share->WoodMRespYr = 0;
	share->CanopyGrossPsnActMo = 0;
	share->FolProdCYr = 0;
	share->FolProdCMo = 0;
	share->FolGRespYr = 0;
	share->RootProdCYr = 0;
	share->RootMRespYr = 0;
	share->RootGRespYr = 0;
	share->SoilRespMo = 0;
	share->SoilRespYr = 0;
	share->OldGDDWoodEff = 0;
	share->WoodProdCYr = 0;
	share->WoodGRespYr = 0;
	share->TotPsn = 0;
	share->MeanSoilMoistEff = 0;
	share->Drainage = 0;
	share->TotDrain = 0;
	share->TotEvap = 0;
	share->TotTrans = 0;
	share->TotPrec = 0;
	share->TotWater = 0;
	share->TotGrossPsn = 0;
	share->NPPFolYr = 0;
	share->NPPWoodYr = 0;
	share->NPPRootYr = 0;
	share->ET = 0;
	share->NEP = 0;
	share->NetCBal = 0;
	share->SoilDecResp = 0;
	share->BudN = 0;
	share->SoilDecRespYr = 0;
	share->WoodDecRespYr = 0;
	share->DelAmax = 0;
	share->DWUE = 0;
	share->CanopyDO3Pot = 0;
	share->DroughtO3Frac = 0;
	share->TotO3Dose = 0;
	share->RootMassN = 0;
	share->TotalLitterMYr = 0;
	share->TotalLitterNYr = 0;
	share->GrossNImmobYr = 0;
	share->GrossNMinYr = 0;
	share->PlantNUptakeYr = 0;
	share->NetNitrYr = 0;
	share->NetNMinYr = 0;
	share->FracDrain = 0;
	share->NDrainYr = 0;
	share->NDrain = 0;
	share->NDrainMgL = 0;
	share->WoodDecResp = 0;
	share->TotalLitterM = 0;
	share->TotalLitterN = 0;
	share->FolN = 0;
	share->FolC = 0;
	share->TotalN = 0;
	share->TotalM = 0;
	share->NO3 = 0;
	share->NH4 = site->NH4;

	share->freeze = 1;
	share->freezeDoy = 0;

	share->SnowLiq = 0;
	share->SnowNH4 = 0;
	share->SnowNO3 = 0;


	share->FolNCon = veg->FolNCon; //the foliar N conc. at the first year is set in input file
	share->FolNConOld = 0;
	share->NdepTot = 0.0;  // ZZX
	share->NdepTot_prev = 10.0;

	//Shared variables with initial conditions
//	share->FolMass=0;   //share->FolMass=veg->FolMassMin;   In PnET-Day only
//	share->BudC=150;  //share->BudC=(veg->FolMassMax-share->FolMass)*veg->CFracBiomass;  In PnET-Day only
	share->Water = site->SWater;
//	share->DeadWoodM=11300;
//	share->WoodC=300;
//	share->PlantC=900;
	share->RootC=share->WoodC/3;
	share->LightEffMin = 1;
//	share->NRatio=1.3993;
//	share->PlantN=1;
//	share->WoodMass=20000;
//	share->RootMass=6;
	///////////////////////////////////////////////
	share->HOM = site->HOM;
	share->HON = site->HON;
	share->RootNSinkEff = .5;
	share->WUEO3Eff = site->O3EffectOnWUE;
	share->WoodMassN = share->WoodMass * veg->WLPctN * share->NRatio;  //
	share->DeadWoodN = share->DeadWoodM * veg->WLPctN * share->NRatio;
	share->NRatioNit = 0.03;//	share->NRatioNit= share->NRatio;
	share->NetNMinLastYr = 10;
	share->DWater = site->WaterStress;

	share->LightEffCBal = 0;
	share->LightEffCBalTot = 0;
	share->LightEffCBalIx = 0;

	for (i = 0; i < 51; i++)share->O3Effect[i] = 0.0;

	share->avgPCBM = 1.0;					// average light effect
	share->AvgDWater = 1.0;					// average water stress

	share->TaveYr = 0;
	share->PARYr = 0;

	share->WoodProdCMo = 0.0;
	share->RootProdCMo = 0.0;
	
	share->Tsoil = 0.0;
	share->WFPS = 0.0;

	share->RnY2 = 0;
	share->RnX2 = 0;
	share->RnY1 = 0;
	share->RnX1 = 0;
	share->RnMax = 0;

	share->FluxN2OYrDe = 0;
	share->FluxNOYrDe = 0;
	share->FluxN2OYr = 0;
	share->FluxNOYr = 0;
	share->FluxN2Yr = 0;

	share->TotCanopyGrossPsnMax = 0;

	share->O3EffectOnPSN = 0;
	share->ifO3EffectOnPSN = 0;
	share->kO3Eff = 0.0026;   // per ppm
	share->CanopyDO3 = 0;
	share->CanopyDO3Avg = 0;			// Avg O3 effect on photosynthesis for the whole canopy for the whole growing season
	share->CanopyDO3Tot = 0;			// Total O3 effect on photosynthesis for the whole canopy for the whole growing season
	share->CanopyDO3TotDay = 0;			// Total days of O3 effect on photosynthesis for the whole canopy for the whole growing season
	share->kO3EffScalar = 1;   //
	
	share->EnvMaxFol = 0; //Linghui Menng 20200507
	share->SppMaxFol = 0;	//Linghui Meng 20200507
	share->FolGDD = 0;
	share->yrspin = 1000;
	share->distid = 0;

	share->ISA0 = 0;
	share->ISA1 = 0;
	share->ISA2 = 0;
	share->In = 0;
	share->ISA = 0;

	// UptakeEff is root N uptake efficiency
	share->UptakeEff = 1;
	share->PlantNUptakeLast = 0;
	share->FUNPotentialRootProdC = 0;
	share->FUNPotentialWoodProdC = 0;
	share->RootExudateCYr = 0;
	share->MycorrhizalCYr = 0;
	share->MycorrhizalNToPlantYr = 0;
	share->FUNFixNYr = 0;
	share->FUNNdemandYr = 0;
	share->FUNNdemandGapYr = 0;
	share->FUNPotentialRootProdCYr = 0;
	share->FUNPotentialWoodProdCYr = 0;
	share->FUNNLimitRootProdCYr = 0;
	share->FUNNLimitWoodProdCYr = 0;
	share->FUNPlantNOverflowYr = 0;
	share->FUNOverflowToNH4Yr = 0;
	share->FUNOverflowToNO3Yr = 0;
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
	share->BGMicrobialCNLast = 0;
	share->AnnualStartTotalN = 0;
	share->AnnualStartSoilOrgN = 0;
	share->AnnualStartMineralN = 0;
	share->AnnualStartPlantStoreN = 0;
	share->AnnualStartVegStructN = 0;
	share->AnnualStartDeadWoodN = 0;
	share->FolLitN = 0;
	share->RootLitM = 0;
	share->RootLitN = 0;
	share->WoodLitM = 0;
	share->WoodLitN = 0;



	// the input replaces the initials
	site->distyrs = 0;
	site->yrurstart = 1900;
	site->yrurstart2 = 2020;



	veg->FastFlowFrac = site->FastFlowFrac;
	veg->f = site->WaterReleaseFact;
	veg->Kho = site->Kho;
	veg->FolRespFrac = 1;  // based on our previous analysis
	veg->SoilRespA = site->SoilRespA;
	veg->SoilRespB = site->SoilRespB;
	veg->SoilMoistFact = site->SoilMoistFact;
	veg->NImmobA = site->NImmobA;
	veg->NImmobB = site->NImmobB;

	InitializeFUN_CORPSEBelowground(site, share);
	share->FolN = share->FolMass * share->FolNCon / 100.0;
	share->AnnualStartSoilOrgN = share->HON;
	share->AnnualStartMineralN = share->NH4 + share->NO3;
	share->AnnualStartPlantStoreN = share->PlantN + share->BudN;
	share->AnnualStartVegStructN = share->FolN + share->WoodMassN + share->RootMassN;
	share->AnnualStartDeadWoodN = share->DeadWoodN;
	share->AnnualStartTotalN = share->AnnualStartSoilOrgN + share->AnnualStartMineralN
		+ share->AnnualStartPlantStoreN + share->AnnualStartVegStructN + share->AnnualStartDeadWoodN;
}


void pnet_model::init_out(int nyrs, out_struct* out)
{
	//allocate memory for the output structure

	int nmonth = nyrs * 12;
	int i;

	
	// intialize arrays
	for (i = 1; i<nmonth + 1; i++)
	{
		out->grosspsn[i] = 0.0;
		out->netpsn[i] = 0.0;
		out->netcbal[i] = 0.0;
		out->vpd[i] = 0.0;
		out->folmass[i] = 0.0;
		out->plantnMo[i] = 0.0;

	}
	
	for (i = 1; i<nyrs + 1; i++)
	{
		/*
		out->grosspsn[i] = 0.0;
		out->netpsn[i] = 0.0;
		out->netcbal[i] = 0.0;
		out->vpd[i] = 0.0;
		out->folmass[i] = 0.0;
		out->plantnMo[i] = 0.0;
		*/
		out->nppfol[i] = 0.0;
		out->nppwood[i] = 0.0;
		out->npproot[i] = 0.0;
		out->nep[i] = 0.0;
		out->gpp[i] = 0.0;
		out->waterstress[i] = 0.0;
		out->trans[i] = 0.0;
		out->soilwater[i] = 0.0;
		out->psn[i] = 0.0;
		out->drain[i] = 0.0;
		out->prec[i] = 0.0;
		out->evap[i] = 0.0;
		out->et[i] = 0.0;
		out->plantc[i] = 0.0;
		out->budc[i] = 0.0;
		out->woodc[i] = 0.0;
		out->rootc[i] = 0.0;
		out->plantnYr[i] = 0.0;
		out->budn[i] = 0.0;
		out->ndrain[i] = 0.0;
		out->netnmin[i] = 0.0;
		out->grossnmin[i] = 0.0;
		out->nplantuptake[i] = 0.0;
		out->grossnimob[i] = 0.0;
		out->bgnetnmin[i] = 0.0;
		out->orggrossnmin[i] = 0.0;
		out->rhizgrossnmin[i] = 0.0;
		out->bulkgrossnmin[i] = 0.0;
		out->mineralgrossnmin[i] = 0.0;
		out->orgnetnmin_layer[i] = 0.0;
		out->rhiznetnmin_layer[i] = 0.0;
		out->bulknetnmin_layer[i] = 0.0;
		out->mineralnetnmin_layer[i] = 0.0;
		out->rootexudatec[i] = 0.0;
		out->mycnplant[i] = 0.0;
		out->funfixn[i] = 0.0;
		out->littern[i] = 0.0;
		out->netnitrif[i] = 0.0;
		out->nratio[i] = 0.0;
		out->foln[i] = 0.0;
		out->litm[i] = 0.0;
		out->litn[i] = 0.0;
		out->rmresp[i] = 0.0;
		out->rgresp[i] = 0.0;
		out->decresp[i] = 0.0;
		out->totalsoilresp[i] = 0.0;
		out->decwresp[i] = 0.0; //Linghui 0729
		out->orgfastc[i] = 0.0;
		out->orgfastn[i] = 0.0;
		out->orgslowc[i] = 0.0;
		out->orgslown[i] = 0.0;
		out->orgdeadc[i] = 0.0;
		out->orgdeadn[i] = 0.0;
		out->rhizfastc[i] = 0.0;
		out->rhizfastn[i] = 0.0;
		out->rhizslowc[i] = 0.0;
		out->rhizslown[i] = 0.0;
		out->rhizdeadc[i] = 0.0;
		out->rhizdeadn[i] = 0.0;
		out->bulkfastc[i] = 0.0;
		out->bulkfastn[i] = 0.0;
		out->bulkslowc[i] = 0.0;
		out->bulkslown[i] = 0.0;
		out->bulkdeadc[i] = 0.0;
		out->bulkdeadn[i] = 0.0;
		out->orgmicn[i] = 0.0;
		out->rhizmicn[i] = 0.0;
		out->bulkmicn[i] = 0.0;
		out->totalmicn[i] = 0.0;
		out->orgprotn[i] = 0.0;
		out->rhizprotn[i] = 0.0;
		out->bulkprotn[i] = 0.0;
		out->totalprotn[i] = 0.0;
		out->protectionnyr[i] = 0.0;
		out->deprotectionnyr[i] = 0.0;
		out->funndemand[i] = 0.0;
		out->funndemandgap[i] = 0.0;
		out->funpotrootc[i] = 0.0;
		out->funpotwoodc[i] = 0.0;
		out->funnlimitrootc[i] = 0.0;
		out->funnlimitwoodc[i] = 0.0;
		out->funplantnoverflow[i] = 0.0;
		out->funoverflowtonh4[i] = 0.0;
		out->funoverflowtono3[i] = 0.0;
		out->nstart_total[i] = 0.0;
		out->nend_total[i] = 0.0;
		out->dn_total[i] = 0.0;
		out->nstart_soilorg[i] = 0.0;
		out->nend_soilorg[i] = 0.0;
		out->dn_soilorg[i] = 0.0;
		out->nstart_mineral[i] = 0.0;
		out->nend_mineral[i] = 0.0;
		out->dn_mineral[i] = 0.0;
		out->nstart_plantstore[i] = 0.0;
		out->nend_plantstore[i] = 0.0;
		out->dn_plantstore[i] = 0.0;
		out->nstart_vegstruct[i] = 0.0;
		out->nend_vegstruct[i] = 0.0;
		out->dn_vegstruct[i] = 0.0;
		out->nstart_deadwood[i] = 0.0;
		out->nend_deadwood[i] = 0.0;
		out->dn_deadwood[i] = 0.0;
		out->ngas[i] = 0.0;
		out->nbalance_resid[i] = 0.0;
		out->folm[i] = 0.0;
		out->deadwoodm[i] = 0.0;
		out->woodm[i] = 0.0;
		out->rootm[i] = 0.0;
		out->hom[i] = 0.0;
		out->hon[i] = 0.0;
		out->ndep[i] = 0.0;

		out->fn2o[i] = 0.0;
		out->fno[i] = 0.0;
		out->fn2[i] = 0.0;
		out->fn2ode[i] = 0.0;
		out->v1[i] = 0.0;
		out->v2[i] = 0.0;
		out->v3[i] = 0.0;



	}




}
