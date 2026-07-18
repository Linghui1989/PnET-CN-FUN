#include "pnet_model.h"

void pnet_model::storeoutput(veg_struct* veg, share_struct* share, out_struct* out, int rstep, int* ystep, int NewYear)
{
	//
	// Add variables to the returned output structure so that the user may work
	// with them (or save them) at the command line after running the model.
	//


	//**************************************************************************
	//ystep is replaced with *ystep during the transformation from Matlab to C
	//**************************************************************************

	// Store iteration step variables (these may be monthly, daily etc. based on
	// the stepping of the input climate data)
	if (NewYear == 0)
	{
		out->grosspsn[rstep] = share->GrsPsnMo;
		out->netpsn[rstep] = share->NetPsnMo;
		out->netcbal[rstep] = share->NetCBal;
		out->vpd[rstep] = share->VPD;
		out->folmass[rstep] = share->FolMass;
		out->plantnMo[rstep] = share->PlantN;
	}

	// Store annual variables at the conclusion of each years run
	if (NewYear == 1)
	{

		out->nppfol[*ystep] = share->NPPFolYr;
		out->nppwood[*ystep] = share->NPPWoodYr;
		out->npproot[*ystep] = share->NPPRootYr;
		out->nep[*ystep] = share->NetCBal;
		out->gpp[*ystep] = share->TotGrossPsn;
		out->psn[*ystep] = share->TotPsn;  // total net psn

		// Water cycle
		out->waterstress[*ystep] = share->AvgDWater;
		out->trans[*ystep] = share->TotTrans;
		out->soilwater[*ystep] = share->TotWater/365; //Linghui 052823
		out->drain[*ystep] = share->TotDrain;
		out->prec[*ystep] = share->TotPrec;
		out->evap[*ystep] = share->TotEvap;
	  //out->et[*ystep] = share->ET;
		out->et[*ystep] = share->TotEvap + share->TotTrans;//ZZX

		// Carbon cycle
		out->plantc[*ystep] = share->PlantC;
		out->budc[*ystep] = share->BudC;
		out->woodc[*ystep] = share->WoodC;
		out->rootc[*ystep] = share->RootC;

		out->folm[*ystep] = share->FolMass;
		out->deadwoodm[*ystep] = share->DeadWoodM;
		out->woodm[*ystep] = share->WoodMass;
		out->rootm[*ystep] = share->RootMass;
		out->hom[*ystep] = share->HOM;
		out->hon[*ystep] = share->HON;
		out->ndep[*ystep] = share->NdepTot;


		// Nitrogen cycle
		out->plantnYr[*ystep] = share->PlantN;
		out->budn[*ystep] = share->BudN;
		out->ndrain[*ystep] = share->NDrainYr;
		out->netnmin[*ystep] = share->NetNMinYr;
		out->grossnmin[*ystep] = share->GrossNMinYr;
		out->nplantuptake[*ystep] = share->PlantNUptakeYr;
		out->grossnimob[*ystep] = share->GrossNImmobYr;
		out->bgnetnmin[*ystep] = share->BGNetNMineralizationYr;
		out->orggrossnmin[*ystep] = share->BGOrgGrossNMineralizationYr;
		out->rhizgrossnmin[*ystep] = share->BGRhizGrossNMineralizationYr;
		out->bulkgrossnmin[*ystep] = share->BGBulkGrossNMineralizationYr;
		out->mineralgrossnmin[*ystep] = share->BGRhizGrossNMineralizationYr + share->BGBulkGrossNMineralizationYr;
		out->orgnetnmin_layer[*ystep] = share->BGOrgNetNMineralizationYr;
		out->rhiznetnmin_layer[*ystep] = share->BGRhizNetNMineralizationYr;
		out->bulknetnmin_layer[*ystep] = share->BGBulkNetNMineralizationYr;
		out->mineralnetnmin_layer[*ystep] = share->BGRhizNetNMineralizationYr + share->BGBulkNetNMineralizationYr;
		out->rootexudatec[*ystep] = share->RootExudateCYr;
		out->mycnplant[*ystep] = share->MycorrhizalNToPlantYr;
		out->funfixn[*ystep] = share->FUNFixNYr;
		out->littern[*ystep] = share->TotalLitterNYr;
		out->netnitrif[*ystep] = share->NetNitrYr;
		out->nratio[*ystep] = share->NRatio;
		out->foln[*ystep] = share->FolNCon;  //

												// TBCA
		out->litm[*ystep] = share->TotalLitterMYr;
		out->litn[*ystep] = share->TotalLitterNYr;
		out->rmresp[*ystep] = share->RootMRespYr;
		out->rgresp[*ystep] = share->RootGRespYr;
		out->decresp[*ystep] = share->SoilDecRespYr;
		out->totalsoilresp[*ystep] = share->RootMRespYr + share->SoilDecRespYr;
		out->decwresp[*ystep] = share->WoodDecRespYr;  //Linghui 0729
		out->orgfastc[*ystep] = share->funBG.organicLayer.unprotC[BG_FAST] + share->funBG.organicLayer.protC[BG_FAST];
		out->orgfastn[*ystep] = share->funBG.organicLayer.unprotN[BG_FAST] + share->funBG.organicLayer.protN[BG_FAST];
		out->orgslowc[*ystep] = share->funBG.organicLayer.unprotC[BG_SLOW] + share->funBG.organicLayer.protC[BG_SLOW];
		out->orgslown[*ystep] = share->funBG.organicLayer.unprotN[BG_SLOW] + share->funBG.organicLayer.protN[BG_SLOW];
		out->orgdeadc[*ystep] = share->funBG.organicLayer.unprotC[BG_DEADMIC] + share->funBG.organicLayer.protC[BG_DEADMIC];
		out->orgdeadn[*ystep] = share->funBG.organicLayer.unprotN[BG_DEADMIC] + share->funBG.organicLayer.protN[BG_DEADMIC];
		out->rhizfastc[*ystep] = share->funBG.mineralRhizosphere.unprotC[BG_FAST] + share->funBG.mineralRhizosphere.protC[BG_FAST];
		out->rhizfastn[*ystep] = share->funBG.mineralRhizosphere.unprotN[BG_FAST] + share->funBG.mineralRhizosphere.protN[BG_FAST];
		out->rhizslowc[*ystep] = share->funBG.mineralRhizosphere.unprotC[BG_SLOW] + share->funBG.mineralRhizosphere.protC[BG_SLOW];
		out->rhizslown[*ystep] = share->funBG.mineralRhizosphere.unprotN[BG_SLOW] + share->funBG.mineralRhizosphere.protN[BG_SLOW];
		out->rhizdeadc[*ystep] = share->funBG.mineralRhizosphere.unprotC[BG_DEADMIC] + share->funBG.mineralRhizosphere.protC[BG_DEADMIC];
		out->rhizdeadn[*ystep] = share->funBG.mineralRhizosphere.unprotN[BG_DEADMIC] + share->funBG.mineralRhizosphere.protN[BG_DEADMIC];
		out->bulkfastc[*ystep] = share->funBG.mineralBulk.unprotC[BG_FAST] + share->funBG.mineralBulk.protC[BG_FAST];
		out->bulkfastn[*ystep] = share->funBG.mineralBulk.unprotN[BG_FAST] + share->funBG.mineralBulk.protN[BG_FAST];
		out->bulkslowc[*ystep] = share->funBG.mineralBulk.unprotC[BG_SLOW] + share->funBG.mineralBulk.protC[BG_SLOW];
		out->bulkslown[*ystep] = share->funBG.mineralBulk.unprotN[BG_SLOW] + share->funBG.mineralBulk.protN[BG_SLOW];
		out->bulkdeadc[*ystep] = share->funBG.mineralBulk.unprotC[BG_DEADMIC] + share->funBG.mineralBulk.protC[BG_DEADMIC];
		out->bulkdeadn[*ystep] = share->funBG.mineralBulk.unprotN[BG_DEADMIC] + share->funBG.mineralBulk.protN[BG_DEADMIC];
		out->orgmicn[*ystep] = share->funBG.organicLayer.livingMicrobeN;
		out->rhizmicn[*ystep] = share->funBG.mineralRhizosphere.livingMicrobeN;
		out->bulkmicn[*ystep] = share->funBG.mineralBulk.livingMicrobeN;
		out->totalmicn[*ystep] = share->funBG.organicLayer.livingMicrobeN
			+ share->funBG.mineralRhizosphere.livingMicrobeN
			+ share->funBG.mineralBulk.livingMicrobeN;
		out->orgprotn[*ystep] = share->funBG.organicLayer.protN[BG_FAST]
			+ share->funBG.organicLayer.protN[BG_SLOW]
			+ share->funBG.organicLayer.protN[BG_DEADMIC];
		out->rhizprotn[*ystep] = share->funBG.mineralRhizosphere.protN[BG_FAST]
			+ share->funBG.mineralRhizosphere.protN[BG_SLOW]
			+ share->funBG.mineralRhizosphere.protN[BG_DEADMIC];
		out->bulkprotn[*ystep] = share->funBG.mineralBulk.protN[BG_FAST]
			+ share->funBG.mineralBulk.protN[BG_SLOW]
			+ share->funBG.mineralBulk.protN[BG_DEADMIC];
		out->totalprotn[*ystep] = out->orgprotn[*ystep] + out->rhizprotn[*ystep] + out->bulkprotn[*ystep];
		out->protectionnyr[*ystep] = share->BGProtectionNYr;
		out->deprotectionnyr[*ystep] = share->BGDeprotectionNYr;
		out->funndemand[*ystep] = share->FUNNdemandYr;
		out->funndemandgap[*ystep] = share->FUNNdemandGapYr;
		out->funpotrootc[*ystep] = share->FUNPotentialRootProdCYr;
		out->funpotwoodc[*ystep] = share->FUNPotentialWoodProdCYr;
		out->funnlimitrootc[*ystep] = share->FUNNLimitRootProdCYr;
		out->funnlimitwoodc[*ystep] = share->FUNNLimitWoodProdCYr;
		out->funplantnoverflow[*ystep] = share->FUNPlantNOverflowYr;
		out->funoverflowtonh4[*ystep] = share->FUNOverflowToNH4Yr;
		out->funoverflowtono3[*ystep] = share->FUNOverflowToNO3Yr;
		out->nstart_total[*ystep] = share->AnnualStartTotalN;
		out->nend_total[*ystep] = share->TotalN;
		out->dn_total[*ystep] = out->nend_total[*ystep] - out->nstart_total[*ystep];
		out->nstart_soilorg[*ystep] = share->AnnualStartSoilOrgN;
		out->nend_soilorg[*ystep] = share->HON;
		out->dn_soilorg[*ystep] = out->nend_soilorg[*ystep] - out->nstart_soilorg[*ystep];
		out->nstart_mineral[*ystep] = share->AnnualStartMineralN;
		out->nend_mineral[*ystep] = share->NH4 + share->NO3;
		out->dn_mineral[*ystep] = out->nend_mineral[*ystep] - out->nstart_mineral[*ystep];
		out->nstart_plantstore[*ystep] = share->AnnualStartPlantStoreN;
		out->nend_plantstore[*ystep] = share->PlantN + share->BudN;
		out->dn_plantstore[*ystep] = out->nend_plantstore[*ystep] - out->nstart_plantstore[*ystep];
		out->nstart_vegstruct[*ystep] = share->AnnualStartVegStructN;
		out->nend_vegstruct[*ystep] = share->FolN + share->WoodMassN + share->RootMassN;
		out->dn_vegstruct[*ystep] = out->nend_vegstruct[*ystep] - out->nstart_vegstruct[*ystep];
		out->nstart_deadwood[*ystep] = share->AnnualStartDeadWoodN;
		out->nend_deadwood[*ystep] = share->DeadWoodN;
		out->dn_deadwood[*ystep] = out->nend_deadwood[*ystep] - out->nstart_deadwood[*ystep];
		out->ngas[*ystep] = share->FluxN2OYr + share->FluxNOYr + share->FluxN2Yr;
		out->nbalance_resid[*ystep] = out->dn_total[*ystep]
			- (share->NdepTot - share->NDrainYr - out->ngas[*ystep]);


		out->fn2o[*ystep] = share->FluxN2OYr;
		out->fno[*ystep] = share->FluxNOYr;
		out->fn2[*ystep] = share->FluxN2Yr;

		out->fn2ode[*ystep] = share->TotalN;
		out->fnode[*ystep] = veg->FolNRetrans;	//Linghui 12142025


		out->v1[*ystep] = share->WoodMassN;
		out->v2[*ystep] = share->DeadWoodN;
		out->v3[*ystep] = veg->WoodTurnover;

		
		if (share->In > 0) {



			out->folm[*ystep] = out->folm[*ystep] * (1 - share->ISA);
			out->woodm[*ystep] = out->woodm[*ystep] * (1 - share->ISA);
			out->rootm[*ystep] = out->rootm[*ystep] * (1 - share->ISA);
			out->deadwoodm[*ystep] = out->deadwoodm[*ystep] * (1 - share->ISA);
			out->hom[*ystep] = out->hom[*ystep] * (1 - share->ISA);
			out->litm[*ystep] = out->litm[*ystep] * (1 - share->ISA);
			out->litn[*ystep] = out->litn[*ystep] * (1 - share->ISA);
			out->nppfol[*ystep] = out->nppfol[*ystep] * (1 - share->ISA);
			out->nppwood[*ystep] = out->nppwood[*ystep] * (1 - share->ISA);
			out->npproot[*ystep] = out->npproot[*ystep] * (1 - share->ISA);
			out->gpp[*ystep] = out->gpp[*ystep] * (1 - share->ISA);
			out->nep[*ystep] = out->nep[*ystep] * (1 - share->ISA);
			out->psn[*ystep] = out->psn[*ystep] * (1 - share->ISA);
			out->hon[*ystep] = out->hon[*ystep] * (1 - share->ISA);

			out->v1[*ystep] = out->v1[*ystep] * (1 - share->ISA);
			out->v2[*ystep] = out->v2[*ystep] * (1 - share->ISA);
			out->v3[*ystep] = out->v3[*ystep] * (1 - share->ISA);

		}
		
		
		//advance year counter

		*ystep = *ystep + 1;

	}
}

