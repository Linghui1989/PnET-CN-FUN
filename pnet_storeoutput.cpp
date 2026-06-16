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
		out->decwresp[*ystep] = share->WoodDecRespYr;  //Linghui 0729


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

