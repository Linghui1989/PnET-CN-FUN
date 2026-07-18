#include "pnet_model.h"

void pnet_model::memset_out(int nyrs, out_struct* out)
{
	//allocate memory for the output structure

	int nmonth = nyrs * 12;
//	int i;

	out->grosspsn = (double*)malloc((nmonth + 1) * sizeof(double)); //use nmonth+1 so that rstep can start from 1
	out->netpsn = (double*)malloc((nmonth + 1) * sizeof(double));
	out->netcbal = (double*)malloc((nmonth + 1) * sizeof(double));
	out->vpd = (double*)malloc((nmonth + 1) * sizeof(double));
	out->folmass = (double*)malloc((nmonth + 1) * sizeof(double));
	out->plantnMo = (double*)malloc((nmonth + 1) * sizeof(double));

	out->nppfol = (double*)malloc((nyrs + 1) * sizeof(double)); //use (nyrs+1)+1 so that ystep can start from 1
	out->nppwood = (double*)malloc((nyrs + 1) * sizeof(double));
	out->npproot = (double*)malloc((nyrs + 1) * sizeof(double));
	out->nep = (double*)malloc((nyrs + 1) * sizeof(double));
	out->gpp = (double*)malloc((nyrs + 1) * sizeof(double));
	out->waterstress = (double*)malloc((nyrs + 1) * sizeof(double));
	out->trans = (double*)malloc((nyrs + 1) * sizeof(double));
	out->soilwater = (double*)malloc((nyrs + 1) * sizeof(double));
	out->psn = (double*)malloc((nyrs + 1) * sizeof(double));
	out->drain = (double*)malloc((nyrs + 1) * sizeof(double));
	out->prec = (double*)malloc((nyrs + 1) * sizeof(double));
	out->evap = (double*)malloc((nyrs + 1) * sizeof(double));
	out->et = (double*)malloc((nyrs + 1) * sizeof(double));
	out->plantc = (double*)malloc((nyrs + 1) * sizeof(double));
	out->budc = (double*)malloc((nyrs + 1) * sizeof(double));
	out->woodc = (double*)malloc((nyrs + 1) * sizeof(double));
	out->rootc = (double*)malloc((nyrs + 1) * sizeof(double));
	out->plantnYr = (double*)malloc((nyrs + 1) * sizeof(double));
	out->budn = (double*)malloc((nyrs + 1) * sizeof(double));
	out->ndrain = (double*)malloc((nyrs + 1) * sizeof(double));
	out->netnmin = (double*)malloc((nyrs + 1) * sizeof(double));
	out->grossnmin = (double*)malloc((nyrs + 1) * sizeof(double));
	out->nplantuptake = (double*)malloc((nyrs + 1) * sizeof(double));
	out->grossnimob = (double*)malloc((nyrs + 1) * sizeof(double));
	out->bgnetnmin = (double*)malloc((nyrs + 1) * sizeof(double));
	out->orggrossnmin = (double*)malloc((nyrs + 1) * sizeof(double));
	out->rhizgrossnmin = (double*)malloc((nyrs + 1) * sizeof(double));
	out->bulkgrossnmin = (double*)malloc((nyrs + 1) * sizeof(double));
	out->mineralgrossnmin = (double*)malloc((nyrs + 1) * sizeof(double));
	out->orgnetnmin_layer = (double*)malloc((nyrs + 1) * sizeof(double));
	out->rhiznetnmin_layer = (double*)malloc((nyrs + 1) * sizeof(double));
	out->bulknetnmin_layer = (double*)malloc((nyrs + 1) * sizeof(double));
	out->mineralnetnmin_layer = (double*)malloc((nyrs + 1) * sizeof(double));
	out->rootexudatec = (double*)malloc((nyrs + 1) * sizeof(double));
	out->mycnplant = (double*)malloc((nyrs + 1) * sizeof(double));
	out->funfixn = (double*)malloc((nyrs + 1) * sizeof(double));
	out->littern = (double*)malloc((nyrs + 1) * sizeof(double));
	out->netnitrif = (double*)malloc((nyrs + 1) * sizeof(double));
	out->nratio = (double*)malloc((nyrs + 1) * sizeof(double));
	out->foln = (double*)malloc((nyrs + 1) * sizeof(double));
	out->litm = (double*)malloc((nyrs + 1) * sizeof(double));
	out->litn = (double*)malloc((nyrs + 1) * sizeof(double));
	out->rmresp = (double*)malloc((nyrs + 1) * sizeof(double));
	out->rgresp = (double*)malloc((nyrs + 1) * sizeof(double));
	out->decresp = (double*)malloc((nyrs + 1) * sizeof(double));
	out->totalsoilresp = (double*)malloc((nyrs + 1) * sizeof(double));
	out->decwresp = (double*)malloc((nyrs + 1) * sizeof(double)); //Linghui 0729
	out->orgfastc = (double*)malloc((nyrs + 1) * sizeof(double));
	out->orgfastn = (double*)malloc((nyrs + 1) * sizeof(double));
	out->orgslowc = (double*)malloc((nyrs + 1) * sizeof(double));
	out->orgslown = (double*)malloc((nyrs + 1) * sizeof(double));
	out->orgdeadc = (double*)malloc((nyrs + 1) * sizeof(double));
	out->orgdeadn = (double*)malloc((nyrs + 1) * sizeof(double));
	out->rhizfastc = (double*)malloc((nyrs + 1) * sizeof(double));
	out->rhizfastn = (double*)malloc((nyrs + 1) * sizeof(double));
	out->rhizslowc = (double*)malloc((nyrs + 1) * sizeof(double));
	out->rhizslown = (double*)malloc((nyrs + 1) * sizeof(double));
	out->rhizdeadc = (double*)malloc((nyrs + 1) * sizeof(double));
	out->rhizdeadn = (double*)malloc((nyrs + 1) * sizeof(double));
	out->bulkfastc = (double*)malloc((nyrs + 1) * sizeof(double));
	out->bulkfastn = (double*)malloc((nyrs + 1) * sizeof(double));
	out->bulkslowc = (double*)malloc((nyrs + 1) * sizeof(double));
	out->bulkslown = (double*)malloc((nyrs + 1) * sizeof(double));
	out->bulkdeadc = (double*)malloc((nyrs + 1) * sizeof(double));
	out->bulkdeadn = (double*)malloc((nyrs + 1) * sizeof(double));
	out->orgmicn = (double*)malloc((nyrs + 1) * sizeof(double));
	out->rhizmicn = (double*)malloc((nyrs + 1) * sizeof(double));
	out->bulkmicn = (double*)malloc((nyrs + 1) * sizeof(double));
	out->totalmicn = (double*)malloc((nyrs + 1) * sizeof(double));
	out->orgprotn = (double*)malloc((nyrs + 1) * sizeof(double));
	out->rhizprotn = (double*)malloc((nyrs + 1) * sizeof(double));
	out->bulkprotn = (double*)malloc((nyrs + 1) * sizeof(double));
	out->totalprotn = (double*)malloc((nyrs + 1) * sizeof(double));
	out->protectionnyr = (double*)malloc((nyrs + 1) * sizeof(double));
	out->deprotectionnyr = (double*)malloc((nyrs + 1) * sizeof(double));
	out->funndemand = (double*)malloc((nyrs + 1) * sizeof(double));
	out->funndemandgap = (double*)malloc((nyrs + 1) * sizeof(double));
	out->funpotrootc = (double*)malloc((nyrs + 1) * sizeof(double));
	out->funpotwoodc = (double*)malloc((nyrs + 1) * sizeof(double));
	out->funnlimitrootc = (double*)malloc((nyrs + 1) * sizeof(double));
	out->funnlimitwoodc = (double*)malloc((nyrs + 1) * sizeof(double));
	out->funplantnoverflow = (double*)malloc((nyrs + 1) * sizeof(double));
	out->funoverflowtonh4 = (double*)malloc((nyrs + 1) * sizeof(double));
	out->funoverflowtono3 = (double*)malloc((nyrs + 1) * sizeof(double));
	out->nstart_total = (double*)malloc((nyrs + 1) * sizeof(double));
	out->nend_total = (double*)malloc((nyrs + 1) * sizeof(double));
	out->dn_total = (double*)malloc((nyrs + 1) * sizeof(double));
	out->nstart_soilorg = (double*)malloc((nyrs + 1) * sizeof(double));
	out->nend_soilorg = (double*)malloc((nyrs + 1) * sizeof(double));
	out->dn_soilorg = (double*)malloc((nyrs + 1) * sizeof(double));
	out->nstart_mineral = (double*)malloc((nyrs + 1) * sizeof(double));
	out->nend_mineral = (double*)malloc((nyrs + 1) * sizeof(double));
	out->dn_mineral = (double*)malloc((nyrs + 1) * sizeof(double));
	out->nstart_plantstore = (double*)malloc((nyrs + 1) * sizeof(double));
	out->nend_plantstore = (double*)malloc((nyrs + 1) * sizeof(double));
	out->dn_plantstore = (double*)malloc((nyrs + 1) * sizeof(double));
	out->nstart_vegstruct = (double*)malloc((nyrs + 1) * sizeof(double));
	out->nend_vegstruct = (double*)malloc((nyrs + 1) * sizeof(double));
	out->dn_vegstruct = (double*)malloc((nyrs + 1) * sizeof(double));
	out->nstart_deadwood = (double*)malloc((nyrs + 1) * sizeof(double));
	out->nend_deadwood = (double*)malloc((nyrs + 1) * sizeof(double));
	out->dn_deadwood = (double*)malloc((nyrs + 1) * sizeof(double));
	out->ngas = (double*)malloc((nyrs + 1) * sizeof(double));
	out->nbalance_resid = (double*)malloc((nyrs + 1) * sizeof(double));


	out->folm = (double*)malloc((nyrs + 1) * sizeof(double));
	out->deadwoodm = (double*)malloc((nyrs + 1) * sizeof(double));
	out->woodm = (double*)malloc((nyrs + 1) * sizeof(double));
	out->rootm = (double*)malloc((nyrs + 1) * sizeof(double));
	out->hom = (double*)malloc((nyrs + 1) * sizeof(double));
	out->hon = (double*)malloc((nyrs + 1) * sizeof(double));
	out->ndep = (double*)malloc((nyrs + 1) * sizeof(double));


	out->fn2o = (double*)malloc((nyrs + 1) * sizeof(double));
	out->fno = (double*)malloc((nyrs + 1) * sizeof(double));
	out->fn2 = (double*)malloc((nyrs + 1) * sizeof(double));
	out->fn2ode = (double*)malloc((nyrs + 1) * sizeof(double));
	out->fnode = (double*)malloc((nyrs + 1) * sizeof(double));


	out->v1 = (double*)malloc((nyrs + 1) * sizeof(double));
	out->v2 = (double*)malloc((nyrs + 1) * sizeof(double));
	out->v3 = (double*)malloc((nyrs + 1) * sizeof(double));




	if (!out->grosspsn || !out->netpsn || !out->netcbal || !out->vpd || !out->folmass || !out->plantnMo || !out->nppfol || !out->nppwood ||
		!out->npproot || !out->nep || !out->gpp || !out->waterstress || !out->trans || !out->soilwater || !out->psn || !out->drain || !out->prec ||
		!out->evap || !out->et || !out->plantc || !out->budc || !out->woodc || !out->rootc || !out->plantnYr || !out->budn || !out->ndrain ||
		!out->netnmin || !out->grossnmin || !out->nplantuptake || !out->grossnimob || !out->bgnetnmin ||
		!out->orggrossnmin || !out->rhizgrossnmin || !out->bulkgrossnmin || !out->mineralgrossnmin ||
		!out->orgnetnmin_layer || !out->rhiznetnmin_layer || !out->bulknetnmin_layer || !out->mineralnetnmin_layer ||
		!out->rootexudatec || !out->mycnplant || !out->funfixn || !out->littern || !out->netnitrif || !out->nratio ||
		!out->foln || !out->litm || !out->litn || !out->rmresp || !out->rgresp || !out->decresp || !out->totalsoilresp || !out->decwresp ||
		!out->orgfastc || !out->orgfastn || !out->orgslowc || !out->orgslown || !out->orgdeadc || !out->orgdeadn ||
		!out->rhizfastc || !out->rhizfastn || !out->rhizslowc || !out->rhizslown || !out->rhizdeadc || !out->rhizdeadn ||
		!out->bulkfastc || !out->bulkfastn || !out->bulkslowc || !out->bulkslown || !out->bulkdeadc || !out->bulkdeadn ||
		!out->orgmicn || !out->rhizmicn || !out->bulkmicn || !out->totalmicn ||
		!out->orgprotn || !out->rhizprotn || !out->bulkprotn || !out->totalprotn ||
		!out->protectionnyr || !out->deprotectionnyr ||
		!out->funndemand || !out->funndemandgap || !out->funpotrootc || !out->funpotwoodc ||
		!out->funnlimitrootc || !out->funnlimitwoodc ||
		!out->funplantnoverflow || !out->funoverflowtonh4 || !out->funoverflowtono3 ||
		!out->nstart_total || !out->nend_total || !out->dn_total ||
		!out->nstart_soilorg || !out->nend_soilorg || !out->dn_soilorg ||
		!out->nstart_mineral || !out->nend_mineral || !out->dn_mineral ||
		!out->nstart_plantstore || !out->nend_plantstore || !out->dn_plantstore ||
		!out->nstart_vegstruct || !out->nend_vegstruct || !out->dn_vegstruct ||
		!out->nstart_deadwood || !out->nend_deadwood || !out->dn_deadwood ||
		!out->ngas || !out->nbalance_resid)  //Linghui 0729
	{
		printf("Unable to allocate memory for out_struct!\n");
		exit(1);

	}

	init_out(nyrs, out); 								// initialize out structure

}

void pnet_model::memset_climate(clim_struct* clim)
{
	//allocate memory for the climate structure

	clim->year = (int*)malloc((clim->length + 1) * sizeof(int)); //use clim->length+1 so that rstep can start from 1
	clim->doy = (int*)malloc((clim->length + 1) * sizeof(int));
	clim->tmax = (double*)malloc((clim->length + 1) * sizeof(double));
	clim->tmin = (double*)malloc((clim->length + 1) * sizeof(double));
	clim->par = (double*)malloc((clim->length + 1) * sizeof(double));
	clim->prec = (double*)malloc((clim->length + 1) * sizeof(double));
	clim->O3 = (double*)malloc((clim->length + 1) * sizeof(double));
	clim->CO2 = (double*)malloc((clim->length + 1) * sizeof(double));
	clim->NH4dep = (double*)malloc((clim->length + 1) * sizeof(double));
	clim->NO3dep = (double*)malloc((clim->length + 1) * sizeof(double));

	if (!clim->year || !clim->doy || !clim->tmax || !clim->tmin || !clim->par || !clim->prec || !clim->O3 || !clim->CO2 || !clim->NH4dep || !clim->NO3dep)
	{
		printf("Unable to allocate memory for clim_struct!\n");
		exit(1);
	}

	if (clim->timestep == 1)  // daily step
	{
		clim->month = (int*)malloc((clim->length + 1) * sizeof(int));
		clim->day = (int*)malloc((clim->length + 1) * sizeof(int));
//		clim->hour = (double*)malloc(sizeof(double));
//		clim->tmean = (double*)malloc(sizeof(double));
//		clim->vpd = (double*)malloc(sizeof(double));
//		clim->tsoil = (double*)malloc(sizeof(double));
		clim->ifleap = (int*)malloc((clim->length + 1) * sizeof(int));

	}

	for (int i = 0; i < clim->length + 1; i++)
	{
		clim->year[i] = 0; //use clim->length+1 so that rstep can start from 1
		clim->doy[i] = 0;
		clim->tmax[i] = 0.0;
		clim->tmin[i] = 0.0;
		clim->par[i] = 0.0;
		clim->prec[i] = 0.0;
		clim->O3[i] = 0.0;
		clim->CO2[i] = 0.0;
		clim->NH4dep[i] = 0.0;
		clim->NO3dep[i] = 0.0;
		clim->month[i] = 0;
		clim->day[i] = 0;
		clim->ifleap[i] = 0;

	}


}


void pnet_model::memfree_climate(clim_struct* clim)
{

	//clim_struct
	free(clim->year);
	free(clim->doy);
	free(clim->tmax);
	free(clim->tmin);
	free(clim->par);
	free(clim->prec);
	free(clim->O3);
	free(clim->CO2);
	free(clim->NH4dep);
	free(clim->NO3dep);

	if (clim->timestep == 1)  // daily step
	{
		clim->month = (int*)malloc((clim->length + 1) * sizeof(int));
		clim->day = (int*)malloc((clim->length + 1) * sizeof(int));
		//		clim->hour = (double*)malloc(sizeof(double));
		//		clim->tmean = (double*)malloc(sizeof(double));
		//		clim->vpd = (double*)malloc(sizeof(double));
		//		clim->tsoil = (double*)malloc(sizeof(double));
		clim->ifleap = (int*)malloc((clim->length + 1) * sizeof(int));

	}


//	free(clim->ifday);
//	free(clim->hour);
//	free(clim->tmean);
//	free(clim->vpd);
//	free(clim->tsoil);


}

void pnet_model::memfree_out(out_struct* out)
{

	//free memory

	//out_struct

	free(out->grosspsn);
	free(out->netpsn);
	free(out->netcbal);
	free(out->vpd);
	free(out->folmass);
	free(out->plantnMo);
	free(out->nppfol);
	free(out->nppwood);
	free(out->npproot);
	free(out->nep);
	free(out->gpp);
	free(out->waterstress);
	free(out->trans);
	free(out->soilwater);
	free(out->psn);
	free(out->drain);
	free(out->prec);
	free(out->evap);
	free(out->et);
	free(out->plantc);
	free(out->budc);
	free(out->woodc);
	free(out->rootc);
	free(out->plantnYr);
	free(out->budn);
	free(out->ndrain);
	free(out->netnmin);
	free(out->grossnmin);
	free(out->nplantuptake);
	free(out->grossnimob);
	free(out->bgnetnmin);
	free(out->orggrossnmin);
	free(out->rhizgrossnmin);
	free(out->bulkgrossnmin);
	free(out->mineralgrossnmin);
	free(out->orgnetnmin_layer);
	free(out->rhiznetnmin_layer);
	free(out->bulknetnmin_layer);
	free(out->mineralnetnmin_layer);
	free(out->rootexudatec);
	free(out->mycnplant);
	free(out->funfixn);
	free(out->littern);
	free(out->netnitrif);
	free(out->nratio);
	free(out->foln);
	free(out->litm);
	free(out->litn);
	free(out->rmresp);
	free(out->rgresp);
	free(out->decresp);
	free(out->totalsoilresp);
	free(out->decwresp); //Linghui 0729
	free(out->orgfastc);
	free(out->orgfastn);
	free(out->orgslowc);
	free(out->orgslown);
	free(out->orgdeadc);
	free(out->orgdeadn);
	free(out->rhizfastc);
	free(out->rhizfastn);
	free(out->rhizslowc);
	free(out->rhizslown);
	free(out->rhizdeadc);
	free(out->rhizdeadn);
	free(out->bulkfastc);
	free(out->bulkfastn);
	free(out->bulkslowc);
	free(out->bulkslown);
	free(out->bulkdeadc);
	free(out->bulkdeadn);
	free(out->orgmicn);
	free(out->rhizmicn);
	free(out->bulkmicn);
	free(out->totalmicn);
	free(out->orgprotn);
	free(out->rhizprotn);
	free(out->bulkprotn);
	free(out->totalprotn);
	free(out->protectionnyr);
	free(out->deprotectionnyr);
	free(out->funndemand);
	free(out->funndemandgap);
	free(out->funpotrootc);
	free(out->funpotwoodc);
	free(out->funnlimitrootc);
	free(out->funnlimitwoodc);
	free(out->funplantnoverflow);
	free(out->funoverflowtonh4);
	free(out->funoverflowtono3);
	free(out->nstart_total);
	free(out->nend_total);
	free(out->dn_total);
	free(out->nstart_soilorg);
	free(out->nend_soilorg);
	free(out->dn_soilorg);
	free(out->nstart_mineral);
	free(out->nend_mineral);
	free(out->dn_mineral);
	free(out->nstart_plantstore);
	free(out->nend_plantstore);
	free(out->dn_plantstore);
	free(out->nstart_vegstruct);
	free(out->nend_vegstruct);
	free(out->dn_vegstruct);
	free(out->nstart_deadwood);
	free(out->nend_deadwood);
	free(out->dn_deadwood);
	free(out->ngas);
	free(out->nbalance_resid);
	free(out->folm);
	free(out->deadwoodm);
	free(out->woodm);
	free(out->rootm);
	free(out->hom);
	free(out->hon);
	free(out->ndep);
	free(out->fn2o);
	free(out->fno);
	free(out->fn2);
	free(out->fn2ode);
	free(out->fnode);
	free(out->v1);
	free(out->v2);
	free(out->v3);




}


void pnet_model::memfree_site(site_struct* site)
{

	//free memory

	//site_struct

	free(site->distintensity);
	free(site->distyear);
	free(site->distremove);
	free(site->distsoilloss);
	free(site->folregen);
	free(site->distRootloss);
//	free(site->distdoy);


}

void pnet_model::pnet_memfree(site_struct* site, clim_struct* clim, share_struct* share, out_struct* out)
{

	//free memory

	//site_struct
//	memfree_site(site);

	//climate_structs
	memfree_climate(clim);

	//out_struct
	memfree_out(out);
}


void pnet_model::memfree_all(site_struct* site, clim_struct* clim, share_struct* share, out_struct* out)
{

	//free memory

	//site_struct
//	memfree_site(site);

	//climate_structs
	memfree_climate(clim);

	//out_struct
	memfree_out(out);
}
