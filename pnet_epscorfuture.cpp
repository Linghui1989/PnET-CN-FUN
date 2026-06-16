#include "pnet_model.h"
#include <iostream>
#include <string>
#include "netcdf.h"
#include <netcdf>
#include <stdio.h>



using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;



void pnet_model::epscor_region_merrimack_future()
{

	/*****************************************************************************************************************
	 NOTE:

	1, use keyword "DEFINE" to find variables that users have to set for different simulations
	2,


	***************************************************************************************************************/
	char inputname[400];
	char note[200];	
	char climtmax[500],climtmin[500],climpar[500],climppt[500];  //input climate file location
	char pathCO2[500], pathO3[500], pathNH4[500], pathNO3[500];	 //path of air chemistries
	int CN_Mode = 1;											 //if running PnET-CN
	int ystep = 1;												 //unlike matlab code, ystep is initialized here instead of in initvars.c
	char ndep[200];
	FILE* fileoutD;		// output daily result
	FILE * ndephb;												 // hb n deposition
	int rstep;
	int ncycle,ryrs;											 //ncycle, number of climate record cycles for spin up; ryrs, number of record years
	
	veg_struct veg;				     //structure to hold veg data
	site_struct site;				 //structure to hold the site data
	clim_struct climday,climhis;//structure to hold the input climate data
	out_struct out;					 //structure to hold the output data,outrgn for deciduous, outrgncnf for conifer
	share_struct share;				 //structure to hold share data

	int dayobs; // days of observations
	int ndays;
	int doyEnd;  // end of year for doy, leap year:366, otherwise 365
	int yr;
	int nv=0;  // number of variables
	int rec, nsteps;

	


//============DEFINE parameters===================
	int ctl_clm,ctl_lu,ctl_out;
	int ctl_co2;

	ctl_clm = 2;	// comtemporary-Merra
	ctl_lu = 0;		// Backyardm,bkam
	ctl_out = 1;
	ctl_co2 = ctl_clm;


	char lu_dir[450];
	sprintf(lu_dir, "%s%s%s%s", exePath, sep, "organized_data", sep); //directory of regional data
	
	char pnetresult[400];
	sprintf(pnetresult, "%s%s%s%s", exePath, sep, "result", sep);

	char fileout_yr[800];
	
	char climdayfile[450];

	sprintf(climdayfile, "%s%s%s", exePath, PathInput, "climateday_preindustrial.txt");

	


	int yrstart, yrend;		// climate starting year and ending year
	int yrendhis,yrspin;    // end of historical data

	if (ctl_clm ==1 ||ctl_clm ==0 || ctl_clm==3)  // GFDL -OR- Huber
	{
		yrstart=1980;//1979;   		// DEFINE by users, starting year of climate
		yrendhis = 2015;// 2014
		yrend = 2015;//2099;		//2014	// DEFINE by users, ending year of climate
		ncycle = 20;//10;			// DEFINE, number of cycling for spin up simulation
	}

	if (ctl_clm ==2  )  // MERRA
	{
		yrstart = 1000;//1979; 		// DEFINE by users, starting year of climate (simulation)
		yrendhis = 1900;// 2014     // End of spin up year
		yrend = 2050;//2014	        // DEFINE by users, ending year of climate
		ncycle = 56;//10;			// DEFINE, number of cycling for spin up simulation
	}

	printf("yrstart %d, yrendhis %d, yrend %d \n", yrstart, yrendhis, yrend);

	yrspin = yrendhis - yrstart +1;   // years of spinup
	ryrs = yrend - yrendhis;        // years of records


	//calculate historical days
	int dayhis = 0; //the length of spinup days

	for (yr=yrstart;yr<=yrendhis;yr++)
	{
		if ((yr % 4 == 0 && yr % 100 != 0) || ( yr % 400 == 0))
		{
			doyEnd = 366;// leap year; doyEnd = 366;
		} //leap year
		else doyEnd = 365;
		dayhis += doyEnd;
	}

	climhis.timestep = 1;								// PnET historical time step
	climhis.length = dayhis;							// for daily running,12784
//	printf("Days of historical data %d \n", dayhis);

	memset_climate(&climhis);							//allocate memory for the climate structure

//	climspin.timestep = 1;								// PnET daily
//	climspin.length = dayhis * ncycle;					// for daily running,12784
//	memset_climate(&climspin);							//allocate memory for the climate structure


//calculate days of simulation
	dayobs = 0; 

	for (yr=yrendhis+1;yr<=yrend;yr++)
	{
		if ((yr % 4 == 0 && yr % 100 != 0) || ( yr % 400 == 0))
		{
			doyEnd = 366;// leap year; doyEnd = 366;
		} //leap year
		else doyEnd = 365;
		dayobs += doyEnd;
	}
//	printf("Days of daily simulation %d \n", dayobs);
	

	climday.timestep = 1;				// PnET daily
	climday.length = ryrs*365;			// for daily running,12784
	memset_climate(&climday);			//allocate memory for the climate structure

	
	printf("set the memory for climate and output \n");


	
	//================set spatial data name and path==========================
		
		if (ctl_clm ==2)
		{
			ctl_co2 = 1;  // b1
//			sprintf(climtmax, "%s%s%s%s", Path, "agg_macav2metdata_tasmax_", scenario, "_2006_2099_CONUS_daily.nc");
			sprintf(climtmax, "%s%s", lu_dir, "tmax_daily.nc");
			sprintf(climtmin, "%s%s", lu_dir, "tmin_daily.nc");
			sprintf(climpar, "%s%s", lu_dir, "par_daily.nc");
			sprintf(climppt, "%s%s", lu_dir, "ppt_daily.nc");
			sprintf(pathCO2, "%s%s", lu_dir, "CO2_daily.nc");
			sprintf(pathO3, "%s%s", lu_dir, "O3_daily.nc");
			sprintf(pathNO3, "%s%s", lu_dir, "NO3_daily.nc");
			sprintf(pathNH4, "%s%s", lu_dir, "NH4_daily.nc");
		}
			
		
		//================     read landcover data    ==========================
		//define the scope of netCDF data
		static const int NLAT = 122;
		static const int NLON = 127;
		static const int NREC = 54750;
		static const int NC_ERR = 2;
		//create NcVar variable
		NcVar Var; //ncvariable to hold regional data
		NcVar latVar, lonVar;
		//create vector to hold the file
		vector<size_t> startp, countp;
		startp.push_back(0);
		startp.push_back(0);
		startp.push_back(0);
		countp.push_back(1);
		countp.push_back(NLAT);
		countp.push_back(NLON);
		//create array to hold landcover and temporate variable
		double nlcdtemp[NLAT][NLON], nlcd[NLAT][NLON], lats[NLAT];
		//define the landcover path

		sprintf(inputname, "%s%s%s", exePath,PathInput, "landcover_final.nc");

		//create dataFile to read netCDF data
		NcFile dataFile(inputname, NcFile::read);
		Var = dataFile.getVar("var");
		latVar = dataFile.getVar("lat");
		latVar.getVar(lats);
		Var.getVar(startp, countp, nlcd);
		dataFile.close();

		printf("read land cover file \n");
		


		ReadClimDay(&climhis, climdayfile); //read historical data
		printf("Start model run ... \n");

	 

		// read HB N deposition
		sprintf(inputname, "%s%s%s", exePath, PathInput, "date.txt");
		ndephb = fopen(inputname, "r");
		fgets(note, 200, ndephb);

		for (rec = 1; rec <= dayobs; rec++)
		{
			fscanf(ndephb, "%d %d \n", &climday.year[rec], &climday.doy[rec]);
			if (is_leapyear(climday.year[rec]))
				climday.ifleap[rec] = 1;
			else climday.ifleap[rec] = 0;

		}
		fclose(ndephb);


		// ============disturbance=============================================
		Disturbance(&site, &veg, &share, 0);  // to allocate memory for disturbances
		Disturbance(&site, &veg, &share, 2);  // default disturbance as HB



		int ptlat, ptlon;   // use for specific sites
		int latzero, latend, lonzero, lonend;  // customize the range of locations
		int lat, lon;


		// whole domain
		latzero = 0;
		latend = NLAT;
		printf("NLAT is: %d \n", NLAT);
		lonzero = 0;
		lonend = NLON;
		printf("NLON is: %d \n", NLON);

	
	
	
	// specific location plot run
	//Harvard Forest
//		ptlat = 59; // HB
//		ptlon = 48;  // HB32

		ptlat = 42; // HB
		ptlon = 94;  // HB32

	// customization, DEFINE
		latzero =ptlat;
		latend =ptlat+1;
		lonzero =ptlon;
		lonend =ptlon+1;


		// ============initialize for historical running========================
		sprintf(inputname, "%s%s%s", exePath, PathInput, "input.txt");
		//printf(inputname, "%s%s%s", exePath, PathInput, "input-hb-hw.txt");// read hardwoods parameters
		ReadInput(&site, &veg, &share, inputname);									// read in input data from input files
		initvars(&site, &veg, &share);

		ryrs = climday.year[climday.length] - climday.year[1] + 1;//Initialize variables
		memset_out(ryrs, &out);	//allocate memory for the output
		




	for (lat = latzero; lat < latend; lat++)  // HB
	{
		for (lon = lonzero; lon < lonend; lon++) // HB
		{

			double nl = nlcd[lat][lon];

			if (!isnan(nl))
			{
				
				printf("executing lat is %d lon is %d \n", lat, lon);

				// read regional tmax data
				dataFile.open(climtmax, NcFile::read);
				Var = dataFile.getVar("tmax");

				for (rec = 0; rec < NREC; rec++)
				{
					// Read the data one record at a time.
					startp[0] = rec;
					Var.getVar(startp, countp, nlcdtemp);
					climday.tmax[rec + 1] = nlcdtemp[lat][lon];
				} // next record
				dataFile.close();

				// read regional tmin data
				dataFile.open(climtmin, NcFile::read);
				Var = dataFile.getVar("tmin");

				for (rec = 0; rec < NREC; rec++)
				{
					// Read the data one record at a time.
					startp[0] = rec;
					Var.getVar(startp, countp, nlcdtemp);
					climday.tmin[rec + 1] = nlcdtemp[lat][lon];
				} // next record
				dataFile.close();

				// read regional par data
				dataFile.open(climpar, NcFile::read);
				Var = dataFile.getVar("par");

				for (rec = 0; rec < NREC; rec++)
				{
					// Read the data one record at a time.
					startp[0] = rec;
					Var.getVar(startp, countp, nlcdtemp);
					climday.par[rec + 1] = nlcdtemp[lat][lon];
				} // next record
				dataFile.close();

				// read regional ppt data
				dataFile.open(climppt, NcFile::read);
				Var = dataFile.getVar("ppt");

				for (rec = 0; rec < NREC; rec++)
				{
					// Read the data one record at a time.
					startp[0] = rec;
					Var.getVar(startp, countp, nlcdtemp);
					climday.prec[rec + 1] = nlcdtemp[lat][lon];
				} // next record
				dataFile.close();

				// read regional CO2 data
				dataFile.open(pathCO2, NcFile::read);
				Var = dataFile.getVar("CO2");

				for (rec = 0; rec < NREC; rec++)
				{
					// Read the data one record at a time.
					startp[0] = rec;
					Var.getVar(startp, countp, nlcdtemp);
					climday.CO2[rec + 1] = nlcdtemp[lat][lon];
				} // next record
				dataFile.close();

				// read regional O3 data
				dataFile.open(pathO3, NcFile::read);
				Var = dataFile.getVar("O3");

				for (rec = 0; rec < NREC; rec++)
				{
					// Read the data one record at a time.
					startp[0] = rec;
					Var.getVar(startp, countp, nlcdtemp);
					climday.O3[rec + 1] = nlcdtemp[lat][lon];
				} // next record
				dataFile.close();

				// read regional NH4 data
				dataFile.open(pathNH4, NcFile::read);
				Var = dataFile.getVar("NH4");

				for (rec = 0; rec < NREC; rec++)
				{
					// Read the data one record at a time.
					startp[0] = rec;
					Var.getVar(startp, countp, nlcdtemp);
					climday.NH4dep[rec + 1] = nlcdtemp[lat][lon];
				} // next record
				dataFile.close();

				// read regional NO3 data
				dataFile.open(pathNO3, NcFile::read);
				Var = dataFile.getVar("NO3");

				for (rec = 0; rec < NREC; rec++)
				{
					// Read the data one record at a time.
					startp[0] = rec;
					Var.getVar(startp, countp, nlcdtemp);
					climday.NO3dep[rec + 1] = nlcdtemp[lat][lon];
				} // next record
				dataFile.close();

				printf("successful read the regional data \n");

				nsteps = climhis.length;
				site.Lat = lats[lat];  // updating latitude


//				write_clim(&climhis, 1);
				write_clim(&climday, 2);

				// ===================        historical running      ==================
				for (rstep = 1; rstep <= nsteps; rstep++)
				{

					AtmEnviron(&site, &climhis, rstep, &share);
					Phenology(&veg, &climhis, rstep, &share, 1);
					Photosyn(&site, &veg, &climhis, rstep, &share);
					Waterbal(&site, &veg, &climhis, rstep, &share);
					AllocateMo(&veg, &share, rstep, CN_Mode);
					Phenology(&veg, &climhis, rstep, &share, 2);
					CNTrans(&site, &veg, &climhis, rstep, &share);
					Decomp(&site, &veg, &climhis, rstep, &share);
					Leach(&share);
					//storeoutput(&veg, &share, &out, rstep, &ystep, 0);

						//End-of-year activity

					if ((climhis.ifleap[rstep] == 0 && climhis.doy[rstep] == 365) || (climhis.ifleap[rstep] == 1 && climhis.doy[rstep] == 366))
					{
						AllocateYr(&site, &veg, &climhis, rstep + 1, &share, CN_Mode); // Note the rstep+1, not rstep
		//				storeoutput(&veg, &share, &out, rstep+1, &ystep, 1);
						YearInit(&share);
						//				printf("       Executing year %ld %lf \n", climhis.year[rstep], out.woodm[rstep]);
						//				printf("       Executing year %d \n", climhis.year[rstep]);
					}

				}// end of the time loop
				 // write out annual results
		//		WriteoutYr(&climhis, &out, fileout_yr);

				printf("successful run historical date \n");

				// ===================       daily running      ==================
				ndays = climday.length;
				ystep = 1;

				for (rstep = 1; rstep <= ndays; rstep++)
				{
					AtmEnviron(&site, &climday, rstep, &share);
					Phenology(&veg, &climday, rstep, &share, 1);
					Photosyn(&site, &veg, &climday, rstep, &share);
					Waterbal(&site, &veg, &climday, rstep, &share);
					AllocateMo(&veg, &share, rstep, CN_Mode);
					Phenology(&veg, &climday, rstep, &share, 2);
					CNTrans(&site, &veg, &climday, rstep, &share);
					Decomp(&site, &veg, &climday, rstep, &share);
					Leach(&share);
					//storeoutput(&veg, &share, &out, rstep, &ystep, 0);
					//WriteoutMo(&site,&climday, &share,rstep,fileoutM,fileout_step_dec);
	//				printf("       Executing year %d doy %d woodm \n", climday.year[rstep], climday.doy[rstep]);

					if (climday.year[rstep] >= site.nYearStart&&climday.year[rstep]<=site.nYearEnd)
						WriteoutDay(&site, &veg, &climday, &share, rstep, fileoutD);  // it is different from the monthly version

					if (climday.doy[rstep] >= 365)
					{
						if (is_leapyear(climday.year[rstep]))
						{
							doyEnd = 366;
							doyEnd = 365;  // all years have 365 days
						} //leap year
						else doyEnd = 365;
						if (climday.doy[rstep] == doyEnd)
						{
							storeoutput(&veg, &share, &out, rstep + 1, &ystep, 1);//Linghui Meng 06302021
							AllocateYr(&site, &veg, &climday, rstep + 1, &share, CN_Mode); // Note the rstep+1, not rstep
							YearInit(&share);
						}
					}
				}

				//define the name of annual result
				sprintf(fileout_yr, "%s%1f%s", pnetresult, nl, "_Output_annual.csv");
				printf("output file name %s \n", fileout_yr);
				WriteoutYr_netcdf(&climday, &out, fileout_yr);

			}
		} // end of loop		
	}

	printf("   ===============================================\n");
	printf("             Model run ends\n");

	Disturbance(&site, &veg, &share, 100);  // to free memory for disturbances

	memfree_out(&out);
	memfree_climate(&climhis);
	memfree_climate(&climday);



	
}


void pnet_model::write_clim(clim_struct* clim, int i)
{

	int j;
	char climname[450];
	FILE* outfile;  // test zzx
	sprintf(climname, "%s%sclim_print_%d.csv", exePath, PathInput, i);
	outfile = fopen(climname, "w");

	fprintf(outfile, " year, month, day, doy, tmax, tmin, par,"
		"prec, O3, CO2, NH4, NO3"
		"\n");

	for (j = 1; j <= clim->length; j++)

	{
		

		fprintf(outfile, "%6d,%6d,%6d,%6d", clim->year[j], clim->month[j], clim->day[j], clim->doy[j]);
		fprintf(outfile, ",%8.3lf,%8.3lf,%8.3lf,%8.3lf,%9.3lf,%9.3lf,%9.6lf,%9.6lf",
			clim->tmax[j], clim->tmin[j], clim->par[j], clim->prec[j], clim->O3[j], clim->CO2[j],
			clim->NH4dep[j], clim->NO3dep[j]);
		//	fprintf(outfile,",%3d",clim->ifleap[j]);
		fprintf(outfile, "\n");
	}


	fclose(outfile);


}