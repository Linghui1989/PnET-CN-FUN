#include "pnet_model.h"

void pnet_model::pnet_daily()
{

	char climname[200];
	char inputname[350];

	char climdayfile[350];	// for daily climate input
	char climhrfile[200];	// for hourly climate input to compare tower flux

	int CN_Mode = 1; //if running PnET-CN
	int i, rstep, nyrs, ndays;
	int ystep = 1; //unlike matlab code, ystep is initialized here instead of in initvars.c 
	int monspin;
	int doyEnd;  // end of year for doy, leap year:366, otherwise 365
	int YrStartSpin, YrEndSpin; // years of starting and end for spinning up


	veg_struct veg;				//structure to hold veg data
	site_struct site;			//structure to hold the site data
	clim_struct clim;			//structure to hold the input climate data
	out_struct out;				//structure to hold the output data
	share_struct share;			//structure to hold share data

	clim_struct climmon, climday, climhr, climdayspin;			//structure to hold the input climate data,monthly,daily,hourly

	FILE* fileoutM;		// output monthly result
	FILE* fileoutD;		// output daily result

	sprintf(inputname, "%s%s%s", exePath, PathInput, "input.txt");  // path of input files
	sprintf(climdayfile, "%s%s%s", exePath, PathInput, "climateday.clim");//daily climate



	clim.timestep = 0;		//  monthly climate
	climday.timestep = 1;	// daily climate
	climdayspin.timestep = 1;	// daily climate


	// read in input data from input files

	ReadInput(&site, &veg, &share, inputname);

	ReadClimDay(&climday, climdayfile);	// read daily climate

	//Initialize variables

	initvars(&site, &veg, &share);


	YrStartSpin = 1001;		//Linghui MENG 20200508
	YrEndSpin = 2023;


	share.yrspin = YrEndSpin - YrStartSpin + 1;		//allocate memory for share file
	climday.length = share.yrspin * 365;             //365 days for each year


	//allocate memory
	memset_out(share.yrspin, &out);

	share.ifO3EffectOnPSN = 1;  // O3 effect
//	share.ifO3EffectOnPSN= 0;  // O3 effect

	share.kO3EffScalar = 1.12;
	share.kO3Eff = share.kO3Eff * share.kO3EffScalar; // 1.12;daily version
	share.WUEO3Eff = 0;

	// ============disturbance=============================================
	Disturbance(&site, &veg, &share, 0);  // to allocate memory for disturbances
	Disturbance(&site, &veg, &share, 1);	// default disturbance

	FILE* fileobs;		// input daily observed gpp
	char inputobs[200];

	printf("   ===============================================\n\n");
	printf("             PnET model starts daily run\n");

	ndays = climday.length;
	share.ISA = 0.00;


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
		Leach(&climday,rstep ,&share);
		//		storeoutput(&veg, &share, &out, rstep, &ystep, 0);
		if (climday.year[rstep] >= site.nYearStart && climday.year[rstep] <= site.nYearEnd)
			WriteoutDay(&site, &veg, &climday, &share, rstep, fileoutD);  // it is different from the monthly version

		//End-of-year activity

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
				AllocateYr(&site, &veg, &climday, rstep + 1, &share, CN_Mode); // Note the rstep+1, not rstep
				storeoutput(&veg, &share, &out, rstep + 1, &ystep, 1);//Linghui Meng 06302021
				YearInit(&share);
			}
		}

		printf("       Executing year %d  doy  %d  \n", climday.year[rstep], climday.doy[rstep]);


	}// end of the run loop


	// write out annual results
	WriteoutYr(&climday, &out);

	//free memory and close files
	memfree_all(&site, &clim, &share, &out);
	Disturbance(&site, &veg, &share, 100);	// default disturbance

	printf("   ===============================================\n");

}






void pnet_model::ReadClimDay(clim_struct* clim, char* climname)
{
	char buf[200];
	int i;

	FILE* fileClim;

	if ((fileClim = fopen(climname, "rt")) == NULL)
	{	
		printf("Unable to open the clim file!\n");
		exit(1) ;
	}

	i=0;
	while(NULL!=fgets(buf,200,fileClim))i++; 

	clim->length=i-1;

	rewind(fileClim);
	

	//allocate memory for the climate structure  // zzx rearrange the code
	memset_climate(clim);


	fgets(buf, 200, fileClim);

	for (i = 1; i <= clim->length; i++)
	{
		fscanf(fileClim, "%d", &clim->year[i]); //year
		fscanf(fileClim, "%d", &clim->month[i]); //mon
		fscanf(fileClim, "%d", &clim->day[i]); //day
		fscanf(fileClim, "%d", &clim->doy[i]); //doy
		fscanf(fileClim, "%lf",&clim->tmax[i]); //tmax
		fscanf(fileClim, "%lf",&clim->tmin[i]); //tmin
		fscanf(fileClim, "%lf",&clim->par[i]); //par
		fscanf(fileClim, "%lf",&clim->prec[i]); //prec
		fscanf(fileClim, "%lf",&clim->O3[i]); //O3
		fscanf(fileClim, "%lf",&clim->CO2[i]); //CO2
		fscanf(fileClim, "%lf",&clim->NH4dep[i]); //NH4dep
		fscanf(fileClim, "%lf",&clim->NO3dep[i]); //NO3dep
//		fscanf(fileClim, "%lf",&clim->tsoil[i]); //NO3dep

	}
	fclose(fileClim);
}



void pnet_model::ReadClimHour(clim_struct* clim, char* climname)
{
	char buf[200];
	int i;

	FILE* fileClim;

	if ((fileClim = fopen(climname, "rt")) == NULL)
	{	
		printf("Unable to open the clim file!\n");
		exit(1) ;
	}

	i=0;
	while(NULL!=fgets(buf,200,fileClim))i++; 

	clim->length=i-1;

	rewind(fileClim);
	

	//allocate memory for the climate structure  // zzx rearrange the code
	clim->year = (int*)malloc((clim->length+1)*sizeof(int)); //use clim->length+1 so that rstep can start from 1
	clim->month = (int*)malloc((clim->length+1)*sizeof(int)); //use clim->length+1 so that rstep can start from 1
	clim->day = (int*)malloc((clim->length+1)*sizeof(int)); //use clim->length+1 so that rstep can start from 1
	clim->hour = (double*)malloc((clim->length+1)*sizeof(double));

	clim->doy = (int*)malloc((clim->length+1)*sizeof(int));
	clim->ifday = (int*)malloc((clim->length+1)*sizeof(int));

	clim->tmax = (double*)malloc((clim->length+1)*sizeof(double));
	clim->tmin = (double*)malloc((clim->length+1)*sizeof(double));
	clim->tmean = (double*)malloc((clim->length+1)*sizeof(double));
	clim->vpd = (double*)malloc((clim->length+1)*sizeof(double));

	clim->par = (double*)malloc((clim->length+1)*sizeof(double));
	clim->prec = (double*)malloc((clim->length+1)*sizeof(double));
	clim->O3 = (double*)malloc((clim->length+1)*sizeof(double));
	clim->CO2 = (double*)malloc((clim->length+1)*sizeof(double));
	clim->NH4dep = (double*)malloc((clim->length+1)*sizeof(double));
	clim->NO3dep = (double*)malloc((clim->length+1)*sizeof(double));
	if (!clim->year || !clim->doy || !clim->tmax || !clim->tmin || !clim->par || !clim->prec || !clim->O3 || !clim->CO2 || !clim->NH4dep || !clim->NO3dep)
	{
		printf("Unable to allocate memory for clim_struct!\n");
		exit(1) ;
	}


	fgets(buf, 200, fileClim);

	for (i = 1; i <= clim->length; i++)
	{
		fscanf(fileClim, "%d",&clim->year[i]); //year
		fscanf(fileClim, "%d", &clim->month[i]); //doy
		fscanf(fileClim, "%d", &clim->day[i]); //doy
		fscanf(fileClim, "%lf", &clim->hour[i]); //doy
		fscanf(fileClim, "%d", &clim->doy[i]); //doy
		fscanf(fileClim, "%d",&clim->ifday[i]); //doy
		fscanf(fileClim, "%lf",&clim->tmean[i]); //tmax
		fscanf(fileClim, "%lf",&clim->vpd[i]); //tmin
		fscanf(fileClim, "%lf",&clim->par[i]); //par
		fscanf(fileClim, "%lf",&clim->prec[i]); //prec
		fscanf(fileClim, "%lf",&clim->O3[i]); //O3
		fscanf(fileClim, "%lf",&clim->CO2[i]); //CO2
		fscanf(fileClim, "%lf",&clim->NH4dep[i]); //NH4dep
		fscanf(fileClim, "%lf\n",&clim->NO3dep[i]); //NO3dep


	}
	fclose(fileClim);
}

void pnet_model::preprocess_climate(site_struct* site,share_struct* share,clim_struct* climday,clim_struct* clim)
{

	int i,j,nyrs,n;
	double num[13];

	clim_struct *climmon;

	climmon = (clim_struct*)malloc(sizeof(clim_struct));

//	climmon->timestep = 0;


//	climmon->length = 13;
	//allocate memory for the  average monthly climate
	memset_climate(climmon);

	nyrs = climday->year[climday->length] - climday->year[1] + 1;

	// set initial values zero
	for ( i=1; i<13;i++)
	{
		num[i] = 0.0;  // number of data in each month

		climmon->year[i] = 1; //year
		climmon->doy[i] = i*30-15; //doy

		climmon->tmax[i] = 0.0;
		climmon->tmin[i] = 0.0;
		climmon->par[i] = 0.0;
		climmon->prec[i] = 0.0;
		climmon->O3[i] = 0.0;
		climmon->CO2[i] = 0.0;
		climmon->NH4dep[i] = 0.0;
		climmon->NO3dep[i] = 0.0;

	}

	j = 1;   // sequence of climday

	for (j=1;j < climday->length+1;j++)
	{

		for ( i=1; i<13;i++) // month of year loop
		{

			if (climday->month[j] == i)
			{
				num[i] += 1;
				climmon->tmax[i] += climday->tmax[j];
				climmon->tmin[i] += climday->tmin[j];
				climmon->par[i] += climday->par[j];
				climmon->prec[i] += climday->prec[j];
				climmon->O3[i] += climday->O3[j];
				climmon->CO2[i] += climday->CO2[j];
				climmon->NH4dep[i] += climday->NH4dep[j];
				climmon->NO3dep[i] += climday->NO3dep[j];

			}

		}

	}

	for ( i=1; i<13;i++)
	{
		climmon->tmax[i] /= num[i];
		climmon->tmin[i] /= num[i];
		climmon->par[i] /= num[i] ;
		climmon->prec[i] /= nyrs;
		climmon->O3[i] /= num[i];
		climmon->CO2[i] /= num[i];
		climmon->NH4dep[i] /= nyrs;
		climmon->NO3dep[i] /= nyrs;

	}

	for (j=1;j < share->yrspin+1;j++)
	{
		for ( i=1; i<13;i++) // month of year loop
		{
			n=(j-1)*12+i;

			clim->year[n] = j; //year
			clim->doy[n] = climmon->doy[i]; //doy

			clim->tmax[n] = climmon->tmax[i];
			clim->tmin[n] = climmon->tmin[i];
			clim->par[n] = climmon->par[i];
			clim->prec[n] = climmon->prec[i];
			clim->O3[n] = climmon->O3[i];

			clim->CO2[n] = 282.23 + exp(clim->year[n]/51.35)*1.03*1.0e-15; //(Franks,2013, New Phytologist, 197:1077-1094)

			clim->NH4dep[n] = climmon->NH4dep[i]/10;  // 0.1 of the observed value
			clim->NO3dep[n] = climmon->NO3dep[i]/10;

		}

	}

	memfree_climate(climmon);




}


/*
void pnet_model::clim_day_mon(site_struct* site,share_struct* share,clim_struct* climday,clim_struct* clim)
{

	int i,j,jj;
	int doy[13]={0,15,46,74,105,135,166,196,227,258,288,319,349};
	double m_tmax,m_tmin,m_preci,m_par,m_co2,m_nh4,m_no3,m_o3;
	int month;

	j=clim->length;
	jj=0;
	month = 12 ;
	m_tmax =0;
	m_tmin =0;
	m_preci =0.0;
	m_par = 0.0;
	m_co2 = 0.0;
	m_nh4=0.0;
	m_no3 = 0.0;
	m_o3 = 0.0;

	for (i = climday->length; i >0; i--)
	{
		if (climday->month[i] != month || i == 1 )
		{

			if ( i==1)
			{
				jj++;
				m_tmax += climday->tmax[i];
				m_tmin += climday->tmin[i];
				m_par += climday->par[i];
				m_preci += climday->prec[i];
				m_o3 += climday->O3[i];
				m_co2 += climday->CO2[i];
				m_nh4 += climday->NH4dep[i];
				m_no3 += climday->NO3dep[i];

			}


			m_tmax = m_tmax/jj;
			m_tmin = m_tmin/jj;
			m_par = m_par/jj;
			m_o3 = m_o3/jj;
			m_co2 = m_co2/jj;


			clim->year[j]=climday->year[i+1];
			clim->doy[j]=doy[month];
			clim->tmin [j]= m_tmin;
			clim->tmax [j]= m_tmax;
			clim->par[j] = m_par;
			clim->prec [j]= m_preci; //

			clim->O3[j] = m_o3;
			clim->CO2[j] = m_co2;
			clim->NH4dep[j] = m_nh4;
			clim->NO3dep[j] = m_no3;

			m_tmax =0;
			m_tmin =0;
			m_preci =0.0;
			m_par = 0.0;

			m_co2 = 0.0;
			m_nh4=0.0;
			m_no3 = 0.0;
			m_o3 = 0.0;

			month = climday->month[i];  // previous month

			j--;	// next monthly record

			m_tmax += climday->tmax[i];
			m_tmin += climday->tmin[i];
			m_par += climday->par[i];
			m_preci += climday->prec[i];
			m_o3 += climday->O3[i];
			m_co2 += climday->CO2[i];
			m_nh4 += climday->NH4dep[i];
			m_no3 += climday->NO3dep[i];

			jj=1;  // reinitialization


		}

		else
		{
			m_tmax += climday->tmax[i];
			m_tmin += climday->tmin[i];
			m_par += climday->par[i];
			m_preci += climday->prec[i];
			m_o3 += climday->O3[i];
			m_co2 += climday->CO2[i];
			m_nh4 += climday->NH4dep[i];
			m_no3 += climday->NO3dep[i];
			jj++;
		}


	}  // end of record loop


	//fill clim
	month = (climday->year[climday->length]-climday->year[1]+1) * 12;  // number of month with observations
	for (j=clim->length - month; j>0; j--)
	{
		clim->year[j]=clim->year[j+12]-1;
		clim->doy[j]=clim->doy[j+month];
		clim->tmin [j]= clim->tmin[j+month];
		clim->tmax [j]= clim->tmax[j+month];
		clim->par[j] = clim->par[j+month];
		clim->prec [j]= clim->prec[j+month]; //

		clim->O3[j] = clim->O3[j+month];
	//	clim->CO2[j] = 282.23 + exp(clim->year[j]/51.35)*1.03*1.0e-15; //(Franks,2013, New Phytologist, 197:1077-1094);
	//	clim->NH4dep[j] = clim->NH4dep[j+month]; //Afshin removed commented out
	//	clim->NO3dep[j] = clim->NO3dep[j+month]; //Afshin removed commented out
	//	clim->CO2[j]=CO2_Daily(clim->year[j], clim->doy[j], 0);
		
	//	if (clim->year[j]<1860)clim->NH4dep[j] = 0.00928; //Afshin commented out
	//	else clim->NH4dep[j] = (4.1176*clim->year[j]- 7508.8)/1000/12*0.742;  //Afshin commented out
	//	clim->NO3dep[j] = clim->NH4dep[j] * 0.34771;  //Afshin commented out


	}


}
*/

void pnet_model::clim_day_fill(site_struct* site,share_struct* share,clim_struct* climday,clim_struct* clim)
{
	// using climday to fill clim
	// climday: input
	// clim: output


	int i,j;

	//fill clim
	i = clim->length - climday->length;

	for (j=i+1; j<=clim->length; j++)
	{
		clim->year[j]=climday->year[j-i];
		clim->doy[j]=climday->doy[j-i];
		clim->month[j]=climday->month[j-i];
		clim->tmin [j]= climday->tmin[j-i];
		clim->tmax [j]= climday->tmax[j-i];
		clim->par[j] = climday->par[j-i];
		clim->prec [j]= climday->prec[j-i]; //
		clim->O3[j] = climday->O3[j-i];
		clim->CO2[j] = climday->CO2[j-i];
		clim->NH4dep[j] = climday->NH4dep[j-i];
		clim->NO3dep[j] = climday->NO3dep[j-i];

	//	clim->CO2[j] = 282.23 + exp(clim->year[j]/51.35)*1.03*1.0e-15; //(Franks,2013, New Phytologist, 197:1077-1094);
	//	clim->NH4dep[j] = clim->NH4dep[j+month];
	//	clim->NO3dep[j] = clim->NO3dep[j+month];
	//	if (clim->year[j]<1860)clim->NH4dep[j] = 0.00928;
	//	else clim->NH4dep[j] = (4.1176*clim->year[j]- 7508.8)/1000/12*0.742;
	//	clim->NO3dep[j] = clim->NH4dep[j] * 0.34771;

		clim->ifleap[j] = climday->ifleap[j-i];


	}

	for (j=i; j>0; j--)
	{
		clim->year[j]=clim->year[j+365]-1;
		clim->doy[j]=clim->doy[j+climday->length];
		clim->month[j]=clim->month[j+climday->length];
		clim->tmin [j]= clim->tmin[j+climday->length];
		clim->tmax [j]= clim->tmax[j+climday->length];
		clim->par[j] = clim->par[j+climday->length];
		clim->prec [j]= clim->prec[j+climday->length]; //

		clim->ifleap[j]= clim->ifleap[j+climday->length]; //

		if (clim->doy[j+1]==1) // first day of following year
		{
			clim->year[j]=clim->year[j+1]-1;
		}
		else
		{
			clim->year[j]=clim->year[j+1];

		}

		clim->O3[j] = clim->O3[j+climday->length];
	/*	clim->CO2[j] = 282.23 + exp(clim->year[j]/51.35)*1.03*1.0e-15; //(Franks,2013, New Phytologist, 197:1077-1094);

		if (clim->year[j]<1860)clim->NH4dep[j] = 0.00928/30.0;
		else clim->NH4dep[j] = (4.1176*clim->year[j]- 7508.8)/1000/12*0.742/30.0;
		clim->NO3dep[j] = clim->NH4dep[j] * 0.34771;
		*/ //Afshin commented out
		clim->CO2[j] = clim->CO2[j+climday->length];
		clim->NO3dep[j] = clim->NO3dep[j+climday->length];
		clim->NH4dep[j] = clim->NH4dep[j+climday->length];


	}


}


