#include "pnet_model.h"

void pnet_model::AtmEnviron(site_struct* site, clim_struct* clim, int rstep, share_struct* share)
{
	//
	// Environmental calculations for the PnET ecosystem model.
	// 

	double LatRad, r, z, decl, z2, h, TA, AC, hr, es, delta, emean;
	double temp;
	double RH_pct;

	if (clim->tmin[rstep] > clim->tmax[rstep]) {
		temp = clim->tmax[rstep];
		clim->tmax[rstep] = clim->tmin[rstep];
		clim->tmin[rstep] = temp;
	}

	share->Tave = (clim->tmin[rstep] + clim->tmax[rstep]) / 2.0;
	share->Tday = (clim->tmax[rstep] + share->Tave) / 2.0;
	share->Tnight = (share->Tave + clim->tmin[rstep]) / 2.0;


	LatRad = site->Lat * (2.0 * PI) / 360.0;
	r = 1 - (0.0167 * cos(0.0172 * (clim->doy[rstep] - 3)));
	z = 0.39785 * sin(4.868961 + 0.017203 * clim->doy[rstep] + 0.033446 * sin(6.224111 + 0.017202 * clim->doy[rstep]));

	if (fabs(z) < 0.7)
	{
		decl = atan(z / (sqrt(1.0 - pow(z,2))));
	}
	else
	{
		decl = PI / 2.0 - atan(sqrt(1 - pow(z,2)) / z);
	}


	if (fabs(LatRad) >= (PI/2.0))
	{
		if (site->Lat < 0)
		{
			LatRad = (-1.0) * (PI/2.0 - 0.01);
		}
		else
		{
			LatRad = (1.0) * (PI/2.0 - 0.01);
		}
	}

	z2 = -tan(decl) * tan(LatRad);

	if (z2 >= 1.0)
	{
		h = 0;
	}
	else if (z2 <= -1.0)
	{
		h = PI;
	}
	else
	{

		TA = fabs(z2);
		if (TA < 0.7)
		{

			AC = 1.570796 - atan(TA / sqrt(1.0 - pow(TA,2)));
		}
		else
		{

			AC = atan(sqrt(1 - pow(TA,2)) / TA);
		}
		if (z2 < 0)
		{
			h = PI-AC;
		}
		else
		{
			h = AC;
		}
	}

	hr = 2.0 * (h*24.0) / (2.0*PI);  // hour
	share->DayLength = (3600 * hr);  // second
	share->NightLength = (3600 * (24.0 - hr));  //second

	es = 0.61078 * exp(17.26939 * share->Tday / (share->Tday + 237.3));
	delta = 4098.0 * es / ((share->Tday + 237.3)*(share->Tday + 237.3));
	if (share->Tday < 0)
	{
		es = 0.61078 * exp(21.87456 * share->Tday / (share->Tday + 265.5));
		delta = 5808.0 * es / ((share->Tday + 265.5)*(share->Tday + 265.5));
	}
	
	emean = 0.61078 * exp(17.26939 * clim->tmin[rstep] / (clim->tmin[rstep] + 237.3));
	if (clim->tmin[rstep] < 0)
	{
	  emean = 0.61078 * exp(21.87456 * clim->tmin[rstep] / (clim->tmin[rstep] + 265.5));
	}

	share->VPD = es - emean;
	// add RH algorithm
	es = 0.61078 * exp(17.26939 * share->Tave / (share->Tave + 237.3)); //recalculate es with Tave for Tw
	if (share->Tave < 0)
	{
		es = 0.61078 * exp(21.87456 * share->Tave / (share->Tave + 265.5));
	}


	share->RH = (emean / es) * 100.0;
	if (share->RH < 0.1) share->RH = 0.1;
	if (share->RH > 100.0) share->RH = 100.0;

	const double RH = 85.0;
	
	//Stull 2011, Wet-bulb temperature from relative humidity and air temperature
	share->Tw = share->Tave * std::atan(0.151977 * std::sqrt(RH + 8.313659))
		          + std::atan(share->Tave + RH) - std::atan(RH - 1.676331)
		          + 0.00391838 * std::pow(RH, 1.5) * std::atan(0.023101 * RH) - 4.686035;

	if (!std::isfinite(share->Tw)) share->Tw = share->Tave;  // safety

	if (share->Tw > share->Tave) share->Tw = share->Tave;
	if (share->Tw < share->Tave - 6.0) share->Tw = share->Tave - 6.0;

	if (share->Tw < -25.0) share->Tw = -25.0;
	if (share->Tw > 10.0) share->Tw = 10.0;

}

int pnet_model::getdays(int year,int doy)  // get number of days in month including doy
{
	int days[13], i ;

	if(doy>367 || doy<1)
	{
		printf("Can not convert the number doy %d to month\n",doy);
		return 0;
	}


	for(i=1; i<13; i++) days[i]=31;
	days[2]=28;
	if ((year % 4 == 0 && year % 100 != 0) || ( year % 400 == 0))days[2]=29; //Feb
	days[4]=days[6]=days[9]=days[11]=30;//Apr,Jun,Sep,Nov

	i=1;
	while(doy>days[i])
	{
		  doy -=days[i];
		  i++;	   
	}
	
	return days[i];

}



int pnet_model::getmonth(int year,int doy)  // get the month of doy
{
	int days[13], i ;

	if(doy>367 || doy<1)
	{
		printf("Can not convert the number to date");
		return 0;
	}


	for(i=1; i<13; i++) days[i]=31;
	days[2]=28;
	if ((year % 4 == 0 && year % 100 != 0) || ( year % 400 == 0))days[2]=29; //Feb
	days[4]=days[6]=days[9]=days[11]=30;//Apr,Jun,Sep,Nov

	i=1;
	while(doy>days[i])
	{
		  doy -=days[i];
		  i++;
	}

	return i;


}


int pnet_model::is_leapyear(int year)  // determine if it is a leap year
{

	if ((year % 4 == 0 && year % 100 != 0) || ( year % 400 == 0)) return 1;
	else return 0;

}
