/*
 * ccp.c
 *
 * Copyright 2021 Joan Antoni Parera Portell
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <sacio.h>
#include <fftw3.h>

#define MAX 5000
#define PI 3.14159265359
#define DEG2RAD 0.017453292520
#define RAD2DEG 57.295779513
#define DEG2KM 111.195
#define KM2DEG 0.00899321
#define EARTHRAD 6371.0

/* Funció per calcular la distància sobre el cercle màxim entre dos
 * punts, així com l'azimut */
void garc(float lat0, float lon0, float lat1, float lon1, float *dist, float *az){
	float dlon, angle, az_angle, x, y;

	lat0 *= DEG2RAD; lon0 *= DEG2RAD; lat1 *= DEG2RAD; lon1 *= DEG2RAD;
	dlon = lon1-lon0;

	/* Distància */
	angle = acos(sin(lat0)*sin(lat1)+cos(lat0)*cos(lat1)*cos(dlon));
	*dist = angle*RAD2DEG*DEG2KM;

	/* Azimut */
	x = sin(dlon)*cos(lat1);
	y = cos(lat0)*sin(lat1)-sin(lat0)*cos(lat1)*cos(dlon);
	az_angle = atan2(x,y);
	*az	= az_angle*RAD2DEG;

	if (*az < 0){
		*az += 360.0;
	}
}

/* Funció per projectar un punt sobre un perfil */
void proj(float dist, float az0, float az1, float *x, float *y){
	float azdiff;

	azdiff = az0-az1;
	azdiff *= DEG2RAD;
	*x = dist*cos(azdiff);
	*y = dist*sin(azdiff);
}

/* Funció principal */
int main(int argc, char **argv){

	int t0, w, n, wsum, nlen, nerr, nlay, ncols, max=MAX, i, j, k, ray;
	float delta, dx, dz, dep, depmin, depmax, beg, p, array[MAX], z,
	inilat, inilon, finlat, finlon, len, azim, stla, stlo, baz, hw,
	az0, az1, dist0, dist1, x, y, xi, yi, ds, dt, a, b, c, up, us, zf, nu,
	acctime, amp, fzr, wl, xamp, dzf, stx, l, q, zv, modn,
	x0, x1, y0, y1, T, freq, den, num;
	char rf[500], outfile[500], model[500], pvar[20];
	FILE *prm_file, *list_file, *mod_file, *out_file;

	/* Arrays de l'FFT */
	fftw_complex fftin[MAX], fftout[MAX], phase;
	fftw_plan plan1, plan2;

	if(argc < 3){
		printf("\nUsage: ccp [par file] [rf list]\n");
		exit(1);
	}

	char *prm = argv[1];
	char *rflist = argv[2];

	/* ---------Comprovar si existeix el fitxer de paràmetres-------- */
	prm_file = fopen(prm, "r");
	if(prm_file == NULL){
		printf("\nParameter file doesn't exist.\n");
		exit(1);
	}

	printf("\n*******************************************************\n");
	printf("******************** CCP STACKING *********************\n");
	printf("*******************************************************\n");
	printf("per Joan A. Parera Portell (2023)\n");
	printf("\nPar. file: %s\n", prm);
	printf("RF file: %s\n", rflist);

	/* ---------------Lectura del fitxer de paràmetres--------------- */
	printf("\n-Parameters-\n");
	w=0;
	wsum=0;
	w=fscanf(prm_file, "%f", &inilat); wsum += w;
	w=fscanf(prm_file, "%f", &inilon); wsum += w;
	w=fscanf(prm_file, "%f", &finlat); wsum += w;
	w=fscanf(prm_file, "%f", &finlon); wsum += w;
	w=fscanf(prm_file, "%d", &t0); wsum += w;
	w=fscanf(prm_file, "%f,%f", &dx, &dz); wsum += w;
	w=fscanf(prm_file, "%f,%f", &depmin, &depmax); wsum += w;
	w=fscanf(prm_file, "%f", &hw); wsum += w;
	w=fscanf(prm_file, "%s", outfile); wsum += w;
	w=fscanf(prm_file, "%s", model); wsum += w;
	w=fscanf(prm_file, "%s", pvar); wsum += w;
	w=fscanf(prm_file, "%f", &zv); wsum += w;
	w=fscanf(prm_file, "%f", &nu); wsum += w;
	w=fscanf(prm_file, "%f", &freq); wsum += w;
	w=fscanf(prm_file, "%d", &ray); wsum += w;
	fclose(prm_file);

	if(wsum != 17){
		printf("Error reading parameter file. Exiting...\n");
		exit(1);
	}

	printf("Initial lat/lon: \t%.2f %.2f deg\n", inilat, inilon);
	printf("Final lat/lon: \t\t%.2f %.2f deg\n", finlat, finlon);
	printf("Initial time: \t\t%d s\n", t0);
	printf("Delta x and delta z: \t%.2f %.2f km\n", dx, dz);
	printf("Min/max depths: \t%.2f %.2f km\n", depmin, depmax);
	printf("Half width: \t\t%.2f km\n", hw);
	printf("Output slice: \t\t%s\n", outfile);
	printf("Earth model: \t\t%s\n", model);
	printf("Ray param. variable: \t%s\n", pvar);
	printf("Depth scaling exp. term:%.2f\n", zv);
	printf("Phase weight exp. term: %.2f\n", nu);
	printf("RF frequency: \t\t%.2f\n", freq);
	printf("FZ (0) or rays (1): \t%d\n", ray);

	list_file = fopen(rflist, "r");
	if(list_file == NULL){
		printf("\nList file doesn't exist.\n");
		exit(1);
	}
	fclose(list_file);

	mod_file = fopen(model, "r");
	if(mod_file == NULL){
		printf("\nModel file doesn't exist.\n");
		exit(1);
	}

	nlay=round(depmax/dz);
	struct earthmodel{
		float dep0[1000], vp0[1000], vs0[1000];
	} emod0;
	struct learthmodel{
		float vp[nlay], vs[nlay];
	} emod;

	w=3;
	n=0;
	while(w==3){
		w=fscanf(mod_file, "%f,%f,%f", &emod0.dep0[n],&emod0.vp0[n],&emod0.vs0[n]);
		n++;
	}
	fclose(mod_file);

	n=0;
	for(j=0; j<nlay; j++){
		z=dz*j;
		if(j==0){
			emod.vp[j]=emod0.vp0[0];
			emod.vs[j]=emod0.vs0[0];
		}
		else{
			while(emod0.dep0[n]<z){n++;}
			if(emod0.dep0[n]==z){
				emod.vp[j]=emod0.vp0[n];
				emod.vs[j]=emod0.vs0[n];
			}
			else{
				x0=emod0.dep0[n-1];
				x1=emod0.dep0[n];
				y0=emod0.vp0[n-1];
				y1=emod0.vp0[n];
				emod.vp[j] = y0 + ((y1-y0)/(x1-x0)) * (z - x0);
				y0=emod0.vs0[n-1];
				y1=emod0.vs0[n];
				emod.vs[j] = y0 + ((y1-y0)/(x1-x0)) * (z - x0);
			}
		}
	}

	/* ------------------Inicialització del perfil------------------- */
	garc(inilat,inilon,finlat,finlon,&len,&azim);
	ncols = round(len/dx);
	printf("\nProfile length: %.2f km\n", len);
	printf("Profile azimuth: %.2f deg\n", azim);
	printf("Profile size: %dx%d\n", ncols, nlay-(int)(depmin/dz));

	float* perfil = calloc(ncols * nlay, sizeof(float));
	float* pnorm  = calloc(ncols * nlay, sizeof(float));
	float* absphase = calloc(ncols * nlay, sizeof(float));
	fftw_complex* pphase = calloc(ncols * nlay, sizeof(fftw_complex));

	/* --------------Start looping through Rfs in list--------------- */
	list_file = fopen(rflist, "r");
	w = 0;
	printf("\nProcessing RFs...");
	while(fgets(rf,500,list_file) != NULL){
		sscanf(rf,"%s",rf);
		rsac1(rf, array, &nlen, &beg, &delta, &max, &nerr, strlen(rf));

		if (nerr != 0){
			printf("\nError reading SAC file: %s\n", rf);
			exit (nerr);
		}

		getfhv(pvar, &p, &nerr, strlen(pvar));
		getfhv("BAZ", &baz, &nerr, strlen("BAZ"));
		getfhv("STLA", &stla, &nerr, strlen("STLA"));
		getfhv("STLO", &stlo, &nerr, strlen("STLO"));

		if(ray==0){
			if(w==0){
				plan1 = fftw_plan_dft_1d(nlen, fftin, fftin, FFTW_FORWARD, FFTW_ESTIMATE);
				plan2 = fftw_plan_dft_1d(nlen, fftin, fftout, FFTW_BACKWARD, FFTW_ESTIMATE);
			}
			for(j=0;j<nlen;j++){
				fftin[j] = array[j];
			}
			fftw_execute(plan1);
			for(j=0;j<nlen;j++){
				if(j>=nlen/2){
					fftin[j] = 0;
				}
			}
			fftw_execute(plan2);
			for(j=0;j<nlen;j++){
				fftout[j] /= (nlen/2);
				if(cabs(fftout[j]) > 0) {
					fftout[j] /= cabs(fftout[j]);
				}
			}
		}
		w++;

		garc(inilat,inilon,stla,stlo,&dist0,&az0);
		proj(dist0,az0,azim,&x,&y);
		stx = x;

		if(beg<0){acctime = -beg;}
		else{acctime = beg;}

		for(j=0; j<nlay; j++){
			up = 1/((EARTHRAD/(EARTHRAD-j*dz))*emod.vp[j]);
			us = 1/((EARTHRAD/(EARTHRAD-j*dz))*emod.vs[j]);
			a = -EARTHRAD*log((EARTHRAD-(j+1)*dz)/EARTHRAD);
			zf = -EARTHRAD*log((EARTHRAD-j*dz)/EARTHRAD);
			dzf = a-zf;
			ds = p*dzf/sqrtf(us*us-p*p);

			a = us*us-p*p;
			b = up*up-p*p;

			if(a>0 && b>0){
				dt = sqrtf(a)*dzf-sqrtf(b)*dzf;
				acctime += dt;
				k = round(acctime/delta);

				if(k >= 0 && k < nlen) {
					amp = array[k];
					proj(ds,baz,azim,&xi,&yi);
					x += xi;
					y += yi;
					if(ray==0){
						phase = fftout[k];
						T = 1/freq;
						wl = T*(1/us);
						fzr = sqrtf(0.5*wl*(dz*j)+0.0625*wl*wl);
						q = powf((dz*j), zv);

						if(fabs(y)<=hw){
							for(n=0;n<ncols;n++){
								l = n*dx-x;
								a = 0.5*1*fzr;
								b = -0.5*(l*l)/(a*a);
								if(b < -20){c=0;}
								else{c = exp(b)/(a*sqrtf(2*PI));}
								xamp = amp*q*c;
								if(xamp != 0){
									perfil[j*ncols+n] += xamp;
									pnorm[j*ncols+n] += 1;
									pphase[j*ncols+n] += phase;
								}
							}
						}
					}
					else{
						n = (int)(x/dx);
						if(n>=0 && x<=len && fabs(y)<=hw){
							perfil[j*ncols+n] += amp;
							pnorm[j*ncols+n] += 1;
						}
					}
				}
			}
		}
	}
	if(ray==0){
		fftw_destroy_plan(plan1);
		fftw_destroy_plan(plan2);
		fftw_cleanup();
	}
	printf("%d!\n", w);
	fclose(list_file);

	/* Aplicació de la fase a les amplituds */
	if(ray==0){
		for(i=0; i<nlay; i++){
			for(w=0; w<ncols; w++){
				if(pnorm[i*ncols+w] > 0){
					perfil[i*ncols+w] /= pnorm[i*ncols+w];
					absphase[i*ncols+w] = powf(cabs(pphase[i*ncols+w])/pnorm[i*ncols+w],nu);
					perfil[i*ncols+w] *= absphase[i*ncols+w];
				}
			}
		}
		for(i=0; i<nlay; i++){
			for(w=0; w<ncols; w++){
				/* Suavitzat per files*/
				if(w>0){
					perfil[i*ncols+w] = 0.5*perfil[i*ncols+w]+(1-0.5)*perfil[i*ncols+w-1];
				}
			}
		}
		for(i=0; i<nlay; i++){
			for(w=ncols-1; w>=0; w--){
				if(w<ncols-1){
					perfil[i*ncols+w] = 0.5*perfil[i*ncols+w]+(1-0.5)*perfil[i*ncols+w+1];
				}
			}
		}
		for(w=0; w<ncols; w++){
			for(i=0; i<nlay; i++){
				if(i>0){
					perfil[i*ncols+w] = 0.5*perfil[i*ncols+w]+(1-0.5)*perfil[(i-1)*ncols+w];
				}
			}
		}
		for(w=0; w<ncols; w++){
			for(i=nlay-1; i>=0; i--){
				if(i<nlay-1){
					perfil[i*ncols+w] = 0.5*perfil[i*ncols+w]+(1-0.5)*perfil[(i+1)*ncols+w];
				}
			}
		}
	}

	/* Escriptura a un fitxer */
	out_file = fopen(outfile, "w");
	fprintf(out_file, "%s,%s,%s,%s,%s,%s,%s\n", "x","z","a","lat","lon","zdeg","d");
	for(i=0; i<nlay; i++){
		for(w=0; w<ncols; w++){
			a = w*dx+dx/2;
			b = i*dz+dz/2;
			stla = asin(sin(inilat * DEG2RAD) * cos(a * KM2DEG * DEG2RAD) + cos(inilat * DEG2RAD) * sin(a * KM2DEG * DEG2RAD) * cos(azim * DEG2RAD)) * RAD2DEG;
			num = sin(a * KM2DEG * DEG2RAD) * sin(azim * DEG2RAD);
			den = cos(inilat * DEG2RAD) * cos(a * KM2DEG * DEG2RAD) - sin(inilat * DEG2RAD) * sin(a * KM2DEG * DEG2RAD) * cos(azim * DEG2RAD);
			stlo = inilon + atan2(num, den) * RAD2DEG;
			c = b*KM2DEG;
			if(i*dz>=depmin){
				fprintf(out_file, "%f,%f,%f,%f,%f,%f,%d\n", a,b,perfil[i*ncols+w],stla,stlo,c,(int)pnorm[i*ncols+w]);
			}
		}
	}
	fclose(out_file);

	/* Alliberar memòria */
	free(perfil);
	free(pnorm);
	free(absphase);
	free(pphase);

	return 0;
}
