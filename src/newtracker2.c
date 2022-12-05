/*

  Frequency and phase tracker developed for use with Transcranial Alternating Current Stimulation
  primarily to phase-cancel Parkinsonian tremor.
  
  This program was developed using the (PCMCIA) NI-DAQCard 1200, which offers fast single channel
  read/write capabilities. A waveform is sampled and placed into a circular buffer. The most
  recent part of the circular buffer is then unwrapeed into a short "data" buffer and multipled
  by a sine and cosine function corresponding to the current tracking frequency. Magnitude and
  phase are then extracted and a phase-locked (and offset) sinusoid output. Frequency tracking
  is accomplished via a similar procedure, where a set of pre-specified sine and cosine sequences
  are multiplied by another separate (hence potentially larger) data buffer. The frequency with
  the highest magnitude response takes over as tracking frequency. To avoid congesting the
  processor, buffering and subsequent testing for this is performed over sequential sampling
  iterations. For instance, with a minimum frequency (fmin) of 4, a maximum (fmax) of 6, and a
  testing resolution (fres) of 0.1, there are 21 pre-computed sine and cosine waves. On the
  first iteration the frequency tracking buffer is unwrapped from the circular buffer, each of
  the 21 frequencies are then tested over the following 21 sampling loops. On the next sample,
  the largest magnitude is assessed and selected as tracking frequency. Frequency tracking can
  be delayed by specifying the trecentre parameter, which corresponds to the number of dead
  iterations before refilling the tracking buffer.
  
  Phase-offset read from Channel 1 ( ChanIn + 1 ), where 0-5 V converts to 0-2*pi phase difference
  
  NeuroConn calibration, 2V input = 4mA output
  
  This program reads parameters from "newtracker.ini" as parameter-value pairs. Default exist
  for all parameters if the file / any parameters are missing. No equals sign should be present,
  merely a space. For instance, "CBUFFER_LEN 65536".
  
  Parameters (with default values and type)
		CBUFFER_LEN	65536	Length of circular buffer. Larger is better to avoid unwrapping.
		DBUFFER_LEN	500		Length of data buffer. Defines the number of data points to assess for phase estimation.
		FBUFFER_LEN	5000	Length of frequency tracking buffer. Longer is more conistent estimate, but takes longer to switch.
		rate		1000	Sampling rate (ideal)
		trecentre	0		Cycles to wait between frequency tracking updates.
		threshold	0.0		Amplitude threshold (sine amplitude)
		flow		4		Lowest frequency to test.
		fhigh		6		Highest frequency to test.
		fres		0.1		Frequency resolution.
		benchmark	0		Benchmark (switch, 0/1). Command line option -b more convenient, only standard test.
								0	Off
								1	Standard benchmark (operating as normal)
								2	Force frequency buffer refill on every loop
								3	Force frequency test (no buffer refil) on every loop (slowest part)
		verbose		1		Verbose mode (switch 0/1).
		ChanIn		0		Input channel specifier
		ChanOut		0		Output channel specifier
		ramptime 	1.0		Ramp time (secs) for onset. Off-time is twice as quick.
		ampout		2.0		Desired output current (NeuroConn)
		outputmode	2		Output mode
								0	Input (straight-through)
								1	Phase (in radians)
								2	Phase-offset sinusoid
		ftracking	0.0		Frequency tracking
								0	Turned ON
								n	Turned OFF; fix frequency at n Hz
		ptracking	1.0		Phase tracking
								0	None - freerun
								1	Tracking (phase from secondary input channel)
							   -n   Forced (tracked) drift with n sec drift rate
		driftdir	1		Forced drift direction (1 or -1)
							   
	Command line options
		-c <filename>		Configuration file. Default "newtracker.ini".
							Uses default values if none specified otherwise.
		-b					Benchmark testing (overrides parameter specification)
  
  Note: It is important, before running, to benchmark the computer and acquisition device since
  the circular buffer is populated sequentially and assumes that a MAXIMUM of one sample of the
  idealised sampling rate has passed. If the setup is too slow, the tracking frequency will not
  be estimated properly. There is a benchmark option in the parameters file, or, for convenience,
  the option -b can be specified at the command line. If your sampling rate is set too high, you
  will be informed at this point.
  
  Developed by John-Stuart Brittain
	   Created Mar-2013
  Last updated Mar-2013

*/

#include "nidaqex.h"
#include <stdio.h>
#include <math.h>
#include <windows.h>		// Make use of high precision timers
#include <conio.h>

#define DEBUG_ITERATIONS	1000
#define PI 					3.142

double fmin( double a, double b ) { return (a<b)?a:b; }
double fmax( double a, double b ) { return (a>b)?a:b; }
double getElapsedTimeMsecs( LARGE_INTEGER t1, LARGE_INTEGER t2, LARGE_INTEGER frequency ) {
	return (t2.QuadPart - t1.QuadPart)*1000.0/frequency.QuadPart;
}

int main( int argc, char* argv[] ) {

	const double twoPI	= 2*PI;

	/* Debug parameters */
	int debug_benchmark	= 0;
	int verbose 		= 1;
	
    /* NI-DAQ parameters */
    i16 iStatus, iRetVal;
    i16 iDevice 		= 1;
    i16 iChanIn 		= 0;
	i16 iChanOut 		= 0;
    i16 iGain 			= 1;
    f64 dVoltageIn		= 0.0;
	f64 dVoltageOut		= 0.0;
	f64 dVoltagePhase   = 0.0;				// Current reading
	f64 dVoltagePhase0  = 0.0;				// Previous (smoothed) reading
    i16 iIgnoreWarning	= 1;
	
	/* Parameters */
	int    rate			= 1000;
	int    trecentre 	= 0;				// Timer for recentring analysis frequency (in samples)
	double threshold 	= 0.0;
	double flow 		= 4.9;
	double fhigh 		= 5.1;
	double fres 		= 0.1;
	int    CBUFFER_LEN	= 65536;
	int    DBUFFER_LEN	= 250;
	int    FBUFFER_LEN	= 5000;
	int    outputmode   = 2;				// 0=input, 1=phase, 2=sine
	double ramptime     = 1.0;				// Ramp-up time
	double ampout 		= 2.0;				// 2 mA output
	double ftracking    = 0.0;				// 0=track, +ve=fixed_at_that_frequency
	double ptracking	= 1.0;				// 0=freerun, 1=tracking, -ve=forced_drift
	int    driftdir     = 1;				// 1=fwd, -1=bwd
	char  *fileini;
	
	/* Internal variables */
    int    j, c, k, startpos;				// Loop variable
	int    cpos;							// Position in circular buffer
	double dt, dmean, dmax;
	double magin, anglein, cossum, sinsum, fcc;
	double phaseoffset, output_correction, ramp, rampinc, phaseinc;
	double *fmagin;
	int fc, fc0, fcount, ftotal;			// Centre analysis frequency (index)
	
	/* Setup buffers (f64 defined as double) */
	f64 *cbuffer, *dbuffer, *fbuffer;		// Buffers
	double **sinwave, **coswave;			// Trigonometric primitives
	int *wavesamples;
	LARGE_INTEGER frequency, td, ti, t0;	// High precision timers (Windows specific)
	FILE *fid;
	char param[80], value[80], paramcomp[80];
	fileini = NULL;
	
	/* Parse command line for initialisation file specifier */
	for (k=1; k<argc; k++)
		if (argv[k][0]=='-') {
			if (argv[k][1]=='c') {
				k+=1;
				if (argc<(k+1)) {
					printf("No configuration file specified.");
					return -1;
				}
				fileini = argv[k];
			}
		}
	
	/* Parse INI file */
	if (fileini == NULL)
		// Cannot free fileini if pointing to argv later, so avoid allocating for default filename
		fid = fopen("newtracker.ini","r");
	else
		fid = fopen(fileini,"r");
	if ( fid != NULL ) {
		printf("Load settings\n");
		while(!feof(fid)) {
			fscanf(fid,"%s %s",param,value);
			printf(" %s = %s\n",param,value);
			strcpy(paramcomp,"CBUFFER_LEN");     if (!strcmp(param,paramcomp)) { CBUFFER_LEN = atoi(value);     continue; };
			strcpy(paramcomp,"DBUFFER_LEN");     if (!strcmp(param,paramcomp)) { DBUFFER_LEN = atoi(value);     continue; };
			strcpy(paramcomp,"FBUFFER_LEN");     if (!strcmp(param,paramcomp)) { FBUFFER_LEN = atoi(value);     continue; };
			strcpy(paramcomp,"rate");    	     if (!strcmp(param,paramcomp)) { rate = atoi(value); 		    continue; };
			strcpy(paramcomp,"trecentre");       if (!strcmp(param,paramcomp)) { trecentre = atoi(value); 	    continue; };
			strcpy(paramcomp,"threshold");       if (!strcmp(param,paramcomp)) { threshold = atof(value); 	    continue; };
			strcpy(paramcomp,"flow");            if (!strcmp(param,paramcomp)) { flow = atof(value); 		    continue; };
			strcpy(paramcomp,"fhigh");           if (!strcmp(param,paramcomp)) { fhigh = atof(value); 		    continue; };
			strcpy(paramcomp,"fres");            if (!strcmp(param,paramcomp)) { fres = atof(value); 		    continue; };
			strcpy(paramcomp,"benchmark"); 		 if (!strcmp(param,paramcomp)) { debug_benchmark = atoi(value); continue; };
			strcpy(paramcomp,"verbose");         if (!strcmp(param,paramcomp)) { verbose = atoi(value); 	    continue; };
			strcpy(paramcomp,"ChanIn");          if (!strcmp(param,paramcomp)) { iChanIn = (i16)atoi(value);    continue; };
			strcpy(paramcomp,"ChanOut");         if (!strcmp(param,paramcomp)) { iChanOut = (i16)atoi(value);   continue; };
			strcpy(paramcomp,"outputmode");      if (!strcmp(param,paramcomp)) { outputmode = atoi(value);      continue; };
			strcpy(paramcomp,"ramptime");        if (!strcmp(param,paramcomp)) { ramptime = atof(value);        continue; };
			strcpy(paramcomp,"ampout");          if (!strcmp(param,paramcomp)) { ampout = atof(value);          continue; };
			strcpy(paramcomp,"ftracking");       if (!strcmp(param,paramcomp)) { ftracking = atof(value);       continue; };
			strcpy(paramcomp,"ptracking");       if (!strcmp(param,paramcomp)) { ptracking = atof(value);       continue; };
			strcpy(paramcomp,"driftdir");        if (!strcmp(param,paramcomp)) { driftdir = atoi(value);        continue; };
			printf("Unknown parameter ignored: %s %s\n   - Press ENTER to confirm -\n",param,value);
			getchar();
		}
		fclose(fid);
	} else {
		if (fileini != NULL) {
			printf("Cannot open configuration file %s\n",fileini);
			return -1;
		}
	}
	
	/* Parse command line options */
	for (k=1; k<argc; k++)
		if (argv[k][0]=='-') {
			// Quick option to benchmark
			if (argv[k][1]=='b') {
				if (strlen(argv[k])>2)
					debug_benchmark = atoi( &argv[k][2] );
				else
					debug_benchmark = 1;
				continue;
			}
			if (argv[k][1]=='c') {
				// Config file (also skip next arg as filename)
				k += 1;
				continue;
			}
			printf("Invalid command line option %s\n",argv[k]);
			return -1;
		} else {
			printf("Invalid command line option %s\n",argv[k]);
			return -1;
		}
	
	/* Set variables after INI read */
	dt = 1000.0/rate;				// Sampling interval (msecs)
	fcount = -trecentre;
	ftotal = (int)ceil((fhigh-flow)/fres+1);
	phaseoffset = 0; fc0 = 0; cpos = 0; ramp = 0;
	if ( ramptime > 0 ) { rampinc = 1/(rate*ramptime); ramp = 0; }
	else				{ rampinc = 0; ramp = 1; }
	output_correction = ampout/2;
	
	/* Dynamic array allocation */
	cbuffer     = (f64*)malloc(CBUFFER_LEN*sizeof(f64));
	dbuffer     = (f64*)malloc(DBUFFER_LEN*sizeof(f64));
	fbuffer     = (f64*)malloc(FBUFFER_LEN*sizeof(f64));
	wavesamples = (int*)malloc(ftotal*sizeof(int));
	fmagin      = (double*)malloc(ftotal*sizeof(double));
	sinwave     = (double**)malloc(ftotal*sizeof(double*));
	coswave     = (double**)malloc(ftotal*sizeof(double*));
	for (k=0; k<ftotal; k++) {
		sinwave[k] = (double*)malloc(DBUFFER_LEN*sizeof(double));
		coswave[k] = (double*)malloc(DBUFFER_LEN*sizeof(double));
	}
	
	/* Get timer conversion */
	QueryPerformanceFrequency(&frequency);
	
	/* Initialise data buffer */
	for ( c=0; c<CBUFFER_LEN; c++ ) {
		cbuffer[c] = 0;
	}
	
	/* Construct trigonometric primitives */
	for ( k=0; k<ftotal; k++ ) {
		for ( c=0; c<DBUFFER_LEN; c++ ) {
			fcc = flow + fres*k;
			coswave[k][c] = cos(twoPI*c*fcc/rate);
			sinwave[k][c] = sin(twoPI*c*fcc/rate);
		}
		// Adjust for wave-length (needed for continuity)
		wavesamples[k] = (int)floor(rate/fcc);
	}
	
	/* Tracking frequency */
	if (ftracking==0)
		// Default to middle frequency
		fc = (int)ceil(ftotal/2);
	else {
		// Find closest frequency
		dmax = 3*fhigh;
		for ( k=0; k<ftotal; k++ )
			// Difference between sine frequency and fixed tracking frequency
			if ( fabs((flow + fres*k) - ftracking) < dmax ) {
				dmax = fabs((flow + fres*k) - ftracking);
				fc = k;
			}
	}
	
	/* Determine phase tracking increment */
	if ( ptracking == 0 ) {
		// No tracking - increment at sampling rate for given frequency (initialiser as also done online)
		phaseinc = (flow + fres*fc)*twoPI/rate;
	} else if ( ptracking < 0 ) {
		if ( ( driftdir != 1 ) & ( driftdir != -1 ) ) {
			printf(" Illegal drift direction ( must be +/- 1 )\n");
			return -1;
		};
		// Forced tracking - determine drift rate
		phaseinc = driftdir*twoPI/(rate*ptracking);
	}
	
	/* Start loop */
	if ( verbose ) {
		printf("Running [press any key to exit] ...\n");
		printf("Tracking frequency: %1.1f Hz\b\b\b",flow+fres*fc);
	}
	QueryPerformanceCounter(&td);
	QueryPerformanceCounter(&ti);
	for ( j=0; j<DEBUG_ITERATIONS; (debug_benchmark?j++:j) ) {
	
		/* Select phase offset method */
		if ( ptracking == 1 ) {
			// Get phase-offset from second channel //
			iStatus = AI_VRead(iDevice, iChanIn+1, iGain, &dVoltagePhase);
			// Check for errors
			iRetVal = NIDAQErrorHandler(iStatus, "AI_VRead", iIgnoreWarning);
			if ( iRetVal != 0 ) {
				printf(" Read error %d\n", iRetVal);
				return iRetVal;
			}
			// Smooth phase input to remove spikes //
			if ( 0 ) {
				dVoltagePhase = 0.1*dVoltagePhase + 0.9*dVoltagePhase0;
				dVoltagePhase0 = dVoltagePhase;
			}
			// Update phase-offset (input 0-5 V)
			phaseoffset = -dVoltagePhase*1.2568; // *(2*PI/5)
		} else if ( ptracking == 0 ) {
			// Increment phase offset at current frequency but do not align to input phase
			phaseoffset = fmod( phaseoffset + (flow + fres*fc)*twoPI/rate, twoPI );
		} else if ( ptracking < 0 ) {
			// Forced drift //
			phaseoffset = fmod( phaseoffset + phaseinc, twoPI );
		} else {		// ( ptracking > 0 ) but not 1; illegal
			// Error checked because it is easy to specify an illegal option //
			printf(" Illegal phase tracking options\n");
			return -2;
		};
		
		/* Acquire input */
		iStatus = AI_VRead(iDevice, iChanIn, iGain, &dVoltageIn);
		// Check for errors
		iRetVal = NIDAQErrorHandler(iStatus, "AI_VRead", iIgnoreWarning);
		if ( iRetVal != 0 ) {
			printf(" Read error %d\n", iRetVal);
			return iRetVal;
		};
		// Store data in circular buffer and increment counter
		cbuffer[cpos++] = dVoltageIn;
		// Reset circular position counter
		if ( cpos == CBUFFER_LEN )
			cpos = 0;
		
		/* Unwrap data into local data buffer */
		c = 0;		// Data buffer position
		startpos = cpos-DBUFFER_LEN;
		if (startpos>=0)
			// Take whole chunk - easy!
			for ( k=startpos; k<(startpos+DBUFFER_LEN+1); k++ )
				dbuffer[c++] = cbuffer[k];
		else {
			// Combine start, end chunks
			for ( k=(CBUFFER_LEN-startpos); k<CBUFFER_LEN; k++ )	// First chunk (back end)
				dbuffer[c++] = cbuffer[k];
			for ( k=0; k<cpos; k++ )								// Second chunk (front end)
				dbuffer[c++] = cbuffer[k];
		}
		/* Calculate and subtract mean */
		dmean = 0;
		for ( c=0; c<DBUFFER_LEN; c++ )
			dmean += dbuffer[c];
		dmean = dmean / DBUFFER_LEN;
		for ( c=0; c<DBUFFER_LEN; c++ )
			dbuffer[c] -= dmean;
		
		//
		/* Single-frequency decomposition */
		//
		
		/* Multiply waveforms by trigonometric functions */
		cossum = 0; sinsum = 0;
		for ( c=0; c<DBUFFER_LEN; c++ ) {
			// Flip data buffer for phase at last point
			cossum += dbuffer[DBUFFER_LEN-c-1]*coswave[fc][c];
			sinsum += dbuffer[DBUFFER_LEN-c-1]*sinwave[fc][c];
		}
	
		// Determine magnitude and angle
		magin = 2*sqrt( pow(cossum/DBUFFER_LEN,2) + pow(sinsum/DBUFFER_LEN,2) );
		if ( magin > 0 ) anglein = atan2( sinsum, cossum );
		else 		     anglein = 0;
		
		//
		/* Recentre analysis frequency */
		//
		
		if ( ftracking == 0 ) {
			
			/* Re-centering (timed by iteration count) */
			if ( fcount == -1 ) {
			
				/* Extract extended epoch from circular buffer */
				c = 0;
				startpos = cpos-FBUFFER_LEN;
				if (startpos>=0)
					/* Take whole chunk - easy! */
					for ( k=startpos; k<cpos; k++ )
						fbuffer[c++] = cbuffer[k];
				else {
					/* Combine start, end chunks */
					for ( k=(CBUFFER_LEN+startpos); k<CBUFFER_LEN; k++ )	// First chunk (back end)
						fbuffer[c++] = cbuffer[k];
					for ( k=0; k<cpos; k++ )								// Second chunk (front end)
						fbuffer[c++] = cbuffer[k];
				}
				/* Calculate and subtract mean */
				dmean = 0;
				for ( c=0; c<FBUFFER_LEN; c++ )
					dmean += fbuffer[c];
				dmean = dmean / FBUFFER_LEN;
				for ( c=0; c<FBUFFER_LEN; c++ )
					fbuffer[c] -= dmean;
			};
			if ( fcount >= 0 ) {
			
				if ( fcount < ftotal ) {
				
					// Single sine/cosine interrogation
					cossum = 0; sinsum = 0;
					for ( c=0; c<FBUFFER_LEN; c++ ) {
						// Multiply and running sum
						cossum += fbuffer[c]*coswave[fcount][c % wavesamples[fcount]];
						sinsum += fbuffer[c]*sinwave[fcount][c % wavesamples[fcount]];
					}
					fmagin[fcount] = 2*sqrt( pow(cossum/FBUFFER_LEN,2) + pow(sinsum/FBUFFER_LEN,2) );
				
				} else {
				
					// Find maximum element
					dmax = 0;
					for ( k=0; k<ftotal; k++ ) {
						if (fmagin[k]>dmax) {
							dmax = fmagin[k];
							fc0 = k;
						}
					}
					// Change frequency if result larger than threshold
					if ( fmagin[fc0] > threshold ) fc = fc0;
					else 						   fc = (int)ceil( ftotal/2 );
					// In-case trecentre is 0; need -1 for buffering but fcount++ below
					fcount = min(-2,-trecentre);
					
					// Debug output
					if ( verbose )
						printf("\b\b\b\b%4.1f",flow+fres*fc);
				}
			}
			fcount++;
			
		}
		
		/* Output phase-offset sinusoid */
		switch ( outputmode ) {
			case 0:			// Input
				dVoltageOut = dVoltageIn;					// Voltage output
				//dVoltageOut = dbuffer[DBUFFER_LEN-1];		// Last element of data buffer
				//dVoltageOut = cbuffer[cpos-1];			// Current element of circular buffer
				break;
			case 1:			// Phase
				dVoltageOut = anglein;
				break;
			case 2:			// Sinusoid
				// Adjust ramp
				if ( magin > threshold ) { if (ramp<1) ramp = fmin( 1, ramp+rampinc );   }
				else					 { if (ramp>0) ramp = fmax( 0, ramp-2*rampinc ); }				
				if ( ptracking != 0 )
					dVoltageOut = ramp*output_correction*cos( anglein + phaseoffset );
				else
					dVoltageOut = ramp*output_correction*cos( phaseoffset );
				break;
			case 3:			// Debug - phase offset
				dVoltageOut = dVoltagePhase;
				break;
		}
		iStatus = AO_VWrite(iDevice, iChanOut, dVoltageOut);
		// Check for errors
		iRetVal = NIDAQErrorHandler(iStatus, "AO_VWrite", iIgnoreWarning);
		if ( iRetVal != 0 ) {
			printf(" Write error %d\n", iRetVal);
			return iRetVal;
		}
		
		/* Check debug status */
		if ( debug_benchmark == 2 )
			fcount = -1;
		else if ( debug_benchmark == 3 )
			fcount = 0;
		
		/* Look for escape */
		if ( kbhit() )
			break;
		
		/* Rate limit loop (but do not give up control of processor) */
		QueryPerformanceCounter(&t0);
		if ( !debug_benchmark )
			if ( t0.QuadPart >= ti.QuadPart )		// Account for wrap-around of the system clock
				while ( getElapsedTimeMsecs( ti, t0, frequency ) < dt ) QueryPerformanceCounter(&t0);
		QueryPerformanceCounter(&ti);
		
    };
	
	/* Reset output to zero */
	iStatus = AO_VWrite(iDevice, iChanOut, 0);
	// Check for errors
	iRetVal = NIDAQErrorHandler(iStatus, "AO_VWrite", iIgnoreWarning);
	if ( iRetVal != 0 ) {
		printf(" Write error %d\n", iRetVal);
		return iRetVal;
	}

	/* Report loop time */
	if (debug_benchmark) {
		QueryPerformanceCounter(&t0);
		printf("\n%f msecs per iteration\n",getElapsedTimeMsecs( td, t0, frequency )/DEBUG_ITERATIONS);
		if ( getElapsedTimeMsecs( td, t0, frequency )/DEBUG_ITERATIONS < dt )
			printf(" This is OK with the current %d Hz sampling rate\n",rate);
		else
			printf("\n\tSTOP! STOP! STOP! STOP! STOP!\n This is NOT OK with the current %d Hz sampling rate\n\tSTOP! STOP! STOP! STOP! STOP!",rate);
		getchar();
	}
	
	/* Deallocate arrays */
	for (k=0; k<ftotal; k++) {
		free((void*)sinwave[k]);
		free((void*)coswave[k]);
	}
	free((void*)sinwave);
	free((void*)coswave);
	free((void*)wavesamples);
	free((void*)fmagin);
    
	/* Pause for keyboard input before exist */
	if ( verbose )
		printf("\nFinished.");
	
	/* Return without errors */
	return 0;
}
