/*
* $Id: rnd.c,v 1.7 2006/07/31 17:23:09 jms Exp $
*
* Revision History
* ===================
* $Log: rnd.c,v $
* Revision 1.7  2006/07/31 17:23:09  jms
* fix to parallelism problem
*
* Revision 1.6  2005/10/25 17:26:38  jms
* check in integration between microsoft changes and baseline code
*
* Revision 1.5  2005/10/14 23:16:54  jms
* fix for answer set compliance
*
* Revision 1.4  2005/09/23 22:29:35  jms
* fix to assume 64b support in the 32b RNG calls. Should speed generation, and corrects a problem with FK between Customer and Orders
*
* Revision 1.3  2005/03/04 21:43:23  jms
* correct segfult in random()
*
* Revision 1.2  2005/01/03 20:08:59  jms
* change line terminations
*
* Revision 1.1.1.1  2004/11/24 23:31:47  jms
* re-establish external server
*
* Revision 1.7  2004/04/08 17:34:15  jms
* cleanup SOLARIS/SUN ifdefs; now all use SUN
*
* Revision 1.6  2004/03/26 20:22:56  jms
* correct Solaris header
*
* Revision 1.5  2004/03/02 20:50:50  jms
* MP/RAS porting changes
*
* Revision 1.4  2004/02/18 16:37:33  jms
* add int32_t for solaris
*
* Revision 1.3  2004/02/18 16:26:49  jms
* 32/64 bit changes for overflow handling needed additional changes when ported back to windows
*
* Revision 1.2  2004/02/18 16:17:32  jms
* add 32bit specific changes to UnifInt
*
* Revision 1.1.1.1  2003/08/08 21:50:34  jms
* recreation after CVS crash
*
* Revision 1.3  2003/08/08 21:35:26  jms
* first integration of rng64 for o_custkey and l_partkey
*
* Revision 1.2  2003/08/07 17:58:34  jms
* Convery RNG to 64bit space as preparation for new large scale RNG
*
* Revision 1.1.1.1  2003/04/03 18:54:21  jms
* initial checkin
*
*
*/
/* 
 * RANDOM.C -- Implements Park & Miller's "Minimum Standard" RNG
 * 
 * (Reference:  CACM, Oct 1988, pp 1192-1201)
 * 
 * NextRand:  Computes next random integer
 * UnifInt:   Yields an long uniformly distributed between given bounds 
 * UnifReal: ields a real uniformly distributed between given bounds   
 * Exponential: Yields a real exponentially distributed with given mean
 * 
 */

#include "config.h"
#include <stdio.h>
#include <math.h>
#ifdef LINUX
#include <stdint.h>
#endif
#ifdef IBM
#include <inttypes.h>
#endif
#ifdef SUN
#include <inttypes.h>
#endif
#ifdef ATT
#include <sys/bitypes.h>
#endif
#ifdef WIN32
#define int32_t	__int32
#endif
#include "dss.h"
#include "rnd.h" 

/*
 *  Extensions to dbgen for generation of skewed data.
 *  default skew is 0 -- i.e. uniform distribution
 */
#define	PI (3.141597)
double	skew = 0;


char *env_config PROTO((char *tag, char *dflt));
void NthElement(DSS_HUGE, DSS_HUGE *);

void
dss_random(DSS_HUGE *tgt, DSS_HUGE lower, DSS_HUGE upper, long stream, double skew, long n)
{
	if( skew == 0 )
	{
		*tgt = UnifInt(lower, upper, stream);
		Seed[stream].usage += 1;
	}
	else
	{
		*tgt = SkewInt(lower, upper, stream, skew, n);
		Seed[stream].usage += 1;
	}

	return;
}

void
row_start(int t)	\
{
	int i;
	for (i=0; i <= MAX_STREAM; i++) 
		Seed[i].usage = 0 ; 
	
	return;
}

void
row_stop(int t)	\
	{ 
	int i;
	
	/* need to allow for handling the master and detail together */
	if (t == ORDER_LINE)
		t = ORDER;
	if (t == PART_PSUPP)
		t = PART;
	
	for (i=0; i <= MAX_STREAM; i++)
		if ((Seed[i].table == t) || (Seed[i].table == tdefs[t].child))
			{ 
			if (set_seeds && (Seed[i].usage > Seed[i].boundary))
				{
				fprintf(stderr, "\nSEED CHANGE: seed[%d].usage = %d\n", 
					i, Seed[i].usage); 
				Seed[i].boundary = Seed[i].usage;
				} 
			else 
				{
				NthElement((Seed[i].boundary - Seed[i].usage), &Seed[i].value);
#ifdef RNG_TEST
				Seed[i].nCalls += Seed[i].boundary - Seed[i].usage;
#endif
				}
			} 
		return;
	}

void
dump_seeds(int tbl)
{
	int i;

	for (i=0; i <= MAX_STREAM; i++)
		if (Seed[i].table == tbl)
#ifdef RNG_TEST
			printf("%d(%ld):\t%ld\n", i, Seed[i].nCalls, Seed[i].value);
#else
			printf("%d:\t%ld\n", i, Seed[i].value);
#endif
	return;
}

/******************************************************************

   NextRand:  Computes next random integer

*******************************************************************/

/*
 * long NextRand( long nSeed )
 */
DSS_HUGE
NextRand(DSS_HUGE nSeed)

/*
 * nSeed is the previous random number; the returned value is the 
 * next random number. The routine generates all numbers in the 
 * range 1 .. nM-1.
 */

{
    nSeed = (nSeed * 16807) % 2147483647;
    return (nSeed);
}

/******************************************************************

   UnifInt:  Yields an long uniformly distributed between given bounds

*******************************************************************/

/*
 * long UnifInt( long nLow, long nHigh, long nStream )
 */
DSS_HUGE
UnifInt(DSS_HUGE nLow, DSS_HUGE nHigh, long nStream)

/*
 * Returns an integer uniformly distributed between nLow and nHigh, 
 * including * the endpoints.  nStream is the random number stream.   
 * Stream 0 is used if nStream is not in the range 0..MAX_STREAM.
 */

{
    double          dRange;
    DSS_HUGE            nTemp,
		nRange;
    int32_t	nLow32 = (int32_t)nLow,
		nHigh32 = (int32_t)nHigh;
	
    if (nStream < 0 || nStream > MAX_STREAM)
        nStream = 0;
	
	if ((nHigh == MAX_LONG) && (nLow == 0))
	{
		dRange = DOUBLE_CAST (nHigh32 - nLow32 + 1);
		nRange = nHigh32 - nLow32 + 1;
	}
	else
	{
		dRange = DOUBLE_CAST (nHigh - nLow + 1);
		nRange = nHigh - nLow + 1;
	}

    Seed[nStream].value = NextRand(Seed[nStream].value);
#ifdef RNG_TEST
	Seed[nStream].nCalls += 1;
#endif
	nTemp = (long) (((double) Seed[nStream].value / dM) * (dRange));
    return (nLow + nTemp);
}



/*
 *  Extensions to dbgen for generation of skewed data.
 */

#define EPSILON  (0.0001)


/******************************************************************

   SkewInt: Yields a random number from a Zipfian distribution with 
   a specified skewVal. Skew is an integer in the range 0..4

*******************************************************************/

/*
 * long SkewInt( long nLow, long nHigh, long nStream )
 */
DSS_HUGE
SkewInt(DSS_HUGE nLow, DSS_HUGE nHigh, long nStream, double skewVal, long n)
{
	double	 zipf;
	double   dRange;
	DSS_HUGE nTemp;
	double	 multiplier;
	double	 Czn;
	DSS_HUGE numDistinctValuesGenerated;

	if (nStream < 0 || nStream > MAX_STREAM)
        	nStream = 0;

	/* check for validity of skewVal */
	if(skewVal < 0 || skewVal > 5)
	{
		zipf = 0; /* assume uniform */
	}
	else if(skewVal==5.0)
	{
		/* special case */
		/* check if any values have been generated for this column */
		if(NumDistinctValuesGenerated[nStream]==0)
		{
			/* pick a skew value to be used for this column*/
			zipf = (int) UnifInt(0,  4, 0);
			ColumnSkewValue[nStream] = zipf;
		}
		else
		{
			/* column skew value has been initialized before */
			zipf = ColumnSkewValue[nStream];
		}
	}
	else
	{
		/* skewVal is between 0 and 4 */
		zipf = skewVal;
	}
	
	/* If no values have been generated for this stream as yet, get multiplier */
	if(NumDistinctValuesGenerated[nStream]==0)
	{
		Multiplier[nStream] = GetMultiplier(n, zipf);
	}
	multiplier = Multiplier[nStream];

	/* 
	 * Check how many copies of the current value
	 * have already been generated for this stream.
	 * If we have generated enough, proceed to
	 * next value, and decide how many copies of
	 * this value should be generated.
	 */

	if(CurrentValueCounter[nStream] == CurrentValueTarget[nStream])
	{
		/* proceed to next value */
		if (nStream < 0 || nStream > MAX_STREAM)
			nStream = 0;
		if (nLow == nHigh)
			nTemp = nLow;
		else
		{
			if (nLow > nHigh)
			{
				nTemp = nLow;
				nLow = nHigh;
				nHigh = nTemp;
			}
			dRange = (double) (nHigh - nLow + 1);
			Seed[nStream].value = NextRand(Seed[nStream].value);
#ifdef RNG_TEST
			Seed[nStream].nCalls += 1;
#endif

			nTemp = (long) (((double) Seed[nStream].value / dM) * (dRange));
			nTemp += nLow;
			CurrentValue[nStream] = nTemp;
		}
	}
	else
	{
		/* return another copy of current value */
		nTemp = CurrentValue[nStream];
		CurrentValueCounter[nStream]++;

		return nTemp;
	}

	/*  
	 * check how many distinct values for this column
	 * have already been generated.
	 */
	numDistinctValuesGenerated = NumDistinctValuesGenerated[nStream] + 1;

	if(n<1)
	{
		n = (nHigh - nLow + 1);
	}

	Czn = n/multiplier;
	CurrentValueTarget[nStream]= 
		(long) (Czn/pow((double)numDistinctValuesGenerated, zipf));
	
	/* ensure that there is at least one value*/
	CurrentValueTarget[nStream] = (CurrentValueTarget[nStream] < 1) ? 1: CurrentValueTarget[nStream];
	CurrentValueCounter[nStream] = 1;
	NumDistinctValuesGenerated[nStream]++;

	return nTemp;
}


/******************************************************************%

   Returns the number of copies of the next value
   in the given distribution. Skew is an integer in the range 0..4

*******************************************************************/
/*
 *
 */
DSS_HUGE 
SkewIntCount(long nStream, double skewVal, long n)
/*
 * 
 */
{
	double	zipf;
	double	multiplier;
	double	Czn;
	long	numDistinctValuesGenerated;
	long	numCopies;

	if (nStream < 0 || nStream > MAX_STREAM)
        	nStream = 0;

	/* check for validity of skewVal */
	if(skewVal < 0 || skewVal > 4)
		zipf = 0; /* assume uniform */
	else if(skewVal==5.0)
	{
		/* column skew value has been initialized before */
		zipf = ColumnSkewValue[nStream];
	}
	else
	{
		/* skewVal is between 0 and 4 */
		zipf = skewVal;
	}

	/*  
	 * check how many distinct values for this column
	 * have already been generated.
	 */
	numDistinctValuesGenerated = NumDistinctValuesGenerated[nStream] + 1;

	multiplier = Multiplier[nStream];

	Czn = n/multiplier;
	numCopies= 
		(long) (Czn/pow((double)numDistinctValuesGenerated, zipf));
	/* ensure that there is at least one value*/
	CurrentValueTarget[nStream] = (CurrentValueTarget[nStream] < 1) ? 1: CurrentValueTarget[nStream];
	NumDistinctValuesGenerated[nStream]++;

	return numCopies;
}


/*
 * Round the number x to the nearest integer.
 */
double
round(double x)
{
	double fx = floor(x);
	double cx = ceil(x);

	if((x-fx) <= (cx-x))
		return fx;
	else
		return cx;
}

/*
 * Perform binary search to find D (number of distinct values for given n, zipf)
 *
 */
double 
SolveForMultipler(long n, double zipf)
{
	long lowD = 1;
	long highD = n;
	long mid;
	double sumMult;
	double nRowsPrime;
	long numRows;
	long i;

	while ( (highD - lowD) > 1)
	{
		mid = (highD + lowD) / 2;
		sumMult = 0.0;

		/* compute multiplier */
		for(i=1;i<=mid;i++)
		{
			sumMult += 1.0 / (pow((double) i, zipf));
		}

		/* compute number of rows generated for this mulitpler */
		numRows = 0;
		for(i=1;i<=mid;i++)
		{
			nRowsPrime = ((n / sumMult) / pow((double) i, zipf));
			numRows += (long) round(nRowsPrime);
		}
		if(((double)(n-numRows))/n  < EPSILON)
		{
			break;
		}

		if(numRows > n)
		{
			/* we overestimated numRows -- we need fewer distinct values */
			highD = mid;
		}
		else
		{
			/* underestimation of numRows -- need lower multiplier */
			lowD = mid;
		}
	}
	return sumMult;
}


double
GetMultiplier(long n, double zipf)
{
	double multiplier;
	double term;
	long i;

	if(zipf == 0.0)
		multiplier = n;
	else if(zipf==1.0)
		multiplier = log(n) + 0.577 ;
	else if(zipf==2.0)
		multiplier = (PI*PI)/6.0;
	else if(zipf==3.0)
		multiplier = 1.202;
	else if(zipf==4.0)
		multiplier = (PI*PI*PI*PI)/90.0;
	else
	{

		/* compute multiplier (must be within given bounds) */
		multiplier = 0;
		for(i=1;i<=n;i++)
		{
			term = 1.0 / pow((double)i, zipf);
			multiplier += term;

			/* if later terms add very little we can stop */
			if(term < EPSILON)
				break;
		}
	}
	return multiplier;
}

