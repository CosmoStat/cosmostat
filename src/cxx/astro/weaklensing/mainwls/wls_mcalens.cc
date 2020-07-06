/******************************************************************************************
**                   Copyright (C) 2009 by CEA  
*******************************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Olivier Fourt
**
**    Date:  26/06/09
**    
**    File:  mrsp_qu_to_eb.cc
**
********************************************************************************************
**
**    DESCRIPTION  convert Healpix polarized map from Q and U fields to E and B fields
**    ----------- 
**                 
*******************************************************************************************/

#include "HealpixClass.h"
#include "WLS_MassMapping.h"
#include "MRS_Sparse.h"
#include "Healpix_PowSpec.h"

// #include "map_alm_qu_eb.h"

char Name_file_Map1[512];
char Name_file_Map2[512];
char Name_file_Map3[512];
char Name_file_Map4[512];
char SignalPowSpec_FileName [512];
char Name_file_Noise_spec[512];
Bool OptS_PS = False;
// Bool OptN_PS = False;

bool input_spectrum = false;
bool NormALM = false;

extern int OptInd;
extern char *OptArg;
extern int GetOpt(int argc, char **argv, char *opts);

bool PolaFastALM = true;
Bool Verbose=False;

float  ZeroPadding=0.0;

float SigmaNoise=0.;
int NbrScale=0;
int Lmax=0;
int ALM_iter=0;

float NSigma=3.;

/***********************************************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options Map_G1_File_input Map_G2_File_input Map_File_output  \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "         [-n Number_of_Scales]\n");
    fprintf(OUTMAN, "             Number of scales in the wavelet transform.  \n");
    fprintf(OUTMAN, "         [-a N_ITER_ALM]\n");
    fprintf(OUTMAN, "             Number of iterations N_ITER_ALM for iterative
    fprintf(OUTMAN, "         [-l lmax]\n");
    fprintf(OUTMAN, "             Default is MIN(%d, 2*nside) \n",ALM_MAX_L);
    fprintf(OUTMAN, "         [-s Nsigma]\n");
    fprintf(OUTMAN, "             Nsigma Detection level.  \n");
    fprintf(OUTMAN, "         [-S Signal_PowSpec_FileName]\n");
    fprintf(OUTMAN, "             Input Signal Power Spectrum file name (for Wiener filtering). Default is automatically estimated. \n");
            
            
  //  fprintf(OUTMAN, "         [-P]\n");
  //  fprintf(OUTMAN, "             Input signal spectrum. Default is no, spectrum estimated from image in spectrum and noise spectrum.\n");
	fprintf(OUTMAN, "         [-v Verbose]\n");
    exit(-1);
}

/************************************************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void sinit( int argc, char *argv[] )
{
    int c;
    // get options
    while( (c = GetOpt( argc,argv,(char *) "S:n:s:i:l:fv" )) != -1 )
    {
		switch (c) 
        {
            case  'S':
                     if (sscanf(OptArg,"%s", SignalPowSpec_FileName) != 1)
                     {
                       fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
                      exit(-1);
                      }
                   OptS_PS = True;
                   break;
             case 'n':
                  if (sscanf(OptArg,"%d",&NbrScale) != 1)
                      {
                        fprintf(OUTMAN, "Error: bad number of scales value: %s\n", OptArg);
                        exit(-1);
                        }
                      break;
           case 's':
               if (sscanf(OptArg,"%f",&NSigma) != 1)
               {
                   fprintf(OUTMAN, "Error: bad number of scales value: %s\n", OptArg);
                   exit(-1);
               }
               break;
            case 'P':
                    input_spectrum = true;
                    if( sscanf(OptArg,"%s", Name_file_Signal_spec) != 1 )
                    {
                        fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
                        exit(-1);
                    }
                    break;
 		    case 'f':
				OptInd--;
           		PolaFastALM = false;
 		        break;
	        case 'l':
                if (sscanf(OptArg,"%d",&Lmax) != 1) 
                {
		            fprintf(OUTMAN, "Error: bad lmax  value: %s\n", OptArg);
		            exit(-1);
                }
                break;
     	    case 'v': Verbose = True; break;
		         	     	
 	     	default:
 	     		usage(argv);
 	     		break;
 		}
	} 

	//get optional input file names from trailing parameters and open files
    
    if( OptInd < argc )
    {
    	strcpy( Name_file_Map1, argv[OptInd] );
    	OptInd++;
    }
    else
    {
    	usage(argv);
    }
    
    if( OptInd < argc )
    {
    	strcpy( Name_file_Map2, argv[OptInd] );
    	OptInd++;
    }
    else
    {
    	usage(argv);
    }
    
    if( OptInd < argc )
    {
    	strcpy( Name_file_Map3, argv[OptInd] );
    	OptInd++;
    }
    else
    {
    	usage(argv);
    }
    
    if( OptInd < argc )
    {
        strcpy( Name_file_Map4, argv[OptInd] );
        OptInd++;
    }
    else
    {
        usage(argv);
    }

	// make sure there are not too many parameters
	if( OptInd < argc )
    {
		fprintf( OUTMAN, "Too many parameters: %s ...\n", argv[OptInd] );
		exit(-1);
	}
}
  
/*********************************************************************/



/*********************************************************************/

int main( int argc, char *argv[] )
{
    int k;
    fitsstruct Header, HD1;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);

	sinit(argc, argv);

    if (Verbose == True)
    { 
        cout << "# PARAMETERS: " << endl;
        if (NbrScale > 0)  cout << "# NbrScale = " <<  NbrScale << endl;
        cout << "# File Name in 1 = " << Name_file_Map1 << endl;
        cout << "# File Name in 2 = " << Name_file_Map2 << endl;
        cout << "# File Name Out 1 = " << Name_file_Map3 << endl;
        cout << "# File Name Out 2 = " << Name_file_Map4 << endl;
        if (PolaFastALM == false )  cout << "# Do not use fast alm calculation. " <<  Lmax << endl;
        if (Lmax > 0)  cout << "# Lmax = " <<  Lmax << endl;
        if (OptS_PS > 0)  cout << "# Signal power spectrum file name = " <<  SignalPowSpec_FileName << endl;

    }
            
	Hdmap Map_G1, Map_G2, Map_E, Map_B;
	
	Map_G1.read( Name_file_Map1 );
	Map_G2.read( Name_file_Map2 );
	
	if( Map_G1.Scheme() != Map_G2.Scheme() )
	{
		fprintf( OUTMAN, "Maps are not in the same scheme in: %s ...\n", argv[0] );
		exit(-1);
	}
    if (Verbose == True) Map_G1.info((char*) "Input G1");

	if( Map_G1.nside() != Map_G2.nside() )
	{
		fprintf( OUTMAN, "Maps have not the same Nside in: %s ...\n", argv[0] );
		exit(-1);
	}
    if (Verbose == True) Map_G2.info((char*) "Input G2");

	int Nside = Map_G1.nside();
    int Lmax_Tot = mrs_get_lmax (Lmax,  Nside,   ZeroPadding);
 	int Mmax = Lmax_Tot;
    // cout << " Lmax_Tot = " << Lmax_Tot << endl;
    
    Map_E.alloc(Nside);
    Map_B.alloc(Nside);
 	
    PowSpec SigPowSpec(1, Lmax_Tot);
    if (OptS_PS == True)
    {
          mrs_read_powspec(SignalPowSpec_FileName, SigPowSpec);
    }
    else
    {
        cout << " Error: a theoretical signal power spectrum is required. " << endl;
        exit(-1);
    }

    
    /*
    if( input_spectrum == true )
    {
        signal_spec.read( Name_file_Signal_spec );
    }
    noise_spec.read( Name_file_Noise_spec );
    */
    WLS_MassMapping MM;
    WLS_Field Kappa, RecG;
    MM.Verbose = Verbose;
    MM.alloc(Map_G1, Map_G2, 0.1, NbrScale);
    if (ALM_iter > 0) MM.set_niter(ALM_iter);
    
    // cout << "G2K" << endl;
    // MM.gamma2eb(MM.GammaData, Kappa);
    MM.sp_kaiser_squires(MM.GammaData, Kappa);

    Kappa.map_G1.info("Kappa");
    Kappa.map_G1.write( "xx_ks.fits" );// Write map E
    int NiterSparse=10;
    MM.sparse_reconstruction(MM.GammaData, Kappa, NSigma, NiterSparse);
    Kappa.map_G1.write( "xx_wt5sig.fits" );// Write map E

 /*
    Hdmap Map;
    //  Map.read(Name_Imag_In);
    // Nside = Map.Nside();
    // int Lmax = mrs_get_lmax (Lmax,  Nside, 0.0);
    // WT.wt_alloc(Nside, NbrScale, Lmax, DEF_MRS_ORDERING);
    // WT.hard_thresholding( Map, NSigma, SigmaNoise, UseMad, KillLastScale);

	Alm<xcomplex<double> > ALM_E(Lmax_Tot,Mmax), ALM_B(Lmax_Tot,Mmax);
	
	// map2alm_pol_iter_QU( Map_G1, Map_G2, ALM_E, ALM_B, 0, weight_T );
    Map_G1.info("Map_G1:");
    Map_G2.info("Map_G2:");
    map2alm_spin(Map_G1,Map_G2,ALM_E,ALM_B,2,weight_T,false);
    alm2map( ALM_E, Map_E );
    alm2map( ALM_B, Map_B );

    // WLS_MassMapping G;
    // ShearAlm AlmG;
    // G.alloc(Map_G1,Map_G2,SigmaNoise);
    
	Map_E.write( Name_file_Map3 );// Write map E
	Map_B.write( Name_file_Map4 );// Write map B
	
    Map_E.info("Map_True");

    Map_E -= Kappa.map_G1;
    Map_E.info("Resi mode E:");
*/
	exit(0);
}

