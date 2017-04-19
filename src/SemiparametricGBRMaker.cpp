// Minimalistic GRBMaker for semiparametric regression
// Based on work by Jean-Baptiste Sauvan and Josh Bendavid

#include "SemiparametricGBRMaker.h"
#include "Utilities.h"
// #include "GBRTrainer.h"
// #include "GBRForest.h"
// #include "GBRApply.h"

#include <TFile.h>
#include <TChain.h>
#include <TCut.h>
#include <TKey.h>
#include "TString.h"

#include "RooArgList.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooConstVar.h"
#include "RooWorkspace.h"
#include "HiggsAnalysis/GBRLikelihood/interface/RooHybridBDTAutoPdf.h"
#include "HiggsAnalysis/GBRLikelihood/interface/RooDoubleCBFast.h"
#include "HiggsAnalysis/GBRLikelihood/interface/HybridGBRForestFlex.h"
#include "CondFormats/EgammaObjects/interface/GBRForestD.h"

#include <iostream>
#include <sstream>
#include <time.h>

// For reading the paramater file
#include <TEnv.h>

// For interpreting the parameters
#include <TSystem.h>
#include <algorithm>

// For debugging
#include <typeinfo>


using namespace std;
using namespace RooFit;



/*****************************************************************/
SemiparametricGBRMaker::SemiparametricGBRMaker()
/*****************************************************************/
{
}

/*****************************************************************/
SemiparametricGBRMaker::~SemiparametricGBRMaker()
/*****************************************************************/
{
}

/*****************************************************************/
void SemiparametricGBRMaker::readParameterFile( const string& parFileName )
/*****************************************************************/
{
    // ########################################
    // Read the parameter file
    // ########################################

    TEnv params;
    int status = params.ReadFile(parFileName.c_str(),EEnvLevel(0));
    if(status!=0) 
    {
        cout << "FATAL: SemiparametricGBRMaker::readParameterFile(): Cannot read file " << parFileName << "\n"; 
    }
    else {
        cout << "INFO: Reading parameters from " << parFileName << " was successful" << endl;
    }

    par_Name            = params.GetValue( "Name",             "" );
    par_OutputDirectory = params.GetValue( "OutputDirectory",  "" );
    par_InputFiles      = params.GetValue( "InputFiles",       "" );
    par_Tree            = params.GetValue( "Tree",             "" );
    par_Options         = params.GetValue( "Options",          "" );
    par_Variables       = params.GetValue( "Variables",        "" );
    par_Target          = params.GetValue( "Target",           "" );
    par_Cut             = params.GetValue( "Cut",              "" );

    par_mu_DownLimit    = params.GetValue( "mu_DownLimit",     0. );
    par_mu_UpLimit      = params.GetValue( "mu_UpLimit",       0. );
    par_sigma_DownLimit = params.GetValue( "sigma_DownLimit",  0. );
    par_sigma_UpLimit   = params.GetValue( "sigma_UpLimit",    0. );
    par_n1_DownLimit    = params.GetValue( "n1_DownLimit",     0. );
    par_n1_UpLimit      = params.GetValue( "n1_UpLimit",       0. );
    par_n2_DownLimit    = params.GetValue( "n2_DownLimit",     0. );
    par_n2_UpLimit      = params.GetValue( "n2_UpLimit",       0. );
    par_alpha1          = params.GetValue( "alpha1",           0. );
    par_alpha2          = params.GetValue( "alpha2",           0. );

    cout << endl;
    cout << "par_Name             = " << par_Name << endl ;
    cout << "par_OutputDirectory  = " << par_OutputDirectory << endl ;
    cout << "par_InputFiles       = " << par_InputFiles << endl ;
    cout << "par_Tree             = " << par_Tree << endl ;
    cout << "par_Options          = " << par_Options << endl ;
    cout << "par_Variables        = " << par_Variables << endl ;
    cout << "par_Target           = " << par_Target << endl ;
    cout << "par_Cut              = " << par_Cut << endl ;
    cout << "par_mu_DownLimit     = " << par_mu_DownLimit    << endl;
    cout << "par_mu_UpLimit       = " << par_mu_UpLimit      << endl;
    cout << "par_sigma_DownLimit  = " << par_sigma_DownLimit << endl;
    cout << "par_sigma_UpLimit    = " << par_sigma_UpLimit   << endl;
    cout << "par_n1_DownLimit     = " << par_n1_DownLimit    << endl;
    cout << "par_n1_UpLimit       = " << par_n1_UpLimit      << endl;
    cout << "par_n2_DownLimit     = " << par_n2_DownLimit    << endl;
    cout << "par_n2_UpLimit       = " << par_n2_UpLimit      << endl;
    cout << "par_alpha1           = " << par_alpha1          << endl;
    cout << "par_alpha2           = " << par_alpha2          << endl;
    cout << endl;


    // ########################################
    // Derived variables
    // ########################################

    cout << "INFO: Determining variables derived from parameter file" << endl;

    // ========================================
    // Fill the chain and list of trees and rootfiles

    chain_ = new TChain( par_Tree.c_str() );

    vector<string> rootFilesVector;
    tokenize( par_InputFiles, rootFilesVector, ":" );
    vector<string>::iterator rootFile = rootFilesVector.begin();
    vector<string>::iterator lastRootFile = rootFilesVector.end();

    for( ; rootFile != lastRootFile ; ++rootFile ) {
        cout << "INFO: Adding " << *rootFile << " to chain" << endl;
        TFile* rootFp = TFile::Open( rootFile->c_str() );
        if( !rootFp || !rootFp->IsOpen() ) {
            cout << "FATAL: SemiparametricGBRMaker::init(): Cannot open input file " << *rootFile << "\n";
            }
        TTree* tree = (TTree*)rootFp->Get( par_Tree.c_str() );
        if( !tree ){
            cout << "FATAL: SemiparametricGBRMaker::init(): Cannot find regression tree " << par_Tree << " in " << *rootFile << "\n";
            }

        // Set class variables
        rootFps_.push_back( rootFp );
        chain_->Add( rootFile->c_str() );
        trees_.push_back( tree );
        }
    

    // ========================================
    // Set outputFile_ to proper path

    stringstream outputFileStream;
    outputFileStream <<  par_OutputDirectory  <<  "/"  <<  par_Name  <<  "_results.root";
    outputFile_ = outputFileStream.str();

    // Create outputFile_ (why here?)
    TFile* outputRootFp = TFile::Open( outputFile_.c_str(), "RECREATE" );
    if(!outputRootFp || !outputRootFp->IsOpen()) {
        cout << "FATAL: SemiparametricGBRMaker::init(): Cannot open output file " << outputFileStream.str() << "\n";
        }
    outputRootFp->Close();


    // ========================================
    // Set variables_ (there has got to be a better way to do this)

    vector<string> variablesNotPointer;
    tokenize( par_Variables , variablesNotPointer, ":" );
    vector<string>::iterator variable  = variablesNotPointer.begin();
    vector<string>::iterator lastVariable = variablesNotPointer.end();
    for( ; variable != lastVariable ; ++variable ){
        cout << "INFO: Adding variable " << *variable << endl;
        variables_.push_back( *variable );
        }
    variablesNotPointer.clear();

}


/*****************************************************************/
void SemiparametricGBRMaker::close()
/*****************************************************************/
{
    vector<TFile*>::iterator rootFp = rootFps_.begin();
    vector<TFile*>::iterator lastRootFp = rootFps_.end();
    for( ; rootFp != lastRootFp ; ++rootFp )
    {
        if(*rootFp)
        {
            (*rootFp)->Close();
            delete *rootFp;
         }
    }
    rootFps_.clear();

    // Maybe close chain / trees / rootfiles here

    }


/*****************************************************************/
void SemiparametricGBRMaker::run()
/*****************************************************************/
{

    // ========================================
    // create RooRealVars for each input variable and add to a RooArgList

    RooArgList vars;
    for (unsigned int ivar=0; ivar < variables_.size(); ++ivar) {
        RooRealVar* var = new RooRealVar(
            TString::Format("var_%i",ivar),
            variables_.at(ivar).c_str(),
            0.);
        vars.addOwned(*var);
        }

    // Clone list (this is input vars only)
    RooArgList condvars(vars);

    RooRealVar* targetvar  = new RooRealVar( "targetvar", par_Target.c_str(), 1. );
    vars.addOwned(*targetvar);


    // ========================================
    // Do weights and cuts and create dataset

    // Retrieve event weight
    string weightStr;
    vector<string> optionTokens;
    tokenize( par_Options, optionTokens, ":");
    for(const auto& token : optionTokens)
    {
        vector<string> tagAndValue;
        tokenize(token, tagAndValue, "=");
        string tag = tagAndValue[0];
        string value = tagAndValue[1];
        if(tag=="EventWeight")
        {
            weightStr = value;
            cout << "INFO: EventWeight = " << weightStr << "\n";
        }
    }

    RooRealVar weight_and_cuts( "weight_and_cuts", "", 1. );
    TCut weightTCut( weightStr.c_str() );
    TCut eventTCut( par_Cut.c_str() );
    weight_and_cuts.SetTitle( weightTCut * eventTCut );
    cout << "INFO: weight_and_cuts = " << weight_and_cuts.GetTitle() << endl;

    cout << "INFO: Creating the dataset" << endl;
    RooDataSet *hdata = RooTreeConvert::CreateDataSet( "hdata", chain_, vars, weight_and_cuts );
    cout << endl;


    // ========================================
    // Define regression parameters and create pdf

    cout << "INFO: Creating the variables and the pdf" << endl;

    // sigma of the CB
    RooRealVar sigmavar(
        "sigmavar",
        "sigmavar",
        0.01
        );
    sigmavar.setConstant(false);

    // mu of the CB
    RooRealVar muvar(
        "muvar",
        "muvar",
        1.00
        );
    sigmavar.setConstant(false);

    // n1 of the CB
    RooRealVar n1var(
        "n1var",
        "n1var",
        3.00
        );
    n1var.setConstant(false);

    // n2 of the CB
    RooRealVar n2var(
        "n2var",
        "n2var",
        3.00
        );
    n2var.setConstant(false);


    // JBS: define non-parametric functions for each regressed parameter
    RooGBRFunctionFlex *sigmafunc = new RooGBRFunctionFlex( "sigmafunc", "" );
    RooGBRFunctionFlex *mufunc    = new RooGBRFunctionFlex( "mufunc",    "" );
    RooGBRFunctionFlex *n1func    = new RooGBRFunctionFlex( "n1func",    "" );
    RooGBRFunctionFlex *n2func    = new RooGBRFunctionFlex( "n2func",    "" );

    // JBS: define mapping of input variables to non-parametric functions
    //   ( in this case trivial since all 4 functions depend on the same inputs,
    //     but this is not a requirement )
    RooGBRTargetFlex *sigmamap = new RooGBRTargetFlex( "sigmamap", "", *sigmafunc, sigmavar, condvars );
    RooGBRTargetFlex *mumap    = new RooGBRTargetFlex( "mumap",    "", *mufunc,    muvar,    condvars );
    RooGBRTargetFlex *n1map    = new RooGBRTargetFlex( "n1map",    "", *n1func,    n1var,    condvars );
    RooGBRTargetFlex *n2map    = new RooGBRTargetFlex( "n2map",    "", *n2func,    n2var,    condvars );

    // JBS: define list of mapped functions to regress
    RooArgList tgts;
    tgts.add( *sigmamap );
    tgts.add( *mumap );
    tgts.add( *n1map );
    tgts.add( *n2map );

    // JBS: define transformations corresponding to parameter bounds for non-parametric outputs  
    RooRealConstraint sigmalim( "sigmalim",  "",  *sigmamap,  par_sigma_DownLimit, par_sigma_UpLimit );
    RooRealConstraint mulim(    "mulim",     "",  *mumap,     par_mu_DownLimit, par_mu_UpLimit );
    RooRealConstraint n1lim(    "n1lim",     "",  *n1map,     par_n1_DownLimit, par_n1_UpLimit );
    RooRealConstraint n2lim(    "n2lim",     "",  *n2map,     par_n2_DownLimit, par_n2_UpLimit );

    // JBS: define pdf, which depends on transformed outputs
    //      (and is intended to be treated as a conditional pdf over the
    //       regression inputs in this case)
    // JBS: The actual pdf below is a double crystal ball,
    //      with crossover points alpha_1 and alpha_2 set constant, but all other
    //      parameters regressed
    RooDoubleCBFast pdf(
        "pdf",
        "",
        *targetvar,
        mulim, sigmalim, RooConst( par_alpha1 ), n1lim, RooConst( par_alpha2 ), n2lim
        );
    cout << endl;


    // ========================================
    // Create RooHybridBDTAutoPdf, something slightly more contrived
    // Probably only needed for programming reasons

    cout << "INFO: Creating RooHybridBDTAutoPdf bdtpdfdiff" << endl;

    // JBS: dummy variable
    RooConstVar etermconst( "etermconst", "", 0. );  

    // JBS: dummy variable
    RooRealVar r( "r", "", 1. );
    r.setConstant();

    // JBS: define list of pdfs
    std::vector<RooAbsReal*> vpdf;
    vpdf.push_back( &pdf );

    // JBS: define list of training datasets
    std::vector<RooAbsData*> vdata;
    vdata.push_back( hdata );

    RooHybridBDTAutoPdf trainer( "bdtpdfdiff", "", tgts, etermconst, r, vdata, vpdf );
    cout << endl;


    // ========================================
    // Apply options on trainer

    int ntrees;
    cout << "INFO: Applying options" << endl;
    for(const auto& token : optionTokens){
        vector<string> tagAndValue;
        tokenize( token, tagAndValue, "=" );
        if(tagAndValue.size()!=2){
            cout << "ERROR: option " << token << " cannot be processed. Should be of the form tag=value.\n";
            continue;
            }
        string tag = tagAndValue[0];
        string value = tagAndValue[1];

        if(tag=="MinEvents") {
            double minEvents = 0;
            fromString(minEvents, value);
            std::vector<double> vMinEvents;
            vMinEvents.push_back(minEvents);
            trainer.SetMinWeights(vMinEvents);
            cout << "INFO: MinEvents = " << minEvents << "\n";
            }
        else if(tag=="Shrinkage") {
            float shrink = 0.;
            fromString(shrink, value);
            trainer.SetShrinkage(shrink);
            cout << "INFO: Shrinkage = " << shrink << "\n";
            }
        else if(tag=="MinSignificance") {
            float sig = 0.;
            fromString(sig, value);
            trainer.SetMinCutSignificance(sig);
            cout << "INFO: MinSignificance = " << sig << "\n";
            }
        else if(tag=="TransitionQuantile") {
            float trans = 0.;
            fromString(trans, value);
            trainer.SetTransitionQuantile(trans);
            cout << "INFO: TransitionQuantile = " << trans << "\n";
            }
        else if(tag=="NTrees") {
            fromString(ntrees, value);
            cout << "INFO: NTrees = " << ntrees << "\n";
            }
        else if(tag=="EventWeight") {
            // Only here so no error message is displayed; option is already applied before
            }
        else{
            cout << "ERROR: SemiparametricGBRMaker::prepareTraining(): Unknown option " << tag << "\n";
            cout << "ERROR: Possibilities are: MinEvents, Shrinkage, MinSignificance, TransitionQuantile, NTrees, EventWeight\n";
            }
        }
    cout << endl;


    // ========================================
    // Run trainer

    cout << "INFO: start training BDT for central value estimation" << endl;
    trainer.TrainForest(ntrees);
    cout << "INFO: end of training BDT" << endl << endl << endl;


    cout << "INFO: Writing results to output " << outputFile_ << endl;


    TFile* outputRootFp = TFile::Open( outputFile_.c_str(), "UPDATE" );
    outputRootFp->WriteObject( &variables_, "varlist" );

    GBRForestD* forestmu;
    forestmu    = new GBRForestD( *mumap->Forest() );
    outputRootFp->WriteObject( forestmu,    "muCB" );

    GBRForestD* forestsigma;
    forestsigma = new GBRForestD( *sigmamap->Forest() );
    outputRootFp->WriteObject( forestsigma, "sigmaCB" );

    GBRForestD* forestn1;
    forestn1    = new GBRForestD( *n1map->Forest() );
    outputRootFp->WriteObject( forestn1,    "n1CB" );

    GBRForestD* forestn2;
    forestn2    = new GBRForestD( *n2map->Forest() );
    outputRootFp->WriteObject( forestn2,    "n2CB" );

    outputRootFp->Close();


    RooWorkspace ws( "w" );
    ws.import( pdf );
    ws.writeToFile( outputFile_.c_str(), false );

}

