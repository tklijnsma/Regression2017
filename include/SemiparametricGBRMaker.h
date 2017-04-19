// Minimalistic GRBMaker for semiparametric regression
// Based on work by Jean-Baptiste Sauvan and Josh Bendavid

#ifndef SEMIPARAMETRICGBRMAKER_H
#define SEMIPARAMETRICGBRMAKER_H

// STD
#include <string>
#include <vector>

class TTree;
class TChain;
class TFile;
class RooHybridBDTAutoPdf;
class RooDoubleCBFast;
// class GBRTrainer;
class GBRForestD;

class SemiparametricGBRMaker
{
    public:
        SemiparametricGBRMaker();
        ~SemiparametricGBRMaker();

        void readParameterFile( const std::string& parFileName );
        void run();
        void close();


    private:

        // Parameter file variables
        std::string par_Name ;
        std::string par_OutputDirectory ;
        std::string par_InputFiles ;
        std::string par_Tree ;
        std::string par_Options ;
        std::string par_Variables ;
        std::string par_Target ;
        std::string par_Cut ;
        float par_mu_DownLimit    ;
        float par_mu_UpLimit      ;
        float par_sigma_DownLimit ;
        float par_sigma_UpLimit   ;
        float par_n1_DownLimit    ;
        float par_n1_UpLimit      ;
        float par_n2_DownLimit    ;
        float par_n2_UpLimit      ;
        float par_alpha1          ;
        float par_alpha2          ;

        // Variables derived from parameter file
        TChain*                   chain_;
        std::vector<TFile*>       rootFps_;
        std::vector<TTree*>       trees_;
        std::vector<std::string>  variables_;
        std::string               outputFile_;

};

#endif