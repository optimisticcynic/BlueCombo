// Standard Library
#include <stdexcept>
#include <iostream>
#include <sstream>
#include "BlueForPhistar/Blue.h"


#include "TSystem.h"
#include <set>
#include <TH1.h>
#include "TMath.h"
#include "TMatrixDSparse.h"
#include "TRandom.h"
#include "TMatrix.h"
#include "TVector.h"
#include "TArray.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include <iomanip>
#include "TLatex.h"
#include <TH1D.h>
#include <TFile.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLine.h>
#include <sstream>
#include <string>
#include <iostream>
#include <sstream>


using namespace std;
const double phistarBins[] = {0.000, 0.004, 0.008, 0.012, 0.016, 0.020, 0.024, 0.029, 0.034, 0.039, 0.045, 0.052, 0.057, 0.064, 0.072, 0.081, 0.091, 0.102, 0.114, 0.128, 0.145, 0.165, 0.189, 0.219, 0.258, 0.312, 0.391, 0.524, 0.695, 0.918, 1.153, 1.496, 1.947, 2.522, 3.277};
size_t nphistar = (sizeof (phistarBins) / sizeof (phistarBins[0])) - 1;
const double yBins[] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
size_t ny = (sizeof (yBins) / sizeof (yBins[0])) - 1;
size_t nbins = nphistar;
static const Int_t NumEst = 68;

void HeapUser(double**& Array, Int_t Dimension) {
    //    Array = (double **) malloc(Dimension * sizeof (double *));
    //    for (Int_t YaxisBin = 0; YaxisBin < Dimension; YaxisBin++) {
    //        Array[YaxisBin] = (double *) malloc(Dimension * sizeof (double));
    //        for (Int_t XaxisBin = 0; XaxisBin < Dimension; XaxisBin++) {
    //            Array[XaxisBin][YaxisBin]=0;
    //        }
    //    }
    Dimension = NumEst;
    Array = new double*[Dimension];
    for (Int_t YaxisBin = 0; YaxisBin < Dimension; YaxisBin++) {
        Array[YaxisBin] = new double[Dimension];
        for (Int_t XaxisBin = 0; XaxisBin < Dimension; XaxisBin++) {

        }
    }
    // cout << "Sanity check " << Array[69] << endl;

}

//void CorrilationMatrixMaker(double** Array, TMatrixD CovMatrix) {
//    for (size_t Ybin = 0; Ybin < nphistar * 2; Ybin++) {
//        for (size_t Xbin = 0; Xbin < nphistar * 2; Xbin++) {
//            if(CovMatrix(Xbin, Xbin)!=0&&CovMatrix(Ybin, Ybin))
//            Array[Xbin][Ybin] = CovMatrix(Xbin, Ybin) / sqrt(CovMatrix(Xbin, Xbin) * CovMatrix(Ybin, Ybin))*(1-.000001*(Xbin!=Ybin));
//            else{
//                if(Xbin==Ybin) Array[Xbin][Ybin]=1;
//                else Array[Xbin][Ybin]=0;
//            }
//            //Array[Ybin][Xbin] = Array[Xbin][Ybin];
//        }
//    }
//}

void CorrilationMatrixMaker(double Array[][NumEst], TMatrixD CovMatrix, bool IsQCD = false) {
    for (size_t Ybin = 0; Ybin < nphistar * (2 - IsQCD); Ybin++) {
        for (size_t Xbin = 0; Xbin < nphistar * (2 - IsQCD); Xbin++) {
            Array[Xbin][Ybin] = CovMatrix(Xbin, Ybin) / sqrt(CovMatrix(Xbin, Xbin) * CovMatrix(Ybin, Ybin))*(1 - .00000001 * (Xbin != Ybin));
            Array[Ybin][Xbin] = Array[Xbin][Ybin];
        }
    }
    if (IsQCD)for (size_t Ybin = 0; Ybin < nphistar * (2 - IsQCD); Ybin++) {
            for (size_t Xbin = Ybin; Xbin < nphistar * (2 - IsQCD); Xbin++) {
                Array[Xbin][Ybin] = 0;
                Array[Ybin][Xbin] = Array[Xbin][Ybin];
            }
            Array[Ybin][Ybin] = 1;
        }
}

void symmetricMaker(double** Array) {
    for (size_t Ybin = 0; Ybin < nphistar * 2; Ybin++) {
        for (size_t Xbin = Ybin; Xbin < nphistar * 2; Xbin++) {
            Array[Ybin][Xbin] = Array[Xbin][Ybin];
        }
    }
}

void symmetricMaker(double Array[][NumEst]) {
    for (size_t Ybin = 0; Ybin < nphistar * 2; Ybin++) {
        for (size_t Xbin = Ybin; Xbin < nphistar * 2; Xbin++) {
            Array[Ybin][Xbin] = Array[Xbin][Ybin];
        }
    }
}

void ElectPlusMuon1D(bool DoNorm = false) {


    static const Int_t NumUnc = 8;
    TString NamUnc[NumUnc];
    Int_t NumObs = nphistar;

    TString NamEst[NumEst];
    TString NamObs[NumObs];
    Int_t IWhichObs[NumEst];

    for (size_t PhiStarBin = 0; PhiStarBin < NumObs; PhiStarBin++) {
        TString ElectronBinName;
        ElectronBinName.Format("Electron_Bin_%d02", PhiStarBin);
        TString MuonBinName;
        MuonBinName.Format("Muon_Bin_%d02", PhiStarBin);
        TString BinName;
        BinName.Format("Bin_%d02", PhiStarBin);
        NamEst[PhiStarBin] = ElectronBinName;
        NamEst[PhiStarBin + NumObs] = MuonBinName;
        NamObs[PhiStarBin] = BinName;
        IWhichObs[PhiStarBin] = PhiStarBin;
        IWhichObs[PhiStarBin + NumObs] = PhiStarBin;
    }
    string OrignialFileLocation;
    if (DoNorm)OrignialFileLocation = "~/work/HomeWork/Phistar/CombineElectWithMu/EPlusMuFileV2/Comb_ForBlue_Norm_1D_Born.root";
    else OrignialFileLocation = "~/work/HomeWork/Phistar/CombineElectWithMu/EPlusMuFileV2/Comb_ForBlue_Abs_1D_Born.root";

    TFile Original(OrignialFileLocation.c_str());
    Original.cd();

    TGraphAsymmErrors* FullPlot = (TGraphAsymmErrors*) Original.Get("Nominal");



    TMatrixD* CovM_tot = (TMatrixD*) Original.Get("CovM_tot");
    if (!CovM_tot) {
        cout << "can't find covM_tot" << endl;
        return;
    }
    TMatrixD* CovM_stat = (TMatrixD*) Original.Get("CovM_stat");
    TMatrixD* CovM_mcstat = (TMatrixD*) Original.Get("Cov_mcstat");
    TMatrixD* CovM_eff = (TMatrixD*) Original.Get("CovM_eff");
    TMatrixD* CovM_bg_tt = (TMatrixD*) Original.Get("CovM_bg_tt");
    TMatrixD* CovM_bg_di = (TMatrixD*) Original.Get("CovM_bg_di");
    TMatrixD* CovM_bg_qcd = (TMatrixD*) Original.Get("CovM_bg_qcd");
    TMatrixD* CovM_pileup = (TMatrixD*) Original.Get("CovM_pileup");
    TMatrixD* CovM_pt = (TMatrixD*) Original.Get("CovM_pt");
    TMatrixD* CovM_lumi = (TMatrixD*) Original.Get("CovM_lumi");
    NamUnc[0].Format("stat");
    NamUnc[1].Format("mcstat");
    NamUnc[2].Format("eff");
    NamUnc[3].Format("bg_tt");
    NamUnc[4].Format("bg_di");
    NamUnc[5].Format("bg_qcd");
    NamUnc[6].Format("pileup");
    NamUnc[7].Format("pt");
    vector<TMatrixD> AllCovs;
    AllCovs.push_back(*CovM_stat);
    AllCovs.push_back(*CovM_mcstat);
    AllCovs.push_back(*CovM_eff);
    AllCovs.push_back(*CovM_bg_tt);
    AllCovs.push_back(*CovM_bg_di);
    AllCovs.push_back(*CovM_bg_qcd);
    AllCovs.push_back(*CovM_pileup);
    AllCovs.push_back(*CovM_pt);
    AllCovs.push_back(*CovM_lumi);

    //double** CorM_statArray;
    //double CorM_statArray[68][68];
    auto CorM_statArray = new double[NumEst][NumEst]();
    auto CorM_mcstatArray = new double[NumEst][NumEst]();
    auto CorM_effArray = new double[NumEst][NumEst]();
    auto CorM_bg_ttArray = new double[NumEst][NumEst]();
    auto CorM_bg_diArray = new double[NumEst][NumEst]();
    auto CorM_bg_qcdArray = new double[NumEst][NumEst]();
    auto CorM_pileupArray = new double[NumEst][NumEst]();
    auto CorM_ptArray = new double[NumEst][NumEst]();




    //    double** CorM_statArray[68][68];
    //    double** CorM_mcstatArray;
    //    double** CorM_effArray;
    //    double** CorM_bg_ttArray;
    //    double** CorM_bg_diArray;
    //    double** CorM_bg_qcdArray;
    //    double** CorM_pileupArray;
    //    double** CorM_ptArray;
    //    double** CorM_lumiArray;
    //    HeapUser(CorM_totArray, NumEst);



    // HeapUser(CorM_statArray, NumEst);
    //HeapUser(CorM_mcstatArray, NumEst);
    //    HeapUser(CorM_effArray, NumEst);
    //    HeapUser(CorM_bg_ttArray, NumEst);
    //    HeapUser(CorM_bg_diArray, NumEst);
    //    HeapUser(CorM_bg_qcdArray, NumEst);
    //    HeapUser(CorM_pileupArray, NumEst);
    //    HeapUser(CorM_ptArray, NumEst);
    //    HeapUser(CorM_lumiArray, NumEst);

    CorrilationMatrixMaker(CorM_statArray, (*CovM_stat));
    CorrilationMatrixMaker(CorM_mcstatArray, (*CovM_mcstat));
    CorrilationMatrixMaker(CorM_effArray, (*CovM_eff));
    CorrilationMatrixMaker(CorM_bg_ttArray, (*CovM_bg_tt));
    CorrilationMatrixMaker(CorM_bg_diArray, (*CovM_bg_di));
    CorrilationMatrixMaker(CorM_bg_qcdArray, (*CovM_bg_qcd), true);
    CorrilationMatrixMaker(CorM_pileupArray, (*CovM_pileup));
    CorrilationMatrixMaker(CorM_ptArray, (*CovM_pt));
    //NoLumi

    Blue *myBlue = new Blue(NumEst, NumUnc, NumObs, &IWhichObs[0]);
    myBlue->PrintStatus();
    // Fill names
    myBlue->FillNamEst(NamEst);
    myBlue->FillNamUnc(NamUnc);
    myBlue->FillNamObs(NamObs);
    static const Int_t LenXEst = NumEst * (NumUnc + 1);
    Double_t XEst[LenXEst];

    for (size_t PhiStarBin = 0; PhiStarBin < NumEst; PhiStarBin++) {
        double x, y;
        FullPlot->GetPoint(PhiStarBin, x, y);
        XEst[PhiStarBin * (NumUnc + 1)] = y;
        XEst[PhiStarBin * (NumUnc + 1) + 1] = sqrt((*CovM_stat)(PhiStarBin, PhiStarBin));
        XEst[PhiStarBin * (NumUnc + 1) + 2] = sqrt((*CovM_mcstat)(PhiStarBin, PhiStarBin));
        XEst[PhiStarBin * (NumUnc + 1) + 3] = sqrt((*CovM_eff)(PhiStarBin, PhiStarBin));
        XEst[PhiStarBin * (NumUnc + 1) + 4] = sqrt((*CovM_bg_tt)(PhiStarBin, PhiStarBin));
        XEst[PhiStarBin * (NumUnc + 1) + 5] = sqrt((*CovM_bg_di)(PhiStarBin, PhiStarBin));
        XEst[PhiStarBin * (NumUnc + 1) + 6] = sqrt((*CovM_bg_qcd)(PhiStarBin, PhiStarBin));
        XEst[PhiStarBin * (NumUnc + 1) + 7] = sqrt((*CovM_pileup)(PhiStarBin, PhiStarBin));
        XEst[PhiStarBin * (NumUnc + 1) + 8] = sqrt((*CovM_pt)(PhiStarBin, PhiStarBin));
    }

    Int_t ind = 0;
    for (Int_t i = 0; i < NumEst; i++) {
        myBlue->FillEst(i, &XEst[ind]);
        ind = ind + NumUnc + 1;
    }


    //    for (Int_t BinX = 0; BinX < NumEst; BinX++) {
    //        for (Int_t BinY = 40; BinY < NumEst; BinY++) {
    //            printf("%.3f\t", CorM_statArray[BinX][BinY]);
    //        }
    //        cout << endl;
    //    }

    symmetricMaker(CorM_statArray);
    myBlue->FillCor(0, CorM_statArray[0]);




    myBlue->FillCor(1, CorM_mcstatArray[0]);
    myBlue->FillCor(2, CorM_effArray[0]);
    myBlue->FillCor(3, CorM_bg_ttArray[0]);
    myBlue->FillCor(4, CorM_bg_diArray[0]);
    myBlue->FillCor(5, CorM_bg_qcdArray[0]);
    myBlue->FillCor(6, CorM_pileupArray[0]);
    myBlue->FillCor(7, CorM_ptArray[0]);

    myBlue->FixInp();
    myBlue->Solve();
    double Results[NumObs][NumUnc + 1];
    TMatrixD* CovarianceResults = new TMatrixD(NumObs, NumObs);
    myBlue->GetResult(Results[0]);
    myBlue->GetCovRes(CovarianceResults);


    TGraphAsymmErrors ElectronPlot(NumObs);
    TGraphAsymmErrors MuonPlot(NumObs);
    TGraphAsymmErrors BlueCombGraph(NumObs);




    double ElectronErrorSquared[NumObs];
    double MuonErrorSquared[NumObs];
    double CombinedErrorSquared[NumObs];

    for (size_t BinNumber = 0; BinNumber < NumObs; BinNumber++) {
        ElectronErrorSquared[BinNumber] = AllCovs[0](BinNumber, BinNumber) + AllCovs[1](BinNumber, BinNumber) + AllCovs[2](BinNumber, BinNumber) + AllCovs[3](BinNumber, BinNumber) +
                AllCovs[4](BinNumber, BinNumber) + AllCovs[5](BinNumber, BinNumber) + AllCovs[6](BinNumber, BinNumber) + AllCovs[7](BinNumber, BinNumber);
        size_t muonBin = BinNumber + NumObs;
        MuonErrorSquared[BinNumber] = AllCovs[0](muonBin, muonBin) + AllCovs[1](muonBin, muonBin) + AllCovs[2](muonBin, muonBin) + AllCovs[3](muonBin, muonBin) +
                AllCovs[4](muonBin, muonBin) + AllCovs[5](muonBin, muonBin) + AllCovs[6](muonBin, muonBin) + AllCovs[7](muonBin, muonBin);
        CombinedErrorSquared[BinNumber] = (*CovarianceResults)(BinNumber, BinNumber);
    }

    for (size_t BinNumber = 0; BinNumber < NumObs; BinNumber++) {
        double x, y;
        FullPlot->GetPoint(BinNumber, x, y);
        double RealX = (phistarBins[BinNumber] + phistarBins[BinNumber + 1]) / 2;
        double XError = -(phistarBins[BinNumber] - phistarBins[BinNumber + 1]) / 2;
        double YerrorElect = sqrt(ElectronErrorSquared[BinNumber] + y * y * .026 * .026 * (!DoNorm));
        ElectronPlot.SetPoint(BinNumber, RealX, y);
        ElectronPlot.SetPointError(BinNumber, XError, XError, YerrorElect, YerrorElect);

        FullPlot->GetPoint(BinNumber + NumObs, x, y);
        double YerrorMuon = sqrt(MuonErrorSquared[BinNumber] + y * y * .026 * .026 * (!DoNorm));
        MuonPlot.SetPoint(BinNumber, RealX, y);
        MuonPlot.SetPointError(BinNumber, XError, XError, YerrorMuon, YerrorMuon);

        BlueCombGraph.SetPoint(BinNumber, RealX, Results[BinNumber][0]);
        double CombError = sqrt(CombinedErrorSquared[BinNumber] + Results[BinNumber][0] * Results[BinNumber][0]*.026 * .026 * (!DoNorm));
        BlueCombGraph.SetPointError(BinNumber, XError, XError, CombError, CombError);

    }

    string filename;
    if (DoNorm)filename = "Results/Comb_Norm_UsingBlue1D.root";
    else filename = "Results/Comb_Abs_UsingBlue1D.root";

    TFile ResultFile(filename.c_str(), "recreate");
    ResultFile.cd();
    ElectronPlot.Write("Electron");
    MuonPlot.Write("Muons");
    BlueCombGraph.Write("h_Comb");
    CovarianceResults->Write("TotalCovarianceMatrix");

    ResultFile.Write();
    delete myBlue;
    delete CovarianceResults;
}

TGraphAsymmErrors * CreateRatio(TGraphAsymmErrors* graph, TGraphAsymmErrors* graphmc, bool isData = false) {
    double x, y, errorl, errorh, xmc, ymc, errorlmc, errorhmc;
    TGraphAsymmErrors* g_ratio = new TGraphAsymmErrors(nbins);
    for (size_t ibin = 0; ibin < nbins; ibin++) {
        graph->GetPoint(ibin, x, y);
        if (!isData) {
            graphmc->GetPoint(ibin, xmc, ymc);
            errorlmc = graphmc->GetErrorYlow(ibin);
            errorhmc = graphmc->GetErrorYhigh(ibin);
            g_ratio->SetPoint(ibin, x, ymc / y);
            g_ratio->SetPointError(ibin, 0, 0, errorlmc / y, errorhmc / y);
        } else {
            errorl = graph->GetErrorYlow(ibin);
            errorh = graph->GetErrorYhigh(ibin);
            g_ratio->SetPoint(ibin, x, 1);
            double DataError = -(phistarBins[ibin] - phistarBins[ibin + 1]) / 2;
            g_ratio->SetPointError(ibin, DataError, DataError, errorl / y, errorh / y);
        }
    }
    return g_ratio;
}

void Plotter(bool DoNorm = true) {

    string FileName;
    if (DoNorm)FileName = "Results/Comb_Norm_UsingBlue1D.root";
    else FileName = "Results/Comb_Abs_UsingBlue1D.root";
    TFile OpenFileName(FileName.c_str());

    TGraphAsymmErrors* ElectronFull = (TGraphAsymmErrors*) OpenFileName.Get("Electron");
    TGraphAsymmErrors* MuonFull = (TGraphAsymmErrors*) OpenFileName.Get("Muons");
    TGraphAsymmErrors* TotalFull = (TGraphAsymmErrors*) OpenFileName.Get("h_Comb");


    TGraphAsymmErrors* ElectronRatioFull = CreateRatio(TotalFull, ElectronFull);
    TGraphAsymmErrors * MuonRatioFull = CreateRatio(TotalFull, MuonFull);
    TGraphAsymmErrors * DataRatioFull = CreateRatio(TotalFull, TotalFull, true);

    TCanvas* FinalPhiRatio = new TCanvas("Powheg", "PowhegPlot", 800, 900);
    FinalPhiRatio->SetLogx();
    gStyle->SetOptStat("");
    FinalPhiRatio->cd();
    //FinalPhiRatio->SetLeftMargin(0);

    DataRatioFull->SetFillColor(kGray);
    DataRatioFull->SetTitle("");
    DataRatioFull->GetYaxis()->SetTitle("Data/Blue");
    DataRatioFull->GetYaxis()->CenterTitle();
    DataRatioFull->Draw("AE2");


    MuonRatioFull->SetMarkerColor(kRed);
    MuonRatioFull->SetLineColor(kRed);
    MuonRatioFull->SetMarkerStyle(21);
    MuonRatioFull->Draw("PE same");
    ElectronRatioFull->SetMarkerColor(kBlue);
    ElectronRatioFull->SetLineColor(kBlue);
    ElectronRatioFull->SetMarkerStyle(23);
    ElectronRatioFull->Draw("PEsame");

    TLegend* leg2 = new TLegend(0.1, 0.9, 0.9, 0.97);
    leg2->SetNColumns(3);
    leg2->SetFillStyle(0);
    //leg2->SetBorderSize(1);
    leg2->SetLineWidth(1);
    leg2->SetTextFont(22);
    double TextSize = 0.02;
    leg2->SetTextSize(TextSize);
    leg2->AddEntry(DataRatioFull, "Blue Combined Z #rightarrow ll", "F");
    leg2->AddEntry(MuonRatioFull, "2012 data Z #rightarrow #mu#mu", "P");
    leg2->AddEntry(ElectronRatioFull, "2012 Data Z #rightarrow ee", "P");
    leg2->Draw();


    string PlotName;
    if (DoNorm)PlotName = "Plots/OneDNorm";
    else PlotName = "Plots/OneDAbs";
    string Withtype = PlotName + ".pdf";
    FinalPhiRatio->Print(Withtype.c_str());

}

int main(int argc, char * argv[]) {
    ElectPlusMuon1D(false);
    ElectPlusMuon1D(true);
    Plotter(true);
    Plotter(false);
    return 1;
}