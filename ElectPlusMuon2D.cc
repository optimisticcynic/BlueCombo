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
#include <TH2.h>


using namespace std;
const double phistarBins[] = {0.000, 0.004, 0.008, 0.012, 0.016, 0.020, 0.024, 0.029, 0.034, 0.039, 0.045, 0.051, 0.057, 0.064, 0.072, 0.081, 0.091, 0.102, 0.114, 0.128, 0.145, 0.165, 0.189, 0.219, 0.258, 0.312, 0.391, 0.524, 0.695, 0.918, 1.153, 1.496, 1.947, 2.522, 3.277};
size_t nphistar = (sizeof (phistarBins) / sizeof (phistarBins[0])) - 1;
const double yBins[] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
size_t ny = (sizeof (yBins) / sizeof (yBins[0])) - 1;
size_t nbins = nphistar*ny;
static const Int_t NumEst = 408;

void CorrilationMatrixMaker(double Array[][NumEst], TMatrixD CovMatrix, bool IsQCD = false) {
    for (size_t Ybin = 0; Ybin < nbins * (2 - IsQCD); Ybin++) {
        for (size_t Xbin = 0; Xbin < nbins * (2 - IsQCD); Xbin++) {
            Array[Xbin][Ybin] = CovMatrix(Xbin, Ybin) / sqrt((CovMatrix(Xbin, Xbin) + 1e-15) * (CovMatrix(Ybin, Ybin) + 1e-15))*(1 - .00000001 * (Xbin != Ybin));
            Array[Ybin][Xbin] = Array[Xbin][Ybin];
        }
    }
    if (IsQCD)for (size_t Ybin = 0; Ybin < nbins * (2 - IsQCD); Ybin++) {
            for (size_t Xbin = Ybin; Xbin < nbins * (2 - IsQCD); Xbin++) {
                Array[Xbin][Ybin] = 0;
                Array[Ybin][Xbin] = Array[Xbin][Ybin];
            }
            Array[Ybin][Ybin] = 1;
        }
}

void symmetricMaker(double** Array) {
    for (size_t Ybin = 0; Ybin < nbins * 2; Ybin++) {
        for (size_t Xbin = Ybin; Xbin < nbins * 2; Xbin++) {
            Array[Ybin][Xbin] = Array[Xbin][Ybin];
        }
    }
}

void symmetricMaker(double Array[][408]) {
    for (size_t Ybin = 0; Ybin < nbins * 2; Ybin++) {
        for (size_t Xbin = Ybin; Xbin < nbins * 2; Xbin++) {
            Array[Ybin][Xbin] = Array[Xbin][Ybin];
        }
    }
}

void ElectPlusMuon2D(bool DoNorm = false, size_t RemoveCorrilation = -1) {


    static const Int_t NumUnc = 8;
    TString NamUnc[NumUnc];
    Int_t NumObs = nbins;

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
        NamEst[PhiStarBin + nbins] = MuonBinName;
        NamObs[PhiStarBin] = BinName;
        IWhichObs[PhiStarBin] = PhiStarBin;
        IWhichObs[PhiStarBin + nbins] = PhiStarBin;
    }
    string OrignialFileLocation;
    if (DoNorm)OrignialFileLocation = "~/work/HomeWork/Phistar/CombineElectWithMu/EPlusMuFileV2/Comb_ForBlue_Norm_2D_Born.root";
    else OrignialFileLocation = "~/work/HomeWork/Phistar/CombineElectWithMu/EPlusMuFileV2/Comb_ForBlue_Abs_2D_Born.root";

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
    //    (*CovM_stat) *= 0;
    for (size_t binx = 0; binx < NumEst; binx++)
        for (size_t biny = 0; biny < NumEst; biny++) {
            if (binx == biny)continue;
            if (RemoveCorrilation == 0)(*CovM_stat)(binx, biny) = 0;
            if (RemoveCorrilation == 1)(*CovM_mcstat)(binx, biny) = 0;
            if (RemoveCorrilation == 2)(*CovM_eff)(binx, biny) = 0;
            if (RemoveCorrilation == 3)(*CovM_bg_tt)(binx, biny) = 0;
            if (RemoveCorrilation == 4)(*CovM_bg_di)(binx, biny) = 0;
            if (RemoveCorrilation == 5)(*CovM_bg_qcd)(binx, biny) = 0;
            if (RemoveCorrilation == 6)(*CovM_pileup)(binx, biny) = 0;
            if (RemoveCorrilation == 7)(*CovM_pt)(binx, biny) = 0;
        }

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

    auto CorM_statArray = new double[NumEst][NumEst]();
    auto CorM_mcstatArray = new double[NumEst][NumEst]();
    auto CorM_effArray = new double[NumEst][NumEst]();
    auto CorM_bg_ttArray = new double[NumEst][NumEst]();
    auto CorM_bg_diArray = new double[NumEst][NumEst]();
    auto CorM_bg_qcdArray = new double[NumEst][NumEst]();
    auto CorM_pileupArray = new double[NumEst][NumEst]();
    auto CorM_ptArray = new double[NumEst][NumEst]();





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
    TMatrixD MatrixResults(NumObs, NumUnc + 1);
    TMatrixD* Weights = new TMatrixD(NumEst, NumObs);
    myBlue->GetWeight(Weights);

    myBlue->GetResult(Results[0]);
    myBlue->GetCovRes(CovarianceResults);
    myBlue->GetResult(&MatrixResults);


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

    for (size_t YBin = 0; YBin < ny; YBin++) {
        for (size_t PhiStarBin = 0; PhiStarBin < nphistar; PhiStarBin++) {
            double x, y;
            size_t BinNumber = nphistar * YBin + PhiStarBin;
            FullPlot->GetPoint(BinNumber, x, y);
            double RealX = (phistarBins[PhiStarBin] + phistarBins[PhiStarBin + 1]) / 2 + YBin * 3.277;
            double XError = -(phistarBins[PhiStarBin] - phistarBins[PhiStarBin + 1]) / 2;
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
    }






    string filename;
    if (RemoveCorrilation == -1) {
        if (DoNorm)filename = "Results/Comb_Norm_UsingBlue2D.root";
        else filename = "Results/Comb_Abs_UsingBlue2Dffcor.root"; //Weird thing getting rid of it
    } else if (RemoveCorrilation == 2) {
        if (DoNorm)filename = "Results/Comb_Norm_UsingBlue2DRemovedEffcor.root";
        else filename = "Results/Comb_Abs_UsingBlue2DRemovedEffcor.root";
    } else {
        cout << "To do doing an option that I shouldn't be yet so haven't asdfasdf " << endl;
        return;
    }
    TFile ResultFile(filename.c_str(), "recreate");
    ResultFile.cd();
    ElectronPlot.Write("Electron");
    MuonPlot.Write("Muons");
    BlueCombGraph.Write("h_Comb");
    CovarianceResults->Write("TotalCovarianceMatrix");
    MatrixResults.Write("ResultsMatrix");
    Weights->Write("Weights");


    TH2D SanityCheckinitial("corrilationHistoInitial", "corrilationHistoInitial", NumEst, 0, NumEst, NumEst, 0, NumEst);
    TH2D SanityCheckEnd("corrilationHistoEnd", "corrilationHistoEnd", NumObs, 0, NumObs, NumObs, 0, NumObs);





    TMatrixD OldCovarianMatrixTotal = AllCovs[7];

    OldCovarianMatrixTotal = AllCovs[0] + AllCovs[1] + AllCovs[2] + AllCovs[3] + AllCovs[4] + AllCovs[5] + AllCovs[6];

    for (size_t binx = 0; binx < NumObs * 2; binx++)
        for (size_t biny = 0; biny < NumObs * 2; biny++) {
            SanityCheckinitial.SetBinContent(binx + 1, biny + 1, OldCovarianMatrixTotal(binx, biny) / sqrt(OldCovarianMatrixTotal(binx, binx) * OldCovarianMatrixTotal(biny, biny)));
        }


    for (size_t binx = 0; binx < NumObs; binx++)
        for (size_t biny = 0; biny < NumObs; biny++) {
            double xx, xy, yx, yy;
            BlueCombGraph.GetPoint(binx, xx, xy);
            BlueCombGraph.GetPoint(biny, yx, yy);
            SanityCheckinitial.SetBinContent(binx + 1, biny + 1, OldCovarianMatrixTotal(binx, biny) / sqrt(OldCovarianMatrixTotal(binx, binx) * OldCovarianMatrixTotal(biny, biny)));
            SanityCheckEnd.SetBinContent(binx + 1, biny + 1, (*CovarianceResults)(binx, biny) / sqrt((*CovarianceResults)(binx, binx)* (*CovarianceResults)(biny, biny)));
        }
    TH1D ChiSquared("ChiSquared", "ChiSquared", 1, 0, 1);
    ChiSquared.SetBinContent(1, myBlue->GetChiq());
    ChiSquared.Write();
    //SanityCheck.Write();
    ResultFile.Write();
    delete myBlue;
    delete CovarianceResults;


}

vector<TGraphAsymmErrors*> SplitGraph(TGraphAsymmErrors* graph, bool isData = false) {
    vector<TGraphAsymmErrors*> v;
    for (uint i = 0; i < ny; i++) {
        TGraphAsymmErrors* g = new TGraphAsymmErrors(nphistar);
        for (uint j = 0; j < nphistar; j++) {
            int bin = i * (nphistar) + j;
            double x, y;
            graph->GetPoint(bin, x, y);
            g->SetPoint(j, (phistarBins[j] + phistarBins[j + 1]) / 2., y);
            double DataError = -(phistarBins[j] - phistarBins[j + 1]) / 2 * isData;
            g->SetPointError(j, DataError, DataError, graph->GetErrorYlow(bin), graph->GetErrorYhigh(bin));
            //if (doXerrors) g->SetPointError(j, (phistarBins[j + 1] - phistarBins[j]) / 2., (phistarBins[j + 1] - phistarBins[j]) / 2., graph->GetErrorYlow(bin), graph->GetErrorYhigh(bin));
            //cout << i << " " << j << " " << bin << " " << (phistarBins[j] + phistarBins[j + 1]) / 2. << " " << (phistarBins[j + 1] - phistarBins[j]) / 2. << " " << y << " " << graph->GetErrorYlow(bin) << " " << graph->GetErrorYhigh(bin) << endl;
        }
        v.push_back(g);
    }
    return v;
}

TGraphAsymmErrors * CreateRatio(TGraphAsymmErrors* graph, TGraphAsymmErrors* graphmc, bool isData) {
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
    if (DoNorm)FileName = "Results/Comb_Norm_UsingBlue2D.root";
    else FileName = "Results/Comb_Abs_UsingBlue2D.root";
    TFile OpenFileName(FileName.c_str());

    TGraphAsymmErrors* ElectronFull = (TGraphAsymmErrors*) OpenFileName.Get("Electron");
    TGraphAsymmErrors* MuonFull = (TGraphAsymmErrors*) OpenFileName.Get("Muons");
    TGraphAsymmErrors* TotalFull = (TGraphAsymmErrors*) OpenFileName.Get("h_Comb");


    TGraphAsymmErrors* ElectronRatioFull = CreateRatio(TotalFull, ElectronFull, false);
    TGraphAsymmErrors * MuonRatioFull = CreateRatio(TotalFull, MuonFull, false);
    TGraphAsymmErrors * DataRatioFull = CreateRatio(TotalFull, TotalFull, true);

    vector<TGraphAsymmErrors*> ElectronCollection = SplitGraph(ElectronRatioFull);
    vector<TGraphAsymmErrors*> MuonCollection = SplitGraph(MuonRatioFull);
    vector<TGraphAsymmErrors*> r_data = SplitGraph(DataRatioFull, true);




    TCanvas* FinalPhiRatio = new TCanvas("FinalPhiRatio", "FinalPhiRatio", 800, 900);
    FinalPhiRatio->SetBottomMargin(0.1);
    FinalPhiRatio->SetFillColor(0);

    //TPad* Info = new TPad("p1", "p1", 0, .03, 1, .9);
    //Info->Draw();
    //Info->cd();
    //Info->Divide(1, 6, 0, 0);

    vector<TPad*> Info;
    Info.push_back((new TPad("p1", "p1", 0, .9 - .87 / 6, 1, .9)));
    Info.push_back((new TPad("p2", "p2", 0, .9 - .87 / 6 * 2, 1, .9 - .87 / 6)));
    Info.push_back((new TPad("p3", "p3", 0, .9 - .87 / 6 * 3, 1, .9 - .87 / 6 * 2)));
    Info.push_back((new TPad("p4", "p4", 0, .9 - .87 / 6 * 4, 1, .9 - .87 / 6 * 3)));
    Info.push_back((new TPad("p5", "p5", 0, .9 - .87 / 6 * 5, 1, .9 - .87 / 6 * 4)));
    Info.push_back((new TPad("p6", "p6", 0, 0, 1, .9 - .87 / 6 * 5)));
    double TextSize = 0.02;
    Info[0]->Draw();
    Info[1]->Draw();
    Info[2]->Draw();
    Info[3]->Draw();
    Info[4]->Draw();
    Info[5]->Draw();

    for (uint i = 0; i < ny; i++) {
        //Info->cd(i + 1);

        Info[i]->cd();
        gPad->SetLeftMargin(.0975);
        if (i != 0)gPad->SetTopMargin(0);
        if (i != ny - 1)gPad->SetBottomMargin(0);
        std::ostringstream strs;
        strs << i;
        std::string gPadName = "p" + strs.str();
        if (i == 5)gPad->SetBottomMargin(0.2169);
        gPad->SetLogx(1);

        if (i == 0)r_data[i]->SetTitle("");
        else r_data[i]->SetTitle("");

        r_data[i]->GetYaxis()->SetTitle("");

        r_data[i]->GetYaxis()->CenterTitle();
        r_data[i]->SetFillColor(kGray);
        r_data[i]->GetYaxis()->SetNdivisions(503);

        if (DoNorm) {
            r_data[i]->GetYaxis()->SetRangeUser(0.94, 1.06);
            if (i == 5)r_data[i]->GetYaxis()->SetRangeUser(0.88, 1.12);
        } else {
            r_data[i]->GetYaxis()->SetRangeUser(0.94, 1.06);
            if (i == 4)r_data[i]->GetYaxis()->SetRangeUser(0.88, 1.12);
            if (i == 5) {
                r_data[i]->GetYaxis()->SetRangeUser(0.88, 1.12);
                r_data[i]->GetYaxis()->SetNdivisions(505);
            }
        }
        r_data[i]->GetYaxis()->SetTitleOffset(0.2); //OFFSET
        r_data[i]->GetYaxis()->SetTitleSize(0.2); //Y TITLE SIZE

        r_data[i]->GetXaxis()->SetTitleSize(0);
        r_data[i]->GetXaxis()->SetLabelSize(0);
        r_data[i]->GetXaxis()->SetTitle("");
        r_data[i]->GetXaxis()->CenterTitle();
        if (i == ny - 1) {
            r_data[i]->GetXaxis()->SetTitleOffset(0.9); //OFFSET
            r_data[i]->GetXaxis()->SetTitleSize(TextSize * 6); //X TITLE SIZE
            r_data[i]->GetXaxis()->SetLabelSize(0.12);
            r_data[i]->GetXaxis()->SetLabelOffset(-0.01);
            r_data[i]->GetXaxis()->SetTitle("#phi*");
            r_data[i]->GetXaxis()->CenterTitle();
        }
        r_data[i]->GetXaxis()->SetRangeUser(.015, 3.28);
        r_data[i]->GetXaxis()->SetTickLength(.12);
        r_data[i]->GetYaxis()->SetTickLength(.02);
        r_data[i]->GetYaxis()->SetLabelSize(0.12);
        if (i == 5) r_data[i]->GetYaxis()->SetLabelSize(0.12 * (13.0 / 15.0));
        r_data[i]->Draw("AE2");

        MuonCollection[i]->SetMarkerColor(kRed);
        MuonCollection[i]->SetLineColor(kRed);
        MuonCollection[i]->SetMarkerStyle(21);
        MuonCollection[i]->Draw("PEsame");

        ElectronCollection[i]->SetMarkerColor(kBlue);
        ElectronCollection[i]->SetLineColor(kBlue);
        ElectronCollection[i]->SetMarkerStyle(23);
        ElectronCollection[i]->Draw("PEsame");
        Info[i]->RedrawAxis();
    }

    FinalPhiRatio->cd(0);
    TLatex mark;
    mark.SetTextSize(TextSize);
    mark.SetNDC(kTRUE);
    mark.DrawLatex(.766, .934, "19.7 fb^{-1} (8 TeV)");
    mark.DrawLatex(0.097, .934, "CMS");
    TLatex mark2;
    mark2.SetTextSize(TextSize);
    mark2.SetTextFont(42);
    mark2.SetNDC(kTRUE);
    if (DoNorm) {
        mark2.DrawLatex(.12, .86, "|y| < 0.4");
        mark2.DrawLatex(.12, .73, "0.4 #leq |y| < 0.8");
        mark2.DrawLatex(.12, .58, "0.8 #leq |y| < 1.2");
        mark2.DrawLatex(.12, .44, "1.2 #leq |y| < 1.6");
        mark2.DrawLatex(.12, .30, "1.6 #leq |y| < 2.0");
        mark2.DrawLatex(.12, .15, "2.0 #leq |y| #leq 2.4");
        mark2.SetTextAngle(90);
        mark2.SetTextSize(TextSize * 1.2);
        mark2.DrawLatex(.05, .4, "Separate/Combined ");
    } else if (true) {
        mark2.DrawLatex(.12, .86, "|y| < 0.4");
        mark2.DrawLatex(.12, .73, "0.4 #leq |y| < 0.8");
        mark2.DrawLatex(.12, .58, "0.8 #leq |y| < 1.2");
        mark2.DrawLatex(.12, .44, "1.2 #leq |y| < 1.6");
        mark2.DrawLatex(.12, .30, "1.6 #leq |y| < 2.0");
        mark2.DrawLatex(.12, .15, "2.0 #leq |y| #leq 2.4");
        mark2.SetTextAngle(90);
        mark2.SetTextSize(TextSize * 1.2);
        mark2.DrawLatex(.05, .4, "Separate/Combined");
    }

    FinalPhiRatio->cd(1);
    TLegend* leg2 = new TLegend(0.0975, 0.8865, 0.9, 0.93);
    leg2->SetNColumns(3);
    leg2->SetFillStyle(0);
    //leg2->SetBorderSize(1);
    leg2->SetLineWidth(1);
    leg2->SetTextFont(22);
    leg2->SetTextSize(TextSize);
    leg2->AddEntry(r_data[0], "BLUE Combined", "F");
    leg2->AddEntry(MuonCollection[0], "2012 data Z #rightarrow #mu#mu", "P");
    leg2->AddEntry(ElectronCollection[0], "2012 Data Z #rightarrow ee", "P");
    leg2->Draw();

    string PlotName;
    if (DoNorm)PlotName = "Plots/TwoDNorm";
    else PlotName = "Plots/TwoDAbs";
    string Withtype = PlotName + ".pdf";
    FinalPhiRatio->cd(0);
    FinalPhiRatio->RedrawAxis();
    FinalPhiRatio->Print(Withtype.c_str());
    Withtype = PlotName + ".png";
    FinalPhiRatio->Print(Withtype.c_str());
}

void Combiner() {//HACK HACKEDY HACK HACK HAck so the story is that We don't like the eff causing a weird offset so instead we are now going to combine 2 results
    Int_t NumObs = nphistar*ny;
    //    TFile WithEffCorFile("Results/Comb_Abs_UsingBlue2Dffcor.root", "read");
    //    TH1D* ChiSquared = (TH1D*) WithEffCorFile.Get("ChiSquared");
    //    TH1D myChiSquared("Chi2", "Chi2", 1, 0, 1);
    //    myChiSquared.SetBinContent(1, ChiSquared->GetBinContent(1));
    //    TGraphAsymmErrors* Electron = (TGraphAsymmErrors*) WithEffCorFile.Get("Electron");
    //    if (!Electron) cout << "Missing electrons" << endl;
    //    TGraphAsymmErrors* Muons = (TGraphAsymmErrors*) WithEffCorFile.Get("Muons");
    //    if (!Muons) cout << "Missing Muons" << endl;
    //    TMatrixD* TotalCovarianceMatrix = (TMatrixD*) WithEffCorFile.Get("TotalCovarianceMatrix");
    //    TMatrixD MatrixResults = *((TMatrixD*) WithEffCorFile.Get("ResultsMatrix"));

    //    cout << ChiSquared->GetBinContent(1);
    //    WithEffCorFile.Close();

    std::string OrignialFileLocation = "~/work/HomeWork/Phistar/CombineElectWithMu/EPlusMuFileV2/Comb_ForBlue_Abs_2D_Born.root"; //Grabbing efficency uncertainty shit
    TFile Original(OrignialFileLocation.c_str());
    TMatrixD CovM_eff = *((TMatrixD*) Original.Get("CovM_eff"));
    TMatrixD EffUncertainty(1, NumObs);





    TFile NoEffCorFile("Results/Comb_Abs_UsingBlue2DRemovedEffcor.root", "read");

    TGraphAsymmErrors* h_Comb = (TGraphAsymmErrors*) NoEffCorFile.Get("h_Comb");
    if (!h_Comb) cout << "Missing h_Comb" << endl;
    TMatrixD* TotalCovarianceMatrix = (TMatrixD*) NoEffCorFile.Get("TotalCovarianceMatrix");
    TMatrixD MatrixResults = *((TMatrixD*) NoEffCorFile.Get("ResultsMatrix"));
    TMatrixD Weights = *((TMatrixD*) NoEffCorFile.Get("Weights"));

    TGraphAsymmErrors* Electron = (TGraphAsymmErrors*) NoEffCorFile.Get("Electron");
    if (!Electron) cout << "Missing electrons" << endl;
    TGraphAsymmErrors* Muons = (TGraphAsymmErrors*) NoEffCorFile.Get("Muons");
    if (!Muons) cout << "Missing Muons" << endl;


    for (size_t Obser = 0; Obser < NumObs; Obser++) {
        double uncertaintySquarred = 0;
        //        cout<<"test 1"<<endl;
        for (size_t EstIndex = 0; EstIndex < NumEst; EstIndex++) {
            if (Obser == 0)cout << "Our uncertainty for EstIndex bin" << EstIndex << " is " << sqrt(CovM_eff(EstIndex, EstIndex));
            if (Obser == 0)cout << "    Our Weight for each bin is " << Weights(EstIndex, Obser);
            uncertaintySquarred += Weights(EstIndex, Obser) * Weights(EstIndex, Obser) * CovM_eff(EstIndex, EstIndex);
            if (Obser == 0)cout << "and our current uncertainty is " << sqrt(uncertaintySquarred) << endl;
        }
        if (Obser == 0)cout << "AND OUR UNCERTAINTY IS " << sqrt(uncertaintySquarred);
        EffUncertainty(0, Obser) = sqrt(uncertaintySquarred);
    }



    for (size_t BinY = 0; BinY < h_Comb->GetN(); BinY++) {
        double x1, y1;
        h_Comb->GetPoint(BinY, x1, y1);
        //        MatrixResults(BinY, 0) = y1;
        MatrixResults(BinY, 3) = EffUncertainty(0, BinY);
        cout << "For bin " << BinY << " our Eff uncertainty is " << EffUncertainty(0, BinY) << endl;

        for (size_t BinX = 0; BinX < h_Comb->GetN(); BinX++) {

            double x2, y2;

            h_Comb->GetPoint(BinX, x2, y2);
            (*TotalCovarianceMatrix)(BinY, BinX) = (*TotalCovarianceMatrix)(BinY, BinX) + y1 * y2 * .026 * .026;
            (*TotalCovarianceMatrix)(BinY, BinX) = (*TotalCovarianceMatrix)(BinY, BinX) + EffUncertainty(0, BinX) * EffUncertainty(0, BinY);
        }
    }

    if (!TotalCovarianceMatrix) cout << "Missing TotalCovarianceMatrix" << endl;
    for (size_t binNumber = 0; binNumber < h_Comb->GetN(); binNumber++) {
        double x, y;
        h_Comb->GetPoint(binNumber, x, y);
        h_Comb->SetPointEYhigh(binNumber, sqrt((*TotalCovarianceMatrix)(binNumber, binNumber)));
        h_Comb->SetPointEYlow(binNumber, sqrt((*TotalCovarianceMatrix)(binNumber, binNumber)));
    }
    TFile* Finalthingy = new TFile("Results/Comb_Abs_UsingBlue2DEffFix.root", "recreate");
    Finalthingy->cd();


    Electron->Write("Electron");
    Muons->Write("Muons");
    h_Comb->Write("h_Comb");
    MatrixResults.Write("ResultsMatrix");
    TotalCovarianceMatrix->Write("TotalCovarianceMatrix");



    Finalthingy->Write();
}

int main(int argc, char * argv[]) {
    ElectPlusMuon2D(false);
    ElectPlusMuon2D(false, 2);
    Combiner();

    ElectPlusMuon2D(true);
    Plotter(false);
    Plotter(true);

    return 1;
}