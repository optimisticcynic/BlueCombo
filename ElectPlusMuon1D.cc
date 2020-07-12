// Combines the electron and muon group results using BLUE
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
size_t nbins = nphistar;
static const Int_t NumEst = 68;



void CorrilationMatrixMaker(double Array[][NumEst], TMatrixD CovMatrix, bool IsQCD = false) { //have covmatix blue needs corrilationmatrix
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

void ElectPlusMuon1D(bool DoNorm = false, size_t RemoveCorrilation = -1) {


    static const Int_t NumUnc = 8;
    TString NamUnc[NumUnc];
    Int_t NumObs = nphistar;

    TString NamEst[NumEst];
    TString NamObs[NumObs];
    Int_t IWhichObs[NumEst];

    for (size_t PhiStarBin = 0; PhiStarBin < NumObs; PhiStarBin++) { //creating names for all observables
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
    TMatrixD* CovM_stat = (TMatrixD*) Original.Get("CovM_stat"); //grabbing all the cov matrix
    TMatrixD* CovM_mcstat = (TMatrixD*) Original.Get("Cov_mcstat");
    TMatrixD* CovM_eff = (TMatrixD*) Original.Get("CovM_eff");
    TMatrixD* CovM_bg_tt = (TMatrixD*) Original.Get("CovM_bg_tt");
    TMatrixD* CovM_bg_di = (TMatrixD*) Original.Get("CovM_bg_di");
    TMatrixD* CovM_bg_qcd = (TMatrixD*) Original.Get("CovM_bg_qcd");
    TMatrixD* CovM_pileup = (TMatrixD*) Original.Get("CovM_pileup");
    TMatrixD* CovM_pt = (TMatrixD*) Original.Get("CovM_pt");
    TMatrixD* CovM_lumi = (TMatrixD*) Original.Get("CovM_lumi");

    for (size_t binx = 0; binx < NumEst; binx++) {
        for (size_t biny = 0; biny < NumEst; biny++) {
            if (binx == biny)continue;
            if (RemoveCorrilation == 0)(*CovM_stat)(binx, biny) = 0;
            if (RemoveCorrilation == 1)(*CovM_mcstat)(binx, biny) = 0;
            if (RemoveCorrilation == 2)(*CovM_eff)(binx, biny) = 1e-20;
            if (RemoveCorrilation == 3)(*CovM_bg_tt)(binx, biny) = 0;
            if (RemoveCorrilation == 4)(*CovM_bg_di)(binx, biny) = 0;
            if (RemoveCorrilation == 5)(*CovM_bg_qcd)(binx, biny) = 0;
            if (RemoveCorrilation == 6)(*CovM_pileup)(binx, biny) = 0;
            if (RemoveCorrilation == 7)(*CovM_pt)(binx, biny) = 0;
        }
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
    //removed the lumi uncertainty from all calculations

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
    for (Int_t i = 0; i < NumEst; i++) { //filling blue
        myBlue->FillEst(i, &XEst[ind]);
        ind = ind + NumUnc + 1;
    }

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
    
    TMatrixD* Weights = new TMatrixD(NumEst,NumObs);
    myBlue->GetWeight(Weights);
    TMatrixD MatrixResults(NumObs, NumUnc + 1);

    TGraphAsymmErrors ElectronPlot(NumObs);
    TGraphAsymmErrors MuonPlot(NumObs);
    TGraphAsymmErrors BlueCombGraph(NumObs);
    myBlue->GetResult(&MatrixResults);



    double ElectronErrorSquared[NumObs];
    double MuonErrorSquared[NumObs];
    double CombinedErrorSquared[NumObs];

    for (size_t BinNumber = 0; BinNumber < NumObs; BinNumber++) { //creating erros from cov matrix
        ElectronErrorSquared[BinNumber] = AllCovs[0](BinNumber, BinNumber) + AllCovs[1](BinNumber, BinNumber) + AllCovs[2](BinNumber, BinNumber) + AllCovs[3](BinNumber, BinNumber) +
                AllCovs[4](BinNumber, BinNumber) + AllCovs[5](BinNumber, BinNumber) + AllCovs[6](BinNumber, BinNumber) + AllCovs[7](BinNumber, BinNumber);
        size_t muonBin = BinNumber + NumObs;
        MuonErrorSquared[BinNumber] = AllCovs[0](muonBin, muonBin) + AllCovs[1](muonBin, muonBin) + AllCovs[2](muonBin, muonBin) + AllCovs[3](muonBin, muonBin) +
                AllCovs[4](muonBin, muonBin) + AllCovs[5](muonBin, muonBin) + AllCovs[6](muonBin, muonBin) + AllCovs[7](muonBin, muonBin);
        CombinedErrorSquared[BinNumber] = (*CovarianceResults)(BinNumber, BinNumber);
    }

    for (size_t BinNumber = 0; BinNumber < NumObs; BinNumber++) { //creating individual electron and muon plot from combined plot
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
    if (RemoveCorrilation == -1) {
        if (DoNorm)filename = "Results/Comb_Norm_UsingBlue1D.root";
        else filename = "Results/Comb_Abs_UsingBlue1Dffcor.root"; //Weird thing getting rid of it
    } else if (RemoveCorrilation == 2) {
        if (DoNorm)filename = "Results/Comb_Norm_UsingBlue1DRemovedEffcor.root";
        else filename = "Results/Comb_Abs_UsingBlue1DRemovedEffcor.root";
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
    Weights->Write("Weights");


    TH2D CorrilationMatrix("CorrilationMatrix", "CorrilationMatrix", nbins, 0, nbins, nbins, 0, nbins);
    for (size_t BinX = 0; BinX < nbins; BinX++) {
        for (size_t BinY = 0; BinY < nbins; BinY++) {
            double Corrilation = (*CovarianceResults)(BinY, BinX) / sqrt((*CovarianceResults)(BinX, BinX)*(*CovarianceResults)(BinY, BinY));
            CorrilationMatrix.SetBinContent(BinX + 1, BinY + 1, Corrilation);
        }
    }

    TH1D ChiSquared("ChiSquared", "ChiSquared", 1, 0, 1);
    ChiSquared.SetBinContent(1, myBlue->GetChiq());
    ResultFile.Write();
    MatrixResults.Write("ResultsMatrix");

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

    TH1D Dumby("DUMBY", "", 1, .001, 3.277);
    if (DoNorm)Dumby.GetYaxis()->SetRangeUser(.95, 1.05);
    else Dumby.GetYaxis()->SetRangeUser(.95, 1.05);

    TCanvas* FinalPhiRatio = new TCanvas("Powheg", "PowhegPlot", 800, 900);

   
    FinalPhiRatio->SetRightMargin(0.01);
    FinalPhiRatio->SetLeftMargin(0.15);
    FinalPhiRatio->SetLogx();
    gStyle->SetOptStat("");
    FinalPhiRatio->cd();
    //FinalPhiRatio->SetLeftMargin(0);

    DataRatioFull->SetFillColor(kGray);
    DataRatioFull->SetTitle("");
    Dumby.GetYaxis()->SetTitle("Separate/Combined");
    Dumby.GetYaxis()->CenterTitle();
    Dumby.GetYaxis()->SetTitleOffset(1.5);
    Dumby.GetXaxis()->SetTitle("#phi*");
    Dumby.GetXaxis()->CenterTitle();
    Dumby.Draw();
    DataRatioFull->Draw("same E2");


    MuonRatioFull->SetMarkerColor(kRed);
    MuonRatioFull->SetLineColor(kRed);
    MuonRatioFull->SetMarkerStyle(21);
    MuonRatioFull->Draw("PE same");
    ElectronRatioFull->SetMarkerColor(kBlue);
    ElectronRatioFull->SetLineColor(kBlue);
    ElectronRatioFull->SetMarkerStyle(23);
    ElectronRatioFull->Draw("PEsame");

    TLegend* leg2 = new TLegend(0.15, 0.9, 0.99, 0.97);
    leg2->SetNColumns(3);
    leg2->SetFillStyle(0);
    //leg2->SetBorderSize(1);
    leg2->SetLineWidth(1);
    leg2->SetTextFont(22);
    double TextSize = 0.02;
    leg2->SetTextSize(TextSize);
    leg2->AddEntry(DataRatioFull, "BLUE Combined ", "F");
    leg2->AddEntry(MuonRatioFull, "2012 data Z #rightarrow #mu#mu", "P");
    leg2->AddEntry(ElectronRatioFull, "2012 Data Z #rightarrow ee", "P");
    leg2->Draw();
    FinalPhiRatio->RedrawAxis();
    TLatex mark;
    mark.SetTextSize(0.02);
    mark.SetNDC(kTRUE);
    mark.DrawLatex(.852, .974, "19.7 fb^{-1} (8 TeV)");
    mark.DrawLatex(0.15, .974, "CMS");
    string PlotName;
    if (DoNorm)PlotName = "Plots/OneDNorm";
    else PlotName = "Plots/OneDAbs";
    string Withtype = PlotName + ".pdf";
    FinalPhiRatio->Print(Withtype.c_str());
    Withtype = PlotName + ".png";
    FinalPhiRatio->Print(Withtype.c_str());
}

void Combiner() { //We don't like the eff causing a weird offset so instead we are now going to combine 2 results
    Int_t NumObs = nphistar;


    TFile NoEffCorFile("Results/Comb_Abs_UsingBlue1DRemovedEffcor.root", "read"); //grabbing

    TGraphAsymmErrors* h_Comb = (TGraphAsymmErrors*) NoEffCorFile.Get("h_Comb"); //combined results
    TGraphAsymmErrors* Electron = (TGraphAsymmErrors*) NoEffCorFile.Get("Electron"); //electron results
    if (!Electron) cout << "Missing electrons" << endl;
    TGraphAsymmErrors* Muons = (TGraphAsymmErrors*) NoEffCorFile.Get("Muons"); //muon results
    if (!Muons) cout << "Missing Muons" << endl;
    TH1D* ChiSquared = (TH1D*) NoEffCorFile.Get("ChiSquared"); //chisquarred
    
    
    TMatrixD* TotalCovarianceMatrix = (TMatrixD*) NoEffCorFile.Get("TotalCovarianceMatrix"); //total covariane matrix
    TMatrixD MatrixResults = *((TMatrixD*) NoEffCorFile.Get("ResultsMatrix"));
    TMatrixD Weights = *((TMatrixD*) NoEffCorFile.Get("Weights"));
    
    std::string OrignialFileLocation = "~/work/HomeWork/Phistar/CombineElectWithMu/EPlusMuFileV2/Comb_ForBlue_Abs_1D_Born.root"; //Grabbing efficency uncertainty shit
    TFile Original(OrignialFileLocation.c_str());
    TMatrixD CovM_eff = *((TMatrixD*) Original.Get("CovM_eff"));
    TMatrixD EffUncertainty(1,NumObs);
    
    for(size_t Obser=0; Obser<NumObs; Obser++){
        double uncertaintySquarred=0;
        for(size_t EstIndex = 0; EstIndex<NumEst; EstIndex++){
            if(Obser==0) cout << "Our uncertainty for EstIndex bin" << EstIndex << " is " <<sqrt(CovM_eff(EstIndex,EstIndex)); //prints all the values for the first event for sanity checking
            if(Obser==0) cout << "    Our Weight for each bin is " << Weights(EstIndex,Obser) ;
            uncertaintySquarred += Weights(EstIndex,Obser)*Weights(EstIndex,Obser)*CovM_eff(EstIndex,EstIndex);
            if(Obser==0) cout << "and our current uncertainty is " << sqrt(uncertaintySquarred) << endl;
        }
        if(Obser==0)cout << "AND OUR UNCERTAINTY IS " << sqrt(uncertaintySquarred);
        EffUncertainty(0,Obser)=sqrt(uncertaintySquarred);
    }
    
    cout << endl << endl << endl;
    
    
    if (!h_Comb) cout << "Missing h_Comb" << endl;
    NoEffCorFile.Close();
    
    if (!TotalCovarianceMatrix) cout << "Missing TotalCovarianceMatrix" << endl;

    for (size_t BinY = 0; BinY < h_Comb->GetN(); BinY++) { //goes over x axis of covarnaince bins
        double x1, y1;
        h_Comb->GetPoint(BinY, x1, y1);
        MatrixResults(BinY, 0) = y1;
        MatrixResults(BinY, 3) = EffUncertainty(0,BinY);
        cout << "For bin " << BinY << " our Eff uncertainty is " << EffUncertainty(0,BinY) << endl;
        
        for (size_t BinX = 0; BinX < h_Comb->GetN(); BinX++) { //goes over y axis
            double x2, y2;
            h_Comb->GetPoint(BinX, x2, y2);
            (*TotalCovarianceMatrix)(BinY, BinX) = (*TotalCovarianceMatrix)(BinY, BinX) + y1 * y2 * .026 * .026;
            (*TotalCovarianceMatrix)(BinY, BinX) = (*TotalCovarianceMatrix)(BinY, BinX) + EffUncertainty(0,BinX)*EffUncertainty(0,BinY);
        }
    }

    for (size_t binNumber = 0; binNumber < h_Comb->GetN(); binNumber++) {
        h_Comb->SetPointEYhigh(binNumber, sqrt((*TotalCovarianceMatrix)(binNumber, binNumber)));
        h_Comb->SetPointEYlow(binNumber, sqrt((*TotalCovarianceMatrix)(binNumber, binNumber)));
    }

    TFile* Finalthingy = new TFile("Results/Comb_Abs_UsingBlue1DNewEff.root", "recreate");
    Finalthingy->cd();
    Electron->Write("Electron");
    Muons->Write("Muons");
    h_Comb->Write("h_Comb");
    MatrixResults.Write("ResultsMatrix");
    TotalCovarianceMatrix->Write("TotalCovarianceMatrix");
    ChiSquared->Write();
    Finalthingy->Write();
}

int main(int argc, char * argv[]) {
    ElectPlusMuon1D(false); //not normalized
    ElectPlusMuon1D(false, 2); //normalized with no eff uncert
    Combiner(); //combines results
    ElectPlusMuon1D(true); //normalized results
    Plotter(true); //plots everything
    Plotter(false);
    return 1;
}
