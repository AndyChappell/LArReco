/**
 *  @file   ProcessValidation.C
 *
 *  @brief  Some helper routines to process the in memory histograms produced
 * 	    by LarReco/validation/Validation.C. That is, currently you need
 *	    to run this in the same session where you compiled Validation.C.
 *
 *  $Log: $
 */
#include "TROOT.h"
#include "TStyle.h"

#include <iostream>
#include "TCanvas.h"
#include "TChain.h"
#include "TCollection.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TPRegexp.h"
#include "../LArReco/validation/Validation.h"

typedef std::map<std::string, TH1F*> InteractionToHistogramMap;

void ProcessValidation(const char* filename)
{
    gROOT->ProcessLine(".L ../LArReco/validation/Validation.C++");
    Parameters parameters;
    parameters.m_applyUbooneFiducialCut = false;
    parameters.m_histogramOutput = true;
    
    Validation(filename, parameters);
    
    // Display the in-memory histograms
    gDirectory->ls("-m");
}

void PlotPurity(const std::string& interaction)
{
    gStyle->SetOptStat(0);
    const std::string str(interaction + "_MUON_Purity");
    TH1F* purity = (TH1F*)gDirectory->Get(str.c_str());
    purity->SetLineColor(kBlack);

    TCanvas c("c", "c", 1024, 768);
    //purity->GetXaxis()->SetTitle("Purity");
    purity->Draw();
    const std::string filename(interaction + "_MUON_Purity.pdf");
    
    TLegend legend(0.15, 0.75, 0.35, 0.85);
    legend.AddEntry(purity, "Purity", "l");
    legend.SetFillStyle(0);
    legend.Draw();
    
    c.SetLogy();
    c.Print(filename.c_str());
}

void PlotCompleteness(const std::string& interaction)
{
    gStyle->SetOptStat(0);
    const std::string str(interaction + "_MUON_Completeness");
    TH1F* completeness = (TH1F*)gDirectory->Get(str.c_str());
    completeness->SetLineColor(kBlack);

    TCanvas c("c", "c", 1024, 768);
    //purity->GetXaxis()->SetTitle("Completeness");
    completeness->Draw();
    const std::string filename(interaction + "_MUON_Completness.pdf");
    
    TLegend legend(0.15, 0.75, 0.35, 0.85);
    legend.AddEntry(completeness, "Completeness", "l");
    legend.SetFillStyle(0);
    legend.Draw();
    
    c.SetLogy();
    c.Print(filename.c_str());
}

void PlotMerge(const InteractionToHistogramMap& map, const std::string& type);
void PlotMergeCumulative(const InteractionToHistogramMap& map, const std::string& type);
void MergeKeys()
{
    gStyle->SetOptStat(0);
    TIter next(gDirectory->GetList());
    std::string reInteraction("\\b(CCQEL\\_MU|CCQEL\\_E|CCRES\\_MU|CCRES\\_E|CCDIS\\_MU|CCDIS\\_E|NCQEL|NCRES|NCDIS|CCCOH|NCCOH|OTHER\\_INTERACTION)");
    std::string reProtons("(\\_P\\_P\\_P\\_P\\_P\\_|\\_P\\_P\\_P\\_P\\_|\\_P\\_P\\_P\\_|\\_P\\_P\\_|\\_P\\_)?\\_?");
    std::string rePrimaries("(PIPLUS|PIMINUS|PIZERO|PHOTON|NULL)?\\_?");
    std::string reParticle("(MUON|ELECTRON|PIPLUS|PIMINUS|PHOTON1|PHOTON2|PHOTON3|PHOTON4|PHOTON5|PROTON1|PROTON2|PROTON3|PROTON4|PROTON5|NULL])?\\_?");
    std::string reMetrics("(Purity|Completeness)\\b");
    std::string reStr(reInteraction + reProtons + rePrimaries + reParticle + reMetrics);
    TPRegexp re(reStr.c_str());

    InteractionToHistogramMap mapI2P, mapI2C;
    while (TObject* obj = next())
    {
        TString name(obj->GetName());
        TObjArray* arr = re.MatchS(name);
        int N = arr->GetEntries();
        if (N > 0)
        {
            std::string interaction = std::string(((TObjString*)arr->At(1))->GetString().Data());
            std::string protons = std::string(((TObjString*)arr->At(2))->GetString().Data());
            std::string primary = std::string(((TObjString*)arr->At(3))->GetString().Data());
            std::string particle = std::string(((TObjString*)arr->At(4))->GetString().Data());
            std::string metric = std::string(((TObjString*)arr->At(5))->GetString().Data());
            if (particle.empty())
            {
                particle = primary;
                primary.clear();
            }
            if (primary == "PIPLUS" || primary == "PIMINUS")
                primary = "PIC";
            std::string interactionKey = interaction +
                (primary.empty() ? "" : "_" + primary) +
                //(protons.empty() ? "" : (protons == "_P_" ? "_P" : "_NP")) +
                "_" + particle;
            if (metric == "Purity")
            {
                auto iterI2P = mapI2P.find(interactionKey);
                if (iterI2P != mapI2P.end())
                {
                    TH1F* current = iterI2P->second;
                    current->Scale(current->GetEntries(), "nosw2");
                    TH1F* addend = static_cast<TH1F*>(gDirectory->Get(name));
                    addend->Scale(addend->GetEntries(), "nosw2");
                    current->Add(addend);
                    current->Scale(1.f / current->GetEntries(), "nosw2");
                }
                else
                {
                    TH1F* hist = static_cast<TH1F*>(gDirectory->Get(name));
                    mapI2P.emplace(interactionKey, hist);
                }
            }
            else
            {
                auto iterI2C = mapI2C.find(interactionKey);
                if (iterI2C != mapI2C.end())
                {
                    TH1F* current = iterI2C->second;
                    current->Scale(current->GetEntries(), "nosw2");
                    TH1F* addend = static_cast<TH1F*>(gDirectory->Get(name));
                    addend->Scale(addend->GetEntries(), "nosw2");
                    current->Add(addend);
                    current->Scale(1.f / current->GetEntries(), "nosw2");
                }
                else
                {
                    TH1F* hist = static_cast<TH1F*>(gDirectory->Get(name));
                    mapI2C.emplace(interactionKey, hist);
                }
            }    
        }
        delete arr;
    }

    //PlotMerge(mapI2C, "Completeness");
    //PlotMerge(mapI2P, "Purity");
    
    PlotMergeCumulative(mapI2C, "Completeness");
    PlotMergeCumulative(mapI2P, "Purity");
}

void PlotMerge(const InteractionToHistogramMap& map, const std::string& type)
{
    std::map<std::string, int> particleColor;
    particleColor.emplace("MUON", kMagenta + 2);
    particleColor.emplace("ELECTRON", kRed + 1);
    particleColor.emplace("PIPLUS", kGreen + 2);
    particleColor.emplace("PIMINUS", kCyan + 2);
    particleColor.emplace("PHOTON1", kOrange);
    particleColor.emplace("PHOTON2", kOrange + 1);
    particleColor.emplace("PHOTON3", kOrange + 2);
    particleColor.emplace("PHOTON4", kOrange + 3);
    particleColor.emplace("PHOTON5", kOrange + 4);
    particleColor.emplace("PROTON1", kBlue);
    particleColor.emplace("PROTON2", kBlue + 2);
    particleColor.emplace("PROTON3", kBlue - 4);
    particleColor.emplace("PROTON4", kBlue - 7);
    particleColor.emplace("PROTON5", kBlue + 4);

    std::map<std::string, TCanvas*> canvases;
    std::map<std::string, TLegend*> legends;
    for (const auto [key, value] : map)
    {
        const int splitPos = key.find_last_of("_");
        std::string interaction = key.substr(0, splitPos);
        std::string particle = key.substr(splitPos + 1);
        auto iter = canvases.find(interaction);
        if (iter != canvases.end())
        {   // Update existing plot
            iter->second->cd();
            TLegend* legend = legends.find(interaction)->second;
            legend->AddEntry(value, particle.c_str(), "l");
            value->SetLineColor(particleColor.find(particle)->second);
            value->SetLineWidth(2);
            value->Draw("hist same");
            legend->Draw();
        }
        else
        {   // Create new plot
            TCanvas* canvas = new TCanvas(interaction.c_str(), interaction.c_str(), 1600, 900);
            canvases.emplace(interaction, canvas);
            TLegend* legend(new TLegend(0.15, 0.65, 0.4, 0.88));
            legend->SetFillStyle(0);
            legend->SetBorderSize(0);
            legend->SetTextSize(0.03);
            legend->AddEntry(value, particle.c_str(), "l");
            legends.emplace(interaction, legend);
            value->SetLineColor(particleColor.find(particle)->second);
            value->SetLineWidth(2);
            // Don't set lower limit of zero, or log scale won't work
            value->GetYaxis()->SetRangeUser(std::numeric_limits<float>::epsilon(), 1);
            canvas->cd();
            gPad->SetLogy();
            value->Draw("hist");
            legend->Draw();
        }
    }
    
    for (const auto [key, canvas] : canvases)
    {
        std::string filename = "pdfs/" + type + "_" + key + ".pdf";
        canvas->Print(filename.c_str());
        delete canvas;
    }
    canvases.clear();
    for (const auto [key, legend] : legends)
    {
        delete legend;
    }
    legends.clear();
}

void PlotMergeCumulative(const InteractionToHistogramMap& map, const std::string& type)
{
    std::map<std::string, int> particleColor;
    particleColor.emplace("MUON", kMagenta + 2);
    particleColor.emplace("ELECTRON", kRed + 1);
    particleColor.emplace("PIPLUS", kGreen + 2);
    particleColor.emplace("PIMINUS", kCyan + 2);
    particleColor.emplace("PHOTON1", kOrange);
    particleColor.emplace("PHOTON2", kOrange + 1);
    particleColor.emplace("PHOTON3", kOrange + 2);
    particleColor.emplace("PHOTON4", kOrange + 3);
    particleColor.emplace("PHOTON5", kOrange + 4);
    particleColor.emplace("PROTON1", kBlue);
    particleColor.emplace("PROTON2", kBlue + 2);
    particleColor.emplace("PROTON3", kBlue - 4);
    particleColor.emplace("PROTON4", kBlue - 7);
    particleColor.emplace("PROTON5", kBlue + 4);

    std::map<std::string, TCanvas*> canvases;
    std::map<std::string, TLegend*> legends;
    for (const auto [key, value] : map)
    {
        const int splitPos = key.find_last_of("_");
        std::string interaction = key.substr(0, splitPos);
        std::string particle = key.substr(splitPos + 1);
        auto iter = canvases.find(interaction);
        if (iter != canvases.end())
        {   // Update existing plot
            iter->second->cd();
            TLegend* legend = legends.find(interaction)->second;
            legend->AddEntry(value, particle.c_str(), "l");
            value->SetLineColor(particleColor.find(particle)->second);
            value->SetLineWidth(2);
            for (int i = 2; i <= value->GetNbinsX(); ++i)
                value->SetBinContent(i, value->GetBinContent(i) + value->GetBinContent(i - 1));
            value->Draw("hist same");
            legend->Draw();
        }
        else
        {   // Create new plot
            TCanvas* canvas = new TCanvas(interaction.c_str(), interaction.c_str(), 1600, 900);
            canvases.emplace(interaction, canvas);
            TLegend* legend(new TLegend(0.15, 0.65, 0.4, 0.88));
            legend->SetFillStyle(0);
            legend->SetBorderSize(0);
            legend->SetTextSize(0.03);
            legend->AddEntry(value, particle.c_str(), "l");
            legends.emplace(interaction, legend);
            value->SetLineColor(particleColor.find(particle)->second);
            value->SetLineWidth(2);
            for (int i = 2; i <= value->GetNbinsX(); ++i)
                value->SetBinContent(i, value->GetBinContent(i) + value->GetBinContent(i - 1));
            // Don't set lower limit of zero, or log scale won't work
            value->GetYaxis()->SetRangeUser(std::numeric_limits<float>::epsilon(), 1);
            canvas->cd();
            gPad->SetLogy();
            value->Draw("hist");
            legend->Draw();
        }
    }
    
    for (const auto [key, canvas] : canvases)
    {
        std::string filename = "pdfs/c_" + type + "_" + key + ".pdf";
        canvas->Print(filename.c_str());
        delete canvas;
    }
    canvases.clear();
    for (const auto [key, legend] : legends)
    {
        delete legend;
    }
    legends.clear();
}
