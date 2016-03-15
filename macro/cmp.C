//This script is trying to compare sebastian's config and 2-D readout

void cmp(const char *file1, const char* file2, const char* file3="", 
	 const char* treename="ep")
{
  TChain *t1 = new TChain(treename);
  t1->Add(file1);

  TChain *t2 = new TChain(treename);
  t2->Add(file2);

  bool usefile3 = false;
  TChain *t3 = new TChain(treename);
  if(strlen(file3)>5)
    {
      t3->Add(file3);
      usefile3=true;
    }

  TCanvas *c11 = new TCanvas("c11","",800,600);
  
  c11->cd(0);
  t1->Draw("P0_p>>hp0_1","");
  t2->Draw("P0_p>>hp0_2","","same");
  if(usefile3) t3->Draw("P0_p>>hp0_3","","same");

  TH1F *hp0_1 = (TH1F*) gROOT->FindObject("hp0_1");
  TH1F *hp0_2 = (TH1F*) gROOT->FindObject("hp0_2");
  TH1F *hp0_3 = 0 ;
  if(usefile3) hp0_3 = (TH1F*) gROOT->FindObject("hp0_3");
  hp0_2->SetLineColor(2);
  if(usefile3) hp0_3->SetLineColor(4);
  hp0_1->SetTitle("file1(black),file2(red)");
  if(usefile3) hp0_1->SetTitle("file1(black),file2(red),file3(blue)");
 
}

void cmpvar(const char *file1, const char* file2, const char* treename,
	    const char *variable, const char *cut="")
{
  TChain *t1 = new TChain(treename);
  t1->Add(file1);

  TChain *t2 = new TChain(treename);
  t2->Add(file2);

  TCanvas *c11 = new TCanvas("c11","",800,600);
  
  c11->cd(0);
  t1->Draw(Form("%s>>h1",variable),cut);
  t1->Draw(Form("%s>>h2",variable),cut,"same");

  TH1F *h1 = (TH1F*) gROOT->FindObject("h1");
  TH1F *h2 = (TH1F*) gROOT->FindObject("h2");
  h2->SetLineColor(2);
  h1->SetTitle("file1(black),file2(red)");

}

void test(){
  TFile *_file1 = TFile::Open("nt_P-10T-10_100ns_0.6mm.root");
  TTree *t1=(TTree*)_file1->Get("ep");
  t1->SetName("tic100");
  TFile *_file2 = TFile::Open("nt_P-10T-10_200ns_0.6mm.root");
  TTree *t2=(TTree*)_file2->Get("ep");
  t2->SetName("tic200");
  
  t1->AddFriend(t2,"tic200");
  
  ((TTreePlayer*)(t1->GetPlayer()))->SetScanRedirect(true);
  ((TTreePlayer*)(t1->GetPlayer()))->SetScanFileName("output_0.6mm.txt");
  t1->Scan("StepS:StepID:StepTDC:HitNum:HitNum_m:tic200.StepTDC:tic200.HitNum_m","HitNum_m>5","", 10,0);
  TCanvas *c1=new TCanvas("c1","",800,600);
  t1->Draw("HitNum_m-tic200.HitNum_m");
  c1->SaveAs("diff_HitNum_m_0.6mm.png");

}
/*
//store scan result into output
((TTreePlayer*)(tree->GetPlayer()))->SetScanRedirect(true);
((TTreePlayer*)(tree->GetPlayer()))->SetScanFileName("output.txt");
tree->Scan(...)
t1->Scan("StepS:StepID:StepTDC:ep2.StepTDC:HitNum:HitNum_m:ep2.HitNum_m")

//there is another way
root> t1.SetScanField(0);  // no page break if Scan is redirected
root > t1.Scan(....); >scan.log
*/
