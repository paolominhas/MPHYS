
# this is just old spare code


#include 

void saveGraph() {
    TFile *f = TFile::Open("../hibeam_g4/outputs/MuonLargeData.root");
    TTree *t = (TTree*)f->Get("hibeam;1");

    Double_t branch_variable; // use correct data type double_t
    t->SetBranchAddress("branch_name", &branch_variable);

    ofstream csvFile("output.csv");
    csvFile << "hibeam_branch" << endl; // optional header

    Long64_t nentries = t->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        t->GetEntry(i);
        csvFile << branch_variable << endl;
    }
    csvFile.close();
}