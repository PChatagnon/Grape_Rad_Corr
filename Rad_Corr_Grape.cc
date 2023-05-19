#include <TF1.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <iomanip>
#include <fstream>
#include <TRandom2.h>
#include <TLorentzVector.h>
#include "RadiativeCorrections.h"

#include <cstdlib>
#include <iostream>

using namespace std;


int main(int argc, char **argv)
{

    // ==================================
    // ==== Reading the input config file
    // ==================================

    ifstream inpconfig("GenOptions.dat");

    map<std::string, std::string> m_Settings;
    if (inpconfig.is_open())
    {
        while (!inpconfig.eof())
        {
            std::string Key;
            std::string Val;
            inpconfig >> Key;
            inpconfig >> Val;
            m_Settings[Key] = Val;
        }
    }
    else
    {
        cout << "Can not open the file GenOptions.dat" << endl;
        cout << "So can not initialize settings " << endl;
        cout << "Exiting" << endl;
        exit(1);
    }

    bool Rad_corr;
    double rad_cut_off_min;
    double rad_cut_off_max;
    bool isLund;

    for (map<std::string, std::string>::iterator it = m_Settings.begin(); it != m_Settings.end(); it++)
    {

        std::string key = (*it).first;
        std::string val = (*it).second;

        if (key.compare("RAD_CORR") == 0)
        {
            Rad_corr = atof(val.c_str());
        }
        else if (key.compare("rad_cut_off_min") == 0)
        {
            rad_cut_off_min = atof(val.c_str());
        }
        else if (key.compare("rad_cut_off_max") == 0)
        {
            rad_cut_off_max = atof(val.c_str());
        }
        else if (key.compare("LUND") == 0)
        {
            isLund = atoi(val.c_str());
        }
    }

    cout << "Rad_corr = " << Rad_corr << endl;
    cout << "rad_cut_off_min = " << rad_cut_off_min << endl;
    cout << "rad_cut_off_max = " << rad_cut_off_max << endl;

    const double Mp = 0.9383;
    const double Me = 0.00051;

    const int nHeader = 10;
    const int npartColumns = 14;

    int nPart;
    double A_Targ;
    int Z_Targ;
    double pol_targ;
    double pol_beam;
    int beamType; // 11 Electron, 22 Photon
    double Eb;
    int InterNuclID; // Interacted Nucleon ID
    int ProcessID;
    double EvWeight;

    const int nMax = 50;

    int index[nMax];
    double t_live[nMax]; // lifetime [ns]
    int type[nMax];      // 1=active
    int pid[nMax];
    int parentInd[nMax];
    int daughtInd[nMax];
    double px[nMax];
    double py[nMax];
    double pz[nMax];
    double E[nMax];
    double m[nMax];
    double vx[nMax];
    double vy[nMax];
    double vz[nMax];

    TLorentzVector L_em, L_ep, L_prot, L_prot_target, L_rad_1, L_rad_2;
    TLorentzVector Lemep;
    L_prot_target.SetPxPyPzE(0., 0., 0.0, Mp);

    double px_prot, py_prot, pz_prot, E_prot;
    double px_ep, py_ep, pz_ep, E_ep;
    double px_em, py_em, pz_em, E_em;
    double px_rad_em, py_rad_em, pz_rad_em, E_rad_em;
    double px_rad_ep, py_rad_ep, pz_rad_ep, E_rad_ep;
    double px_rad, py_rad, pz_rad, E_rad;
    double Theta_rad, Phi_rad, Angle_g_lep;
    double E_rad_cm, Theta_rad_cm, Phi_rad_cm, Angle_g_lep_cm;
    double Inv_Mass;
    double vertex_z;
    double Eg, t, Q2;

    TFile *file_out;
    if (!isLund)
    {
        file_out = new TFile("tcs_gen.root", "Recreate");
    }

    TTree *tr1 = new TTree("tr1", "TCS MC events");
    tr1->Branch("Eg", &Eg, "Eg/D");
    tr1->Branch("Q2", &Q2, "Q2/D");
    tr1->Branch("t", &t, "t/D");

    tr1->Branch("px_prot", &px_prot, "px_prot/D");
    tr1->Branch("py_prot", &py_prot, "py_prot/D");
    tr1->Branch("pz_prot", &pz_prot, "pz_prot/D");
    tr1->Branch("E_prot", &E_prot, "E_prot/D");
    tr1->Branch("px_ep", &px_ep, "px_ep/D");
    tr1->Branch("py_ep", &py_ep, "py_ep/D");
    tr1->Branch("pz_ep", &pz_ep, "pz_ep/D");
    tr1->Branch("E_ep", &E_ep, "E_ep/D");
    tr1->Branch("px_em", &px_em, "px_em/D");
    tr1->Branch("py_em", &py_em, "py_em/D");
    tr1->Branch("pz_em", &pz_em, "pz_em/D");
    tr1->Branch("E_em", &E_em, "E_em/D");
    tr1->Branch("px_rad", &px_rad, "px_rad/D");
    tr1->Branch("py_rad", &py_rad, "py_rad/D");
    tr1->Branch("pz_rad", &pz_rad, "pz_rad/D");
    tr1->Branch("px_rad_em", &px_rad_em, "px_rad_em/D");
    tr1->Branch("py_rad_em", &py_rad_em, "py_rad_em/D");
    tr1->Branch("pz_rad_em", &pz_rad_em, "pz_rad_em/D");
    tr1->Branch("px_rad_ep", &px_rad_ep, "px_rad_ep/D");
    tr1->Branch("py_rad_ep", &py_rad_ep, "py_rad_ep/D");
    tr1->Branch("pz_rad_ep", &pz_rad_ep, "pz_rad_ep/D");
    tr1->Branch("E_rad", &E_rad, "E_rad/D");
    tr1->Branch("E_rad_cm", &E_rad_cm, "E_rad_cm/D");
    tr1->Branch("Inv_Mass", &Inv_Mass, "Inv_Mass/D");

    tr1->Branch("Theta_rad", &Theta_rad, "Theta_rad/D");
    tr1->Branch("Phi_rad", &Phi_rad, "Phi_rad/D");
    tr1->Branch("Angle_g_lep", &Angle_g_lep, "Angle_g_lep/D");

    tr1->Branch("Theta_rad_cm", &Theta_rad_cm, "Theta_rad_cm/D");
    tr1->Branch("Phi_rad_cm", &Phi_rad_cm, "Phi_rad_cm/D");
    tr1->Branch("Angle_g_lep_cm", &Angle_g_lep_cm, "Angle_g_lep_cm/D");

    fstream Lund_in;
    fstream Lund_out;

    TRandom2 rand;
    rand.SetSeed(0);

    /////////////Set Rad Corr Parameters//////////////
    RadiativeCorrections Rad_corr_1(rad_cut_off_min, rad_cut_off_max);
    /////////////////////////////////////////////////

    for (int i = 1; i < argc; i++)
    {
        string input_lund_file = argv[i];
        cout << input_lund_file << endl;
        Lund_in.open(input_lund_file, ios::in | ios::app);
        Lund_out.open(input_lund_file.erase(input_lund_file.size() - 4) + "_radiated.txt", ofstream::out);

        if (Lund_in.is_open())
        {

            std::string line;

            while (std::getline(Lund_in, line))
            {

                std::vector<std::string> result;
                std::istringstream iss(line);

                for (std::string s; iss >> s;)
                {
                    result.push_back(s);
                }

                if (result.size() == nHeader)
                {
                    nPart = atoi(result.at(0).c_str());
                    A_Targ = atoi(result.at(1).c_str());
                    Z_Targ = atoi(result.at(2).c_str());
                    pol_targ = atof(result.at(3).c_str());
                    pol_beam = atof(result.at(4).c_str());
                    beamType = atoi(result.at(5).c_str());
                    Eb = atof(result.at(6).c_str());
                    InterNuclID = atoi(result.at(7).c_str());
                    ProcessID = atoi(result.at(8).c_str());
                    EvWeight = atof(result.at(9).c_str());

                    for (int ipart = 0; ipart < nPart; ipart++)
                    {
                        result.clear();
                        result.shrink_to_fit();

                        std::getline(Lund_in, line);

                        std::istringstream iss(line);

                        for (std::string s; iss >> s;)
                        {
                            result.push_back(s);
                        }

                        index[ipart] = atoi(result.at(0).c_str());
                        t_live[ipart] = atof(result.at(1).c_str());
                        type[ipart] = atoi(result.at(2).c_str());
                        pid[ipart] = atoi(result.at(3).c_str());
                        parentInd[ipart] = atoi(result.at(4).c_str());
                        daughtInd[ipart] = atoi(result.at(5).c_str());
                        px[ipart] = atof(result.at(6).c_str());
                        py[ipart] = atof(result.at(7).c_str());
                        pz[ipart] = atof(result.at(8).c_str());
                        E[ipart] = atof(result.at(9).c_str());
                        m[ipart] = atof(result.at(10).c_str());
                        vx[ipart] = atof(result.at(11).c_str());
                        vy[ipart] = atof(result.at(12).c_str());
                        vz[ipart] = atof(result.at(13).c_str());
                    }

                    ////////////Assign 4 momenta//////////////////
                    L_em.SetPxPyPzE(px[3], py[3], pz[3], E[3]);
                    L_ep.SetPxPyPzE(px[2], py[2], pz[2], E[2]);
                    L_prot.SetPxPyPzE(px[0], py[0], pz[0], E[0]);
                    Lemep = L_em + L_ep;

                    Eg = (L_em+L_ep+L_prot-L_prot_target).E();
                    t = (L_prot_target-L_prot).M2();
                    Q2 = (L_em+L_ep).M2();
                    /*cout<<"check CoM before rad"<<endl;
                    cout<<(L_em+L_ep).Px()<<"  "<<(L_em+L_ep).Py()<<"  "<<(L_em+L_ep).Pz()<<endl;
                    cout<<(L_em).Px()<<"  "<<(L_em).Py()<<"  "<<(L_em).Pz()<<endl;
                    cout<<(L_ep).Px()<<"  "<<(L_ep).Py()<<"  "<<(L_ep).Pz()<<endl;*/

                    ////////////Boost in CoM//////////////////////
                    L_em.Boost(-Lemep.BoostVector());
                    L_ep.Boost(-Lemep.BoostVector());
                    L_prot.Boost(-Lemep.BoostVector());
                    //////////////////////////////////////////////

                    /////////////////Radiative correction must be performed in the CoM frame///////////////////
                    bool in_rad_tail = (rand.Uniform(0, 1) > Rad_corr_1.Compute_cs_correction_factor((L_em + L_ep).M())); // randomly choose if the photon is above cut_off_min

                    if (Rad_corr && in_rad_tail)
                        Rad_corr_1.Soft_Photon_Emission(L_em, L_ep, L_rad_1, L_rad_2);
                    else
                    {
                        L_rad_1.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
                        L_rad_2.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
                    }
                    //////////////////////////////////////////////////////////

                    E_rad_cm = (L_rad_1 + L_rad_2).E();
                    Theta_rad_cm = L_rad_1.Theta();
                    Phi_rad_cm = L_rad_1.Phi();
                    Angle_g_lep_cm = L_rad_1.Angle(L_em.Vect());

                    ////////////Boost back in lab//////////////////////
                    L_em.Boost(Lemep.BoostVector());
                    L_ep.Boost(Lemep.BoostVector());
                    L_prot.Boost(Lemep.BoostVector());
                    L_rad_1.Boost(Lemep.BoostVector());
                    L_rad_2.Boost(Lemep.BoostVector());
                    //////////////////////////////////////////////
                    /*if(in_rad_tail){cout<<"check CoM after rad"<<endl;
                    cout<<(L_em+L_ep).Px()<<"  "<<(L_em+L_ep).Py()<<"  "<<(L_em+L_ep).Pz()<<endl;
                    cout<<(L_em).Px()<<"  "<<(L_em).Py()<<"  "<<(L_em).Pz()<<endl;
                    cout<<(L_ep).Px()<<"  "<<(L_ep).Py()<<"  "<<(L_ep).Pz()<<endl;}*/

                    px_prot = L_prot.Px();
                    py_prot = L_prot.Py();
                    pz_prot = L_prot.Pz();
                    E_prot = L_prot.E();
                    px_ep = L_ep.Px();
                    py_ep = L_ep.Py();
                    pz_ep = L_ep.Pz();
                    E_ep = L_ep.E();
                    px_em = L_em.Px();
                    py_em = L_em.Py();
                    pz_em = L_em.Pz();
                    E_em = L_em.E();
                    px_rad = (L_rad_1 + L_rad_2).Px();
                    py_rad = (L_rad_1 + L_rad_2).Py();
                    pz_rad = (L_rad_1 + L_rad_2).Pz();
                    E_rad = (L_rad_1 + L_rad_2).E();
                    px_rad_em = L_rad_1.Px();
                    py_rad_em = L_rad_1.Py();
                    pz_rad_em = L_rad_1.Pz();
                    E_rad_em = L_rad_1.E();
                    px_rad_ep = L_rad_2.Px();
                    py_rad_ep = L_rad_2.Py();
                    pz_rad_ep = L_rad_2.Pz();
                    E_rad_ep = L_rad_2.E();
                    vertex_z = vz[0];

                    Inv_Mass = (L_em + L_ep).M();

                    Theta_rad = L_rad_1.Theta();
                    Phi_rad = L_rad_1.Phi();
                    Angle_g_lep = L_rad_1.Angle(L_em.Vect());

                    if (!isLund)
                    {
                        tr1->Fill();
                    }

                    // Writing Header
                    Lund_out << 6 << setw(5) << 1 << setw(5) << 1 << setw(5) << 0 << " " << setw(5) << "  " << 0 << " " << setw(15) << 0 << setw(15) << 0 << setw(15) << 0 << setw(5) << 0 << setw(5) << " " << 0 << "\n"; // Writing Proton
                    Lund_out << 1 << setw(5) << 1 << setw(5) << 1 << setw(7) << 2212 << setw(5) << 0 << setw(5) << 0 << setw(15) << px_prot << setw(15) << py_prot << setw(15) << pz_prot;
                    Lund_out << setw(15) << E_prot << setw(15) << Mp << setw(15) << 0. << setw(15) << 0. << setw(15) << vertex_z << "\n";
                    // Writing Beam
                    Lund_out << 2 << setw(5) << -1 << setw(5) << 1 << setw(7) << 11 << setw(5) << 0 << setw(5) << 0 << setw(15) << px[1] << setw(15) << py[1] << setw(15) << pz[1];
                    Lund_out << setw(15) << E[1] << setw(15) << Me << setw(15) << 0. << setw(15) << 0. << setw(15) << vertex_z << "\n";
                    // Writing Positron
                    Lund_out << 3 << setw(5) << 1 << setw(5) << 1 << setw(7) << -11 << setw(5) << 0 << setw(5) << 0 << setw(15) << px_ep << setw(15) << py_ep << setw(15) << pz_ep;
                    Lund_out << setw(15) << E_ep << setw(15) << Me << setw(15) << 0. << setw(15) << 0. << setw(15) << vertex_z << "\n";
                    // Writing Electron
                    Lund_out << 4 << setw(5) << -1 << setw(5) << 1 << setw(7) << 11 << setw(5) << 0 << setw(5) << 0 << setw(15) << px_em << setw(15) << py_em << setw(15) << pz_em;
                    Lund_out << setw(15) << E_em << setw(15) << Me << setw(15) << 0. << setw(15) << 0. << setw(15) << vertex_z << "\n";
                    // Writing Photons
                    Lund_out << 5 << setw(5) << 0 << setw(5) << 1 << setw(7) << 22 << setw(5) << 0 << setw(5) << 0 << setw(15) << px_rad_em << setw(15) << py_rad_em << setw(15) << pz_rad_em;
                    Lund_out << setw(15) << E_rad_em << setw(15) << 0.0 << setw(15) << 0. << setw(15) << 0. << setw(15) << vertex_z << "\n";
                    Lund_out << 6 << setw(5) << 0 << setw(5) << 1 << setw(7) << 22 << setw(5) << 0 << setw(5) << 0 << setw(15) << px_rad_ep << setw(15) << py_rad_ep << setw(15) << pz_rad_ep;
                    Lund_out << setw(15) << E_rad_ep << setw(15) << 0.0 << setw(15) << 0. << setw(15) << 0. << setw(15) << vertex_z << "\n";
                }
            }
        }

        Lund_in.close();
        Lund_out.close();
    }

    tr1->Write();
    file_out->Close();

    return 0;
}
