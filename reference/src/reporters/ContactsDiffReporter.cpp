#include "reporters/ContactsDiffReporter.hpp"
#include "utils/Units.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;
using namespace cg::reference;

using CType = QuasiAdiabatic::ContactType;

ContactsDiffReporter::ContactsDiffReporter(Results const &results,
                                           string filePath):
    results(results) {

    ifstream input(filePath, ifstream::in);
    assert (input);

    string token, ignore;
    unsigned step_nr;

    while (input >> token) {
        if (token != "K") continue;

        int step;
        input >> step;
        
        QuasiAdiabatic::Contact contact;
        contact.adiabCoeff = 0;

        input >> contact.i >> contact.j;
        contact.i--; contact.j--;

        input >> ignore;

        int type;
        input >> type;

        if (type > 0) { // skip native contacts
            if (type == 4) contact.type = CType::BB;
            else if (type == 5) contact.type = CType::SS;
            else if (type == 6) contact.type = CType::BS;
            else if (type == 7) contact.type = CType::SB;
            else assert(0);
            
            refContacts[step].push_back(contact);
        }
    }
}

string cTypeToName(CType type) {
    switch (type) {
    case CType::NONE:
        return "NONE";
    case CType::BB:
        return "BB";
    case CType::BS:
        return "BS";
    case CType::SB:
        return "SB";
    case CType::SS:
        return "SS";
    }
    return "";
}

void ContactsDiffReporter::report(int step) {
    auto it = refContacts.find(step);
    if (!results.qaContacts) {
        cerr << "Information about qa contacts wasn't dumped in step " << step << endl;
        return;
    }
    if (it != refContacts.end()) {
        assert (results.qaContacts);

        map<pair<int, int>, CType> contacts;
        set<pair<int, int>> refContacts;

        for (QuasiAdiabatic::Contact const &cnt : *results.qaContacts) {
            contacts[{cnt.i, cnt.j}] = cnt.type;
        }

        cout << "Step #" << setw(8) << step << endl;

        for (QuasiAdiabatic::Contact const &refCnt : it -> second) {
            refContacts.insert({refCnt.i, refCnt.j});
            CType type = CType::NONE;

            auto jt = contacts.find({refCnt.i, refCnt.j});
            if (jt != contacts.end()) type = jt -> second;

            if (type != refCnt.type) {
                cout << "Bad contact type. " << "i = " << refCnt.i << ", j = " << refCnt.j
                     << ", type = " << cTypeToName(type)
                     << ", refType = " << cTypeToName(refCnt.type) << endl;
            }
        }

        for (QuasiAdiabatic::Contact const &cnt : *results.qaContacts) {
            auto jt = refContacts.find({cnt.i, cnt.j});
            if (jt == refContacts.end()) {
                cout << "Contact shouldn't exist: " << "i = " << cnt.i << ", j = " << cnt.j
                     << ", type = " << cTypeToName(cnt.type) << endl;
            }
        }
    }
    else {
        cerr << "No information about contacts in step " << step << endl;
    }
}
