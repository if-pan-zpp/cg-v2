#include "forces/nonlocal/QuasiAdiabatic.hpp"
#include "data/Results.hpp"
using namespace cg::reference;

QuasiAdiabatic::QuasiAdiabatic(PseudoAtoms const &pseudoAtoms,
                               Topology const &top,
                               Neighborhood const &verletList,
                               toolkit::Parameters const &parameters,
                               SharedData &sharedData):
    pseudoAtoms(pseudoAtoms),
    top(top),
    verletList(verletList),
    sharedData(sharedData) {


    vector<bool> isSet(NUM_AMINO_ACIDS);
    for (auto const &[aminoAcid, specificity] : parameters.specificities) {
        AACode aaCode = aaCodeFromName(string(aminoAcid));
        CoordNumber crdNum;

        crdNum.backbone = 2;
        if (aaCode == PRO) crdNum.backbone = 1;
        
        crdNum.sidechain = specificity.coordinationNumber;
        crdNum.ssHydrophobic = specificity.hydrophobicCoordinationNumber;
        crdNum.ssPolar = specificity.polarCoordinationNumber;

        coordBounds[aaCode] = crdNum;
        // TODO: make this conversion better
        aaTypes[aaCode] = (QuasiAdiabatic::AAType) (int) specificity.polarity;
            
        isSet[aaCode] = true;
    }
    
    for (size_t code = 0; code < NUM_AMINO_ACIDS; ++code) {
        if (!isSet[code]) {
            assert(false);
            // TODO: add some message
        }
    }

    for (auto const &[aminoAcidPair, r_min] : parameters.ssMinimalDistances) {
        AACode ai = aaCodeFromName(string(aminoAcidPair.first));
        AACode aj = aaCodeFromName(string(aminoAcidPair.second));
        ss_lj_r_min[ai][aj] = r_min;
    }

    totalCutoff = *max_element(lj_r_min, lj_r_min + NUM_CONTACT_TYPES);
    for (size_t ai = 0; ai < NUM_AMINO_ACIDS; ++ai) {
        for (size_t aj = 0; aj < NUM_AMINO_ACIDS; ++aj) {
            totalCutoff = max(totalCutoff, ss_lj_r_min[ai][aj]);
        }
    }
    totalCutoff *= cutoffCoeff;
    totalCutoffSquared = totalCutoff * totalCutoff;
}

uint8_t QuasiAdiabatic::getSSBound(CoordNumber const &crdNum, AACode aaCode) const {
    switch (aaTypes[aaCode]) {

    case GLY_OR_PRO:
        return crdNum.sidechain;

    case HYDROPHOBIC:
        return crdNum.ssHydrophobic;

    case POLAR:
    case CHARGED_NEG:
    case CHARGED_POS:
        return crdNum.ssPolar;

    default:
        __builtin_unreachable();
        return 0;
    };
}

Real3 QuasiAdiabatic::hVector(int i) const {
    assert (i > 0 && i + 1 < pseudoAtoms.n);
    Integers const &chainId = pseudoAtoms.chainId;
    assert (chainId[i - 1] == chainId[i + 1]);

    Reals3 const &pos = pseudoAtoms.pos;
    Real3 va = pos.col(i + 1) - pos.col(i);
    Real3 vb = pos.col(i) - pos.col(i - 1);
    return va.cross(vb).normalized();
}

Real3 QuasiAdiabatic::nVector(int i) const {
    assert (i > 0 && i + 1 < pseudoAtoms.n);
    Integers const &chainId = pseudoAtoms.chainId;
    assert (chainId[i - 1] == chainId[i + 1]);

    Reals3 const &pos = pseudoAtoms.pos;
    Real3 v0 = pos.col(i + 1) - pos.col(i);
    Real3 v1 = pos.col(i) - pos.col(i - 1);
    return (v0 - v1).normalized();
}

bool QuasiAdiabatic::isDisabled(ContactType type, int i, int j) const {
    CoordNumber const &crd_i = sharedData.coordNumbers[i];
    CoordNumber const &crd_j = sharedData.coordNumbers[j];
    
    AACode const aaCode_i = pseudoAtoms.aminoAcidCode[i];
    AACode const aaCode_j = pseudoAtoms.aminoAcidCode[j];

    CoordNumber const &crd_bound_i = coordBounds[aaCode_i];
    CoordNumber const &crd_bound_j = coordBounds[aaCode_j];

	Integers const &chainId = pseudoAtoms.chainId;

    // TODO: add neimin feature

	switch (type) {
	case BB:
        if (crd_i.backbone >= crd_bound_i.backbone) return true;
        if (crd_j.backbone >= crd_bound_j.backbone) return true;
		if ((abs(i - j) == 4) and chainId[i] == chainId[j]) return true;
		return false;

    // TODO: BS and SB don't check the coordination number for
    // backbone, is this a bug in cg.f?
	case BS:
        if (crd_j.sidechain >= crd_bound_j.sidechain) return true;
        return false;

	case SB:
        if (crd_i.sidechain >= crd_bound_i.sidechain) return true;
        return false;

	case SS:
        if (getSSBound(crd_i, aaCode_i) >= getSSBound(crd_bound_j, aaCode_i)) return true;
        if (getSSBound(crd_j, aaCode_j) >= getSSBound(crd_bound_i, aaCode_j)) return true;

        // TODO: add chainId condition here
        if (aaCode_i == TRP && aaCode_j == TRP && abs(i - j) == 3) return true; 

        // TODO: electrostatic interactions condition
        return false;

	default:
		return false;
	}
}

QuasiAdiabatic::ContactType QuasiAdiabatic::getContactType(int i, int j) const {
    Real3 const r_ij = (pseudoAtoms.pos.col(j) - pseudoAtoms.pos.col(i)).normalized();
    
    Real3 const h_i = hVector(i);
    Real3 const h_j = hVector(j);

    // Check for BB contact
    if (not isDisabled(BB, i, j)             and
        abs(h_i.dot(r_ij)) > backbone_2_min  and 
        abs(h_j.dot(r_ij)) > backbone_2_min  and
        abs(h_i.dot(h_j))  > backbone_1_min) {
        return BB;
    }

    // Check for SB contact
    Real3 const n_i = nVector(i);
    if (not isDisabled(SB, i, j)             and
        n_i.dot(r_ij) < sidechain_max        and
        abs(h_j.dot(r_ij)) > backbone_2_min) {
        return SB;
    }

	// Check for BS contact
    Real3 const n_j = nVector(j);
	if (not isDisabled(BS, i, j)             and
		n_j.dot(-r_ij) < sidechain_max       and
		abs(h_i.dot(r_ij)) > backbone_2_min) {
		return BS;
	}
		
	// Check for SS contact
	if (not isDisabled(SS, i, j)        and
		n_i.dot(r_ij)  < sidechain_max  and
		n_j.dot(-r_ij) < sidechain_max) {
		return SS;
	}

    return NONE;
}

Real QuasiAdiabatic::getCutoff(AACode i, AACode j, ContactType type) const {
    switch (type) {
    case SS:
        return ss_lj_r_min[i][j];
    case NONE:
        assert(false);
    default:
        return lj_r_min[type];
    };
}

Real QuasiAdiabatic::getMaxCutoff(AACode i, AACode j) const {
    static const Real max_r_min = max(lj_r_min[BB],
                                  max(lj_r_min[BS],
                                      lj_r_min[SB]));
    return max(max_r_min, ss_lj_r_min[i][j]);
}

void QuasiAdiabatic::compute(Reals3 &forces) {
    Reals3 const &pos = pseudoAtoms.pos;
    Integers const &chainId = pseudoAtoms.chainId;
    vector<AACode> const &aaCodes = pseudoAtoms.aminoAcidCode;

    for (auto const &[i, j] : verletList.pairs) {
        Real3 diff_vec = top.offset(pos.col(i), pos.col(j));
        Real sq_dist = diff_vec.squaredNorm();

        if (takenPairs.find({i, j}) != takenPairs.end()) continue;
        if (i == 0 || j == 0 || i + 1 == pseudoAtoms.n || j + 1 == pseudoAtoms.n) continue;
        if (chainId[i - 1] != chainId[i + 1] || chainId[j - 1] != chainId[j + 1]) continue;

        if (sq_dist >= totalCutoffSquared) continue;

        Real dist = sqrt(sq_dist);
        if (dist < cutoffCoeff * getMaxCutoff(aaCodes[i], aaCodes[j])) {

            ContactType type = getContactType(i, j);
            if (type != NONE && dist < cutoffCoeff * getCutoff(aaCodes[i], aaCodes[j], type)) {
                activeContacts.push_back(Contact {
                        .i = (unsigned) i,
                        .j = (unsigned) j,
                        .type = type,
                        .adiabCoeff = 0
                    });
                takenPairs.insert({i, j});
            }
        }
    }
}

void QuasiAdiabatic::dumpResults(Results &results) {
    results.potEnergy += energy;
    results.qaContacts = &activeContacts;
}
