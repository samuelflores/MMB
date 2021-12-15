// vim: expandtab ts=4 sts=4 sw=4 :
/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */
#include "NtCForces.h"
#include <SimTKcommon/Scalar.h>
#include <SimTKcommon/internal/SpatialAlgebra.h>
#include <cstddef>
#include <array>
#include <ios>
#include <cmath>
#include "MMBLogger.h"
#include "ParameterReader.h"

#define _R2D(x) ((x) * (180.0 / SimTK::Pi))
#define _D2R(x) ((x) * (SimTK::Pi / 180.0))
#define _2PI (2.0 * SimTK::Pi)

static const double BF_SIGN[4] = { 1.0, 1.0, -1.0, -1.0 };

#ifdef NTC_DEBUG_CALC
std::ofstream ntcdump{};
size_t dump_counter{0};
#endif // NTC_DEBUG_CALC

template <typename T>
int int_round(const T &v) {
    T vs = v + T(0.5) * (v > T(0)) - T(0.5) * (v < T(0));
    return int(vs);
}

static
double angle_rad(const Vec3 &cross_1, const Vec3 &cross_2, const Vec3 &cross_3, const Vec3 &d_d2) {
    Vec3 direction;
    direction[0] = cross_3[0] * d_d2[0];
    direction[1] = cross_3[1] * d_d2[1];
    direction[2] = cross_3[2] * d_d2[2];

    double scal = dot(cross_1, cross_2);
    scal = scal > 1.0 ? 1.0 : (scal < -1.0 ? -1.0 : scal);

    double angle = std::acos(scal);
    if (direction[0] < 0.0 && direction[1] < 0.0 && direction[2] < 0.0)
        angle = -angle;
    angle += _2PI * (angle < 0.0);

    return angle;
}

static
double dist_angle_rad(double angle, double rotationAngle) {
    double dist = fabs(angle - rotationAngle);
    if (dist > SimTK::Pi)
        dist = _2PI - dist;

    const double from = rotationAngle;
    const double to = from > SimTK::Pi ? from - SimTK::Pi : from + SimTK::Pi;

    double flip = 0.0;
    if (from < to) {
        flip = (from < angle && angle <= to) ? 1.0 : -1.0;
    } else {
        flip = (to < angle && angle <= from) ? -1.0 : 1.0;
    }

    return flip * dist;
}

#ifdef NTC_DEBUG_CALC
static
auto operator<<(std::ostream &os, const NTC_Classes &ntc) -> std::ostream & {
    os << "[" << ntc.NtC_Class_String << ", " << ntc.FirstBPResidue.outString() << ", " << ntc.SecondBPResidue.outString() << ", " << ntc.NtC_step_ID << "]";
    return os;
}

template <typename T, size_t N>
auto operator<<(std::ostream &os, const std::array<T, N> &arr) -> std::ostream & {
    os << "[";
    for (size_t idx{0}; idx < N - 1; idx++)
        os << arr[idx] << ", ";
    os << arr[N - 1] << "]";
    return os;
}

template <typename Arg>
auto _dump(std::ostream &os, Arg &&a) {
    os << a << "\n";
}

template <typename Arg, typename ...Args>
auto _dump(std::ostream &os, Arg &&a, Args &&...args) {
    os << a << "; ";
    _dump(os, std::forward<Args>(args)...);
}
#endif // NTC_DEBUG_CALC

NTC_Torque::NTC_Torque(SimbodyMatterSubsystem &matter,
                       ParameterReader &myParameterReader,
                       NTC_PAR_Class &myNTC_PAR_Class,
                       BiopolymerClassContainer &myBiopolymerClassContainer)
    : matter(matter), myParameterReader(myParameterReader),
      myNTC_PAR_Class(myNTC_PAR_Class),
      myBiopolymerClassContainer(myBiopolymerClassContainer)
{}


void NTC_Torque::calcForce(const State &state, Vector_<SpatialVec> &bodyForces,
                           Vector_<Vec3> &particleForces,
                           Vector &mobilityForces) const {

#ifdef NTC_DEBUG_CALC
    if (!ntcdump.is_open())
        ntcdump.open("ntcdump.txt");
#endif // NTC_DEBUG_CALC

    const auto SCALE_FACTOR = myParameterReader.NtCForceScaleFactor;
    Vec3 states[4];
    for (int r = 0; r < myParameterReader.ntc_class_container.numNTC_Torsions(); r++) {
        const auto &ntc = myParameterReader.ntc_class_container.getNTC_Class(r);
        const auto &bondRow = myNTC_PAR_Class.myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[ntc.NTC_PAR_BondRowIndex];
        const auto &indices = ntc.atomIndices;
        auto bp = &myBiopolymerClassContainer.updBiopolymerClass(ntc.NtC_FirstBPChain);

        if (bondRow.bondLength[0] == 0.0) {
            for (std::size_t idx = 0; idx < 4; idx++) {
                states[idx] = bp->calcAtomLocationInGroundFrame(state, indices[idx]);
            }

            Vec3 d_d1 = states[1] - states[0];
            Vec3 d_d2 = states[2] - states[1];
            Vec3 d_d3 = states[3] - states[2];

            Vec3 cross_1 = d_d1 % d_d2;
            Vec3 cross_2 = d_d2 % d_d3;

            cross_1 = cross_1 / cross_1.norm();
            cross_2 = cross_2 / cross_2.norm();

            Vec3 cross_3 = cross_1 % cross_2;

            double angle = angle_rad(cross_1, cross_2, cross_3, d_d2);
            double dist_ang = dist_angle_rad(angle, bondRow.rotationAngle);

            static const std::set<std::string> consideredTorsions = {
                "tau0a", "tau4a", "tau0b", "tau4b",
                "delta", "epsilon", "zeta", "alpha1", "beta1", "gamma1", "delta1",
                "chi", "chi1"
            };

            if (ntc.meta == 0) {
                double force = 0;
                auto calculator = consideredTorsions.find(bondRow.dihedraltype);
                if (calculator != consideredTorsions.cend())
                    force = SCALE_FACTOR * dist_ang;

                Vec3 torque = d_d2 / d_d2.norm() * force;

#ifdef NTC_DEBUG_CALC
                if (dump_counter % 10 == 0)
                    _dump(ntcdump, r, bondRow.dihedraltype, _R2D(angle), dist_ang, (100 * dist_ang / SimTK::Pi), force, torque.real());
#endif // NTC_DEBUG_CALC

                for (std::size_t idx = 0; idx < 4; idx++) {
                    const auto &body = bp->updAtomMobilizedBody(matter, indices[idx]);
                    bodyForces[body.getMobilizedBodyIndex()] += BF_SIGN[idx] * SpatialVec(torque, Vec3(0));
                }
            } else if (ntc.meta == 1) {
                MMBLOG_PLAIN(CRITICAL, "NtC Meta mode is unavailable in this build of MMB");
                /*
                double bias = 0;
                int value = -1;

                angle *= K_ANGLE;

                if (isfinite(angle) == 1) {
                    int i = int_round(angle);
                    const auto count = ntc.count;

                    myBiopolymerClassContainer.hist[count][i] += 1.0;
                    value = i;

                    myBiopolymerClassContainer.counter[count] += 1.0;

                    if (myBiopolymerClassContainer.counter[count] > 0.0) {
                        myBiopolymerClassContainer.prob[count][i] =
                        myBiopolymerClassContainer.hist[count][i] / myBiopolymerClassContainer.counter[count];
                    }

                    if (value > 0 &&
                        myBiopolymerClassContainer.prob[count][value] > 1e-3 &&
                        myBiopolymerClassContainer.prob[count][value + 1] > 1e-3 &&
                        value < 360) {
                        bias = 2.479 * (log(myBiopolymerClassContainer.prob[count][value] / (1e-3)) -
                                        log(myBiopolymerClassContainer.prob[count][value + 1] / (1e-3)));
                    }

                    if (myBiopolymerClassContainer.prob[count][value + 1] > 1e-3 &&
                        value == 360) {
                        bias = 2.479 * (log(myBiopolymerClassContainer.prob[count][value] / (1e-3)) -
                                        log(myBiopolymerClassContainer.prob[count][1] /(1e-3)));
                    }

                    double pot_angle = torqueConstant * (-sin((dist_ang + 180.0) / K_ANGLE)) *
                                       (360.0 / K_ANGLE + 1.0) /
                                       (1.0 + bondRow.CONFALVALUE) / (360.0 / K_ANGLE);

                    Vec3 torque = d_d2 / d_d2.norm() * (pot_angle) / (1.0 + ntc.weight2) * ntc.weight;

                    for (std::size_t idx = 0; idx < 4; idx++) {
                        const auto &body = bp->updAtomMobilizedBody(matter, indices[idx]);
                        bodyForces[body.getMobilizedBodyIndex()] += BF_SIGN[idx] * SpatialVec(torque, Vec3(0));
                    }

                    if (isfinite(log(myBiopolymerClassContainer.prob[count][value] / (1e-3))) &&
                        isfinite(log(myBiopolymerClassContainer.prob[count][value + 1] / (1e-3)))) {

                        if (isfinite(bias) && bias != 0.0) {
                            const auto absBias = std::abs(bias);
                            Vec3 torque = -d_d2 / d_d2.norm() * bias * std::abs(pot_angle) / absBias * ntc.weight2 * ntc.weight;

                            for (std::size_t idx = 0; idx < 4; idx++) {
                                const auto &body = bp->updAtomMobilizedBody(matter, indices[idx]);
                                bodyForces[body.getMobilizedBodyIndex()] += BF_SIGN[idx] * SpatialVec(torque, Vec3(0));
                            }
                        }
                    }
                }
                */
            }
            // end real torsions
        } else { // bonds
            // No bond lengths are calculated
        }
    }
}

Real NTC_Torque::calcPotentialEnergy(const State &state) const {
    double energy = 0.0;
    double rms = 0.0;
    double rmsTorsionAngleForThisNtCAndDinucleotide = 0.0;
    String oldNtCClassString = "ZZZZZZ";

    const auto SCALE_FACTOR = myParameterReader.NtCForceScaleFactor;
    Vec3 states[4];
    for (int r = 0; r < myParameterReader.ntc_class_container.numNTC_Torsions(); r++) {
        const auto &ntc = myParameterReader.ntc_class_container.getNTC_Class(r);
        const auto &indices = ntc.atomIndices;
        auto bp = &myBiopolymerClassContainer.updBiopolymerClass(ntc.NtC_FirstBPChain);
        // If we have changed our NtC class type, meaning we are computing a new NtC
        if (ntc.NtC_Class_String != oldNtCClassString) {
            if (r > 0) {
                MMBLOG_FILE_FUNC_LINE(
                    INFO, "RMSD Angle sum for NtC of type "
                      << ntc.NtC_step_ID
                      << ntc.NtC_INDEX
                      << " " << oldNtCClassString << " is "
                      << rmsTorsionAngleForThisNtCAndDinucleotide << endl);
            }

            rmsTorsionAngleForThisNtCAndDinucleotide = 0.;
        }

        const auto &chainId1 = ntc.NtC_FirstBPChain;
        const auto &bondRow = myNTC_PAR_Class.myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[ntc.NTC_PAR_BondRowIndex];

        MMBLOG_FILE_FUNC_LINE(DEBUG, "doing base pair #" << r << endl);

        if (bondRow.bondLength[0] == 0.0) {
            for (std::size_t idx = 0; idx < 4; idx++) {
                states[idx] = bp->calcAtomLocationInGroundFrame(state, indices[idx]);
            }

            Vec3 d_d1 = states[1] - states[0];
            Vec3 d_d2 = states[2] - states[1];
            Vec3 d_d3 = states[3] - states[2];

            Vec3 cross_1 = d_d1 % d_d2;
            Vec3 cross_2 = d_d2 % d_d3;

            cross_1 = cross_1 / cross_1.norm();
            cross_2 = cross_2 / cross_2.norm();

            Vec3 cross_3 = cross_1 % cross_2;

            double angle = return_angle(cross_1, cross_2, cross_3, d_d2);
            double dist_ang = return_dist_ang(angle, bondRow.rotationAngle);

            double e = 0.5 * SCALE_FACTOR * dist_ang * dist_ang;

            energy += e;

            MMBLOG_FILE_FUNC_LINE(
                DEBUG, " NTC sampling - CHAIN ID = "
                       << chainId1 << ", residuenumbers " << ntc.FirstBPResidue.ResidueNumber << "/" << ntc.SecondBPResidue.ResidueNumber
                       << " difference-angle = " << dist_ang
                       << " , CONFALVALUE = " << bondRow.CONFALVALUE << " , "
                       << _R2D(angle) << " = angle at time t for atoms  = "
                       << bondRow.residue1Atom[0] << " , "
                       << bondRow.residue1Atom[1] << " , "
                       << bondRow.residue1Atom[2] << " , "
                       << bondRow.residue1Atom[3] << " , "
                       << _R2D(bondRow.rotationAngle)
                       << " = angle_0 from  input , "
                       << "energy = " << energy << endl);

            auto daSq = std::abs(dist_ang);
            rms += daSq;
            rmsTorsionAngleForThisNtCAndDinucleotide += daSq;
          // end real torsions
        } else { // bonds
            // Bond lengths are disregarded
        }

        oldNtCClassString = ntc.NtC_Class_String;
    }

    MMBLOG_PLAIN_NOSEV(INFO, rms << " = RMSD Angle sum " << endl);

    return energy;
}

bool NTC_Torque::dependsOnlyOnPositions() const { return true; }

Real NTC_Torque::return_dist_ang(double angle, double rotationAngle) const {
    double ang_diff = _R2D(angle - rotationAngle);
    double dist_ang = 180.0 - abs(180.0 - abs(ang_diff));

    double angle_1 = _R2D(angle);
    double angle_2 = _R2D(rotationAngle);

    double interval_begin = angle_2;
    double interval_end = (interval_begin + 180.0);
    if (interval_end > 360.0)
        interval_end -= 360.0;

    if (interval_end > interval_begin) {
        if (angle_1 < interval_begin || angle_1 > interval_end) {
            dist_ang = -dist_ang;
        }
    } else {
        if (angle_1 < interval_begin && angle_1 > interval_end) {
            dist_ang = -dist_ang;
        }
    }

    return dist_ang;
}

Real NTC_Torque::return_angle(const Vec3 &cross_1, const Vec3 &cross_2, const Vec3 &cross_3, const Vec3 &d_d2) const {
    double angle;

    Vec3 direction;

    direction[0] = cross_3[0] * d_d2[0];
    direction[1] = cross_3[1] * d_d2[1];
    direction[2] = cross_3[2] * d_d2[2];

    double scalar_product = dot(cross_1, cross_2);

    if (scalar_product > 1.0)
        scalar_product = 1.0;
    else if (scalar_product < -1.0)
        scalar_product = -1.0;

    angle = std::acos(scalar_product) * 180.0 / SimTK::Pi;
    double flipper = 1.0 - 2.0 * (direction[0] < 0.0 && direction[1] < 0.0 && direction[2] < 0.0);
    angle = flipper * angle;

    angle += 360.0 * (angle < 0.0);

    return _D2R(angle);
}
