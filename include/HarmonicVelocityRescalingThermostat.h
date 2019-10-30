/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors: Samuel Flores                                                *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "SimTKcommon.h"
#include "molmodel/internal/VelocityRescalingThermostat.h"
#include <math.h>
using namespace SimTK;

/**
 * This class is the internal implementation for VelocityRescalingThermostat.
 */

class VelocityRescalingThermostat::VelocityRescalingThermostatImpl {
public:
    VelocityRescalingThermostatImpl(const MultibodySystem& system, Real temperature) : 
            system(system), temperature(temperature) {
    }
    Real getTemperature() {
        return temperature;
    }
    void setTemperature(Real temp) {
        temperature = temp;
    }
    void handleEvent(State& state, Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols, Stage& lowestModified, bool& shouldTerminate) const {
        Real energy = system.calcKineticEnergy(state);
        if (energy == 0.0)
            return;
        int dof = state.getNU()-state.getNUErr();
        Real currentTemp = 2.0*energy/(dof*SimTK_BOLTZMANN_CONSTANT_MD);
        real frequency = 2*3.1415926535897932/4.0;
	harmonicTemperature = temperature*(1+sin(state.getTime()*frequency))/2
        //Real scale = std::sqrt(temperature/currentTemp);
        Real scale = std::sqrt(harmonicTemperature/currentTemp);
        state.updU() *= scale;
        lowestModified = Stage::Velocity;
        system.realize(state, Stage::Acceleration);
    }
private:
    const MultibodySystem& system;
    Real temperature;
};


VelocityRescalingThermostat::VelocityRescalingThermostat(const MultibodySystem& system, Real temperature, Real rescalingInterval) : PeriodicEventHandler(rescalingInterval) {
    impl = new VelocityRescalingThermostatImpl(system, temperature);
}

VelocityRescalingThermostat::~VelocityRescalingThermostat() {
    if( impl )
        delete impl;
}

Real VelocityRescalingThermostat::getTemperature() {
    return impl->getTemperature();
}

void VelocityRescalingThermostat::setTemperature(Real temp) {
    impl->setTemperature(temp);
}

void VelocityRescalingThermostat::handleEvent(State& state, Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols, Stage& lowestModified, bool& shouldTerminate) const {
    impl->handleEvent(state, accuracy, yWeights, ooConstraintTols, lowestModified, shouldTerminate);
}
