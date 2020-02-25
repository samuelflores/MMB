
%pythonprepend OpenMM::MonteCarloBarostat::setDefaultPressure(double pressure) %{
    if unit.is_quantity(pressure):
        pressure = pressure.value_in_unit(unit.bar)
%}
%pythonprepend OpenMM::MonteCarloBarostat::setDefaultTemperature(double temp) %{
    if unit.is_quantity(temp):
        temp = temp.value_in_unit(unit.kelvin)
%}
%pythonprepend OpenMM::MonteCarloMembraneBarostat::MonteCarloMembraneBarostat %{
    args = list(args)
    if unit.is_quantity(args[0]):
        args[0] = args[0].value_in_unit(unit.bar)
    if unit.is_quantity(args[1]):
        args[1] = args[1].value_in_unit(unit.bar*unit.nanometer)
    if unit.is_quantity(args[2]):
        args[2] = args[2].value_in_unit(unit.kelvin)
%}
%pythonprepend OpenMM::MonteCarloMembraneBarostat::setDefaultPressure(double pressure) %{
    if unit.is_quantity(pressure):
        pressure = pressure.value_in_unit(unit.bar)
%}
%pythonprepend OpenMM::MonteCarloMembraneBarostat::setDefaultSurfaceTension(double surfaceTension) %{
    if unit.is_quantity(surfaceTension):
        surfaceTension = surfaceTension.value_in_unit(unit.bar*unit.nanometer)
%}
%pythonprepend OpenMM::MonteCarloMembraneBarostat::setDefaultTemperature(double temp) %{
    if unit.is_quantity(temp):
        temp = temp.value_in_unit(unit.kelvin)
%}
%pythonprepend OpenMM::AndersenThermostat::setDefaultTemperature(double temperature) %{
    if unit.is_quantity(temperature):
        temperature = temperature.value_in_unit(unit.kelvin)
%}
%pythonprepend OpenMM::CustomCVForce::addCollectiveVariable(const std::string &name, Force *variable) %{
    if not variable.thisown:
        s = ("the %s object does not own its corresponding OpenMM object"
             % self.__class__.__name__)
        raise Exception(s)
%}
%pythonprepend OpenMM::CustomCVForce::addTabulatedFunction(const std::string &name, TabulatedFunction *function) %{
    if not function.thisown:
        s = ("the %s object does not own its corresponding OpenMM object"
             % self.__class__.__name__)
        raise Exception(s)
%}
%pythonprepend OpenMM::CustomGBForce::addTabulatedFunction(const std::string &name, TabulatedFunction *function) %{
    if not function.thisown:
        s = ("the %s object does not own its corresponding OpenMM object"
             % self.__class__.__name__)
        raise Exception(s)
%}
%pythonprepend OpenMM::NoseHooverIntegrator::addThermostat(double temperature, double collisionFrequency, int chainLength, int numMTS, int numYoshidaSuzuki) %{
    if unit.is_quantity(temperature):
        temperature = temperature.value_in_unit(unit.kelvin)
    if unit.is_quantity(collisionFrequency):
        collisionFrequency = collisionFrequency.value_in_unit(unit.picosecond**-1)
%}
%pythonprepend OpenMM::NoseHooverIntegrator::addSubsystemThermostat(const std::vector< int > &thermostatedParticles, const std::vector< std::pair< int, int > > &thermostatedPairs, double temperature, double collisionFrequency, double relativeTemperature, double relativeCollisionFrequency, int chainLength=3, int numMTS=3, int numYoshidaSuzuki=3) %{
    if unit.is_quantity(temperature):
        temperature = temperature.value_in_unit(unit.kelvin)
    if unit.is_quantity(collisionFrequency):
        collisionFrequency = collisionFrequency.value_in_unit(unit.picosecond**-1)
    if unit.is_quantity(relativeTemperature):
        relativeTemperature = relativeTemperature.value_in_unit(unit.kelvin)
    if unit.is_quantity(relativeCollisionFrequency):
        relativeCollisionFrequency = relativeCollisionFrequency.value_in_unit(unit.picosecond**-1)
%}
%pythonprepend OpenMM::NoseHooverIntegrator::setTemperature(double temperature, int chainID=0) %{
    if unit.is_quantity(temperature):
        temperature = temperature.value_in_unit(unit.kelvin)
%}
%pythonprepend OpenMM::NoseHooverIntegrator::setRelativeTemperature(double temperature, int chainID=0) %{
    if unit.is_quantity(temperature):
        temperature = temperature.value_in_unit(unit.kelvin)
%}
%pythonprepend OpenMM::NoseHooverIntegrator::setCollisionFrequency(double frequency, int chainID=0) %{
    if unit.is_quantity(frequency):
        frequency = frequency.value_in_unit(unit.picosecond**-1)
%}
%pythonprepend OpenMM::NoseHooverIntegrator::setRelativeCollisionFrequency(double frequency, int chainID=0) %{
    if unit.is_quantity(frequency):
        frequency = frequency.value_in_unit(unit.picosecond**-1)
%}
%pythonprepend OpenMM::NoseHooverIntegrator::setMaximumPairDistance(double distance) %{
    if unit.is_quantity(distance):
        distance = distance.value_in_unit(unit.nanometer)
%}
%pythonprepend OpenMM::DrudeNoseHooverIntegrator::setMaxDrudeDistance(double distance) %{
    if unit.is_quantity(distance):
        distance = distance.value_in_unit(unit.nanometer)
%}
%pythonprepend OpenMM::RPMDMonteCarloBarostat::setDefaultPressure(double pressure) %{
    if unit.is_quantity(pressure):
        pressure = pressure.value_in_unit(unit.bar)
%}
%pythonprepend OpenMM::CustomManyParticleForce::addTabulatedFunction(const std::string &name, TabulatedFunction *function) %{
    if not function.thisown:
        s = ("the %s object does not own its corresponding OpenMM object"
             % self.__class__.__name__)
        raise Exception(s)
%}
%pythonprepend OpenMM::CompoundIntegrator::addIntegrator(Integrator *integrator) %{
    if not integrator.thisown:
        s = ("the %s object does not own its corresponding OpenMM object"
             % self.__class__.__name__)
        raise Exception(s)
%}
%pythonprepend OpenMM::MonteCarloAnisotropicBarostat::setDefaultPressure(const Vec3 &pressure) %{
    if unit.is_quantity(pressure):
        pressure = pressure.value_in_unit(unit.bar)
%}
%pythonprepend OpenMM::MonteCarloAnisotropicBarostat::setDefaultTemperature(double temp) %{
    if unit.is_quantity(temp):
        temp = temp.value_in_unit(unit.kelvin)
%}
%pythonprepend OpenMM::CustomHbondForce::addTabulatedFunction(const std::string &name, TabulatedFunction *function) %{
    if not function.thisown:
        s = ("the %s object does not own its corresponding OpenMM object"
             % self.__class__.__name__)
        raise Exception(s)
%}
%pythonprepend OpenMM::CustomIntegrator::addTabulatedFunction(const std::string &name, TabulatedFunction *function) %{
    if not function.thisown:
        s = ("the %s object does not own its corresponding OpenMM object"
             % self.__class__.__name__)
        raise Exception(s)
%}
%pythonprepend OpenMM::CustomNonbondedForce::addTabulatedFunction(const std::string &name, TabulatedFunction *function) %{
    if not function.thisown:
        s = ("the %s object does not own its corresponding OpenMM object"
             % self.__class__.__name__)
        raise Exception(s)
%}
%pythonprepend OpenMM::CustomNonbondedForce::addInteractionGroup(const std::set< int > &set1, const std::set< int > &set2) %{
    set1 = list(set1)
    set2 = list(set2)
%}
%pythonprepend OpenMM::CustomNonbondedForce::setInteractionGroupParameters(int index, const std::set< int > &set1, const std::set< int > &set2) %{
    set1 = list(set1)
    set2 = list(set2)
%}
%pythonprepend OpenMM::System::setVirtualSite(int index, VirtualSite *virtualSite) %{
    if not virtualSite.thisown:
        s = ("the %s object does not own its corresponding OpenMM object"
             % self.__class__.__name__)
        raise Exception(s)
%}
%pythonprepend OpenMM::System::addForce(Force *force) %{
    if not force.thisown:
        s = ("the %s object does not own its corresponding OpenMM object"
             % self.__class__.__name__)
        raise Exception(s)
%}
%pythonprepend OpenMM::Platform::registerPlatform(Platform *platform) %{
    if not platform.thisown:
        s = ("the %s object does not own its corresponding OpenMM object"
             % self.__class__.__name__)
        raise Exception(s)
%}
%pythonprepend OpenMM::CustomCompoundBondForce::addTabulatedFunction(const std::string &name, TabulatedFunction *function) %{
    if not function.thisown:
        s = ("the %s object does not own its corresponding OpenMM object"
             % self.__class__.__name__)
        raise Exception(s)
%}