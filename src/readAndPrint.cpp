#include "../include/readAndPrint.hpp"

int
readInput(ddpDomain_type & domain, ddpGrid_type & grid,
          ddpProblemInfo_type & problem,
          ddpCarrierConstants_type & carrierConstants,
          char const * nameOfFile)
{
    // COULD DO THIS A SHORTER WAY PASSING DOMAIN AND GRIDS
    // DATA TYPE ADDRESS BUT I THINK THIS IS A LITTLE SAFER/EASIER
    // TO READ.  Temporary memory anyway.

    // Physical inputs
    double
    xLeftEndPoint,
    xRightEndPoint,
    timeInitial,
    timeFinal,
    temperature,
    electronCharge,
    vacuumPermittivity,
    semiCondRelativePerm,
    electrolyteRelativePerm,
    BoltzmannConst,
    characteristicLength,
    characteristicTime,
    characteristicDensity,
    intrinsicDensity,
    appliedBias,
    absorptionCoeff,
    photonFlux,
    builtInBias,
    increase_time_step_factor;

    // electron and hole properties
    double
    electronMobility,
    electronRecomboTime,
    holeMobility,
    holeRecomboTime,
    reductantMobility,
    oxidantMobility,
    electronTransferRate,
    electronRecombinatioVelocity,
    holeTransferRate,
    holeRecombinatioVelocity;

    std::string
    ElecFieldCouplingStatus,
    IlluminationStatus,
    electronChargeSign,
    holeChargeSign,
    reductantChargeSign,
    oxidantChargeSign,
    ivp_type;

    // Problem dependent inputs
    int
    maxOrderMX,
    maxOrderDG,
    numElements,
    numTimeStamps,
    numBoundaryElements,
    GaussLegendreNumPoints;

    double boundaryLayerWidth;


    temperature    = 300;
    electronCharge = 1.6e-19;
    BoltzmannConst = 1.3792e-23;
    vacuumPermittivity = 8.85e-14;


///////////////////////////////////////////////////////////////////////////////
// READ IN INPUTS FROM BOOST INPUT PARSER
///////////////////////////////////////////////////////////////////////////////
    boost::property_tree::ptree ifl;
    boost::property_tree::ini_parser::read_ini(nameOfFile, ifl);


    // NOTE: If the values of the variables read in from file are the same as
    // the default values, grvy will print to screen that it is using the
    // pre-registered values
    // The title before the / and the variable is the section of the input file
    // Read from the computational constants section with default values given


    increase_time_step_factor = ifl.get<double>("computational.timeStepFactor");
    maxOrderMX = ifl.get<int>("computational.maxOrderMX");
    maxOrderDG = ifl.get<int>("computational.maxOrderDG");
    numElements = ifl.get<int>("computational.numElements");
    numBoundaryElements = ifl.get<int>("computational.numBoundaryElements");
    boundaryLayerWidth = ifl.get<double>("computational.boundaryLayerWidth");
    numTimeStamps = ifl.get<int>("computational.numTimeStamps");
    GaussLegendreNumPoints = ifl.get<int>("computational.GaussLegendreNumPoints");


    // Read from the Physical constants section with default values given
    xLeftEndPoint = ifl.get<double>("physical.xLeftEndPoint");
    xRightEndPoint = ifl.get<double>("physical.xRightEndPoint");
    timeInitial = ifl.get<double>("physical.timeInitial");
    timeFinal = ifl.get<double>("physical.timeFinal");
    appliedBias = ifl.get<double>("physical.appliedBias");
    semiCondRelativePerm = ifl.get<double>("physical.semiCondRelativePerm");
    electrolyteRelativePerm = ifl.get<double>("physical.electrolyteRelativePerm");
    characteristicLength = ifl.get<double>("physical.characteristicLength");
    characteristicTime = ifl.get<double>("physical.characteristicTime");
    characteristicDensity = ifl.get<double>("physical.characteristicDensity");
    absorptionCoeff = ifl.get<double>("physical.absorptionCoeff");
    photonFlux = ifl.get<double>("physical.photonFlux");
    intrinsicDensity = ifl.get<double>("physical.intrinsicDensity");
    builtInBias = ifl.get<double>("physical.builtInBias");


    // Read in the electron properties section with default values given
    electronMobility = ifl.get<double>("electrons.Mobility");
    electronRecomboTime = ifl.get<double>("electrons.recombinationTime");
    electronChargeSign = ifl.get<std::string>("electrons.ChargeSign");
    electronTransferRate = ifl.get<double>("electrons.TransferRate");
    electronRecombinatioVelocity = ifl.get<double>("electrons.RecombinationVelocity");

// Read in the electron properties section with default values given
    holeMobility = ifl.get<double>("holes.Mobility");
    holeRecomboTime = ifl.get<double>("holes.recombinationTime");
    holeChargeSign = ifl.get<std::string>("holes.ChargeSign");
    holeTransferRate = ifl.get<double>("holes.TransferRate");
    holeRecombinatioVelocity = ifl.get<double>("holes.RecombinationVelocity");
    // Read in the reductant properties section with default values given


    reductantMobility = ifl.get<double>("reductants.Mobility");
    reductantChargeSign = ifl.get<std::string>("reductants.ChargeSign");

// Read in the oxidant properties section with default values given
    oxidantMobility = ifl.get<double>("oxidants.Mobility");
    oxidantChargeSign = ifl.get<std::string>("oxidants.ChargeSign");


    // Read in the coupling status
    ElecFieldCouplingStatus = ifl.get<std::string>("couplingStatus.couplingToPoisson");

    // Read in the Illumination status
    IlluminationStatus = ifl.get<std::string>("illuminationStatus.illumination");


    // Some GSl functions aren't happy if GuassLegendreNumPoints is higher than 16
    assert(GaussLegendreNumPoints <= 16);

    //  The integration won't work properly if we don't ensure the following
    assert(2 * maxOrderDG <= GaussLegendreNumPoints);
    assert(2 * maxOrderMX <= GaussLegendreNumPoints);


///////////////////////////////////////////////////////////////////////////////
// Push everything on to grid and files.
///////////////////////////////////////////////////////////////////////////////

    // Prepare Grid

    // Set the SCALED grid
    domain.LeftEndPoint = xLeftEndPoint;
    domain.RightEndPoint = xRightEndPoint;

    // Assure number of boundary layer element dont exceed the number of elements
    assert(numElements >  numBoundaryElements);

    // Set the grids domain and set the grid info
// The + 1 accounts for the ghost element in the mixed FEM part for poisson.
    grid.Domain = domain;
    grid.OrderDGMax = maxOrderDG;
    grid.NumElementsNoGhost = numElements;
    grid.NumElementsWithGhost = numElements + 1;
    grid.NumBoundaryElements = numBoundaryElements;
    grid.BoundaryLayerWidth = boundaryLayerWidth;

    // Set problems constants from read in values

    // physical left end point
    problem.PhysicalStartEndPT = xLeftEndPoint * characteristicLength;

    // physical right end point
    problem.diffInEndPTS = (xRightEndPoint -xLeftEndPoint)
                           * characteristicLength;

    problem.characteristicLength = characteristicLength;

    problem.MaxOrderDG = maxOrderDG;
    problem.MaxOrderMX = maxOrderMX;
    problem.NumFrames = numTimeStamps;
    problem.GaussLegendreNumPoints = GaussLegendreNumPoints;
    problem.TimeFinal = timeFinal;
    problem.TimeInitial = timeInitial;
    problem.appliedBias = appliedBias;
    problem.temperature = temperature;
    problem.electronCharge = electronCharge;
    problem.vacuumPermittivity = vacuumPermittivity;
    problem.semiCondRelativePerm = semiCondRelativePerm;
    problem.electrolyteRelativePerm = electrolyteRelativePerm;
    problem.BoltzmannConst = BoltzmannConst;
    problem.thermalVoltage = BoltzmannConst*temperature/electronCharge;
    problem.characteristicTime = characteristicTime;
    problem.characteristicDensity = characteristicDensity;
    problem.intrinsicDensity = intrinsicDensity;
    problem.Absorption_Coeff = absorptionCoeff;
    problem.Photon_Flux = photonFlux;
    problem.BuiltInBias = builtInBias;
    problem.increase_time_step_factor = increase_time_step_factor;

    // Apply map to convert ElecFieldCoupingStatus from string
    // to ddpElecFielCoupling_type
    ddpElecFieldCouplingMap CouplingMap;

    CouplingMap.insert(std::pair<std::string,
                       ddpElecFieldCoupling_type>("On",On));

    CouplingMap.insert(std::pair<std::string,
                       ddpElecFieldCoupling_type>("Off",Off));

    problem.ElecFieldCouplingStatus = CouplingMap[ElecFieldCouplingStatus];

    // Apply map to convert charge sign from string to ddpChargeSign_type
    ddpChargeSignMap_type ChargeSignMap;

    ChargeSignMap.insert(std::pair<std::string,
                         ddpChargeSign_type>("Positive",Positive));

    ChargeSignMap.insert(std::pair<std::string,
                         ddpChargeSign_type>("Negative",Negative));


    // Apply map to convert illumination from string to
    // ddpIlluminationStatus_type
    ddpIlluminationMap IlluminationMap;
    IlluminationMap.insert(std::pair<std::string,
                           ddpIlluminationStatus_type>("On", Illuminated));

    IlluminationMap.insert(std::pair<std::string,
                           ddpIlluminationStatus_type>("Off", Dark));

    problem.IlluminationStatus = IlluminationMap[IlluminationStatus];


    // Assign read in values to carrierConsts
    carrierConstants.electron_Mobility = electronMobility;
    carrierConstants.electron_RecombinationTime = electronRecomboTime;
    carrierConstants.electron_ChargeSign = ChargeSignMap[electronChargeSign];
    carrierConstants.electron_TransferRate = electronTransferRate;
    carrierConstants.electron_RecombinationVelocity = electronRecombinatioVelocity;

    carrierConstants.hole_Mobility = holeMobility;
    carrierConstants.hole_RecombinationTime = holeRecomboTime;
    carrierConstants.hole_ChargeSign = ChargeSignMap[holeChargeSign];
    carrierConstants.hole_TransferRate = holeTransferRate;
    carrierConstants.hole_RecombinationVelocity = holeRecombinatioVelocity;

    carrierConstants.reductant_Mobility = reductantMobility;
    carrierConstants.reductant_ChargeSign = ChargeSignMap[reductantChargeSign];

    carrierConstants.oxidant_Mobility = oxidantMobility;
    carrierConstants.oxidant_ChargeSign = ChargeSignMap[oxidantChargeSign];

    return 0;
}


int ddpPrintGrid(ddpGrid_type const & grid)
{
    // This function will print out all the information about the grid;

    std::cout << "The grid has " << grid.NumElementsWithGhost
              << " elements.\n" << std::endl;
    std::cout << "The grid has a boundary layer of width W = " << grid.BoundaryLayerWidth
              << " at each end of the domain.\n" << std::endl;
    std::cout << "There are " << grid.NumBoundaryElements
              << " elements in each boundary layer. \n"  << std::endl;
    std::cout << "We list the information each element has: " << std::endl;



    for(int i = 0; i < grid.NumElementsWithGhost; ++i)
    {
        std::cout << "element[" << i << "]" << std::endl;
        std::cout << "The left endpoint of element[" << i << "] is "
                  << grid.ElementList[i].Left << std::endl;
        std::cout << "The right endpoint of element[" << i << "] is "
                  << grid.ElementList[i].Right << std::endl;
    }

    std::cout << "Remember last element is a ghost element" << std::endl;

    return 0;
}

int ddpPrintTimeStamps(std::vector< double > const & timeStamps)
{
    // prints the number of times stamps that were put into the
    // code using the input file and also will print out to the screen
    // an array of the actual times it will print
    std::cout << "We request a total of "  << timeStamps.size()
              << " snap shots to be taken " << std::endl;
    std::cout << "We request a \"time stamp\" be taken at the following times "
              << std::endl;

    for(int i = 0; i < timeStamps.size(); ++i)
    {
        std::cout << "time = " << timeStamps[i] << std::endl;
    }
    return 0;
}


int
ddpPrintBij_Backward(ddpBijBack_type const & bij)
{
    // This will print the backward bijection to the screen.
    // It will print the [element number, polynomial order]
    // and the global degree of freedom it gets mapped to
    //
    //
    // [elem, order] - > [index]


    int bijLength = bij.size();
    int globalCount;
    std::pair<int, int > orderedPair;
    std::cout << "bijLength = " << bijLength << std::endl;
    int a;
    int b;
    for(int i = 0; i < bijLength; ++i)
    {
        globalCount = i;
        orderedPair = bij.find(i)->second;
        a = orderedPair.first;
        b = orderedPair.second;
        std::cout << "The ordered pair ["<< a <<","<< b<< "] gets mapped to " << i
                  << ". Line number " << __LINE__  << std::endl;
    }
    return 0;
}




