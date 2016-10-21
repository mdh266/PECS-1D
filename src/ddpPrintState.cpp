#include "../include/ddpPrintState.hpp"

///////////////////////////////////////////////////////////////////////////////
// Unipolar Printing
///////////////////////////////////////////////////////////////////////////////


int
ddpPrintState(ddpProblemInfo_type const & problem,
              ddpCarrier_type const & carrier1,
              ddpPoisson_type const & poisson,
              int const & timeStamp
             )
{
    // THIS WILL PRINT THE STATE OF THE SIMULATION AT TIME STAMPS
    // THE FILE WILL BE NAMED "State(TimeStamp#).dat" and will contain
    // Space delimited values:
    // pt uPTval qPTval ElectricFieldPTVal PotentialPTVal JPTval

    // this is to deal with the name of the file
    ofstream ptr;
    std::string prefix = "State";
    std::string extension = ".dat";
    std::string filename;
    std::stringstream ss;

    char snapShotNumberString[3];

    sprintf (snapShotNumberString, "%4.4d", timeStamp);

    // Use stringstream toconcatenate the strings using C++ including '\0'
    // character to signify end of string.
    ss << prefix << snapShotNumberString << extension << '\0';

    // convert the string stream into a C++ string and then int a C string.
    filename = ss.str ();
    char str[filename.size ()];

    // convert the C++ string to a C string by copying by value.
    for (unsigned int i = 0; i < filename.size (); ++i)
    {
        str[i] = filename[i];
    }

    // Now can finally use the string str = "uUvals(snapShotNumber).txt\0"
    ptr.open (str);


    ddpDenseVector_type
    XPTS = carrier1.PTSDense,
    UPTS = carrier1.VandeMondeMatrices.globalVandeMondeDG
           * carrier1.carrierState.uDof,
           QPTS = carrier1.VandeMondeMatrices.globalVandeMondeDG
                  * carrier1.carrierState.qDof,
                  ElecPTS = carrier1.VandeMondeMatrices.globalVandeMondeMX
                            * poisson.PoissonState.elecDof,
                            POTPTS = carrier1.VandeMondeMatrices.globalVandeMondeDG
                                     * poisson.PoissonState.potDof;

    double Sign = carrier1.carrierProps.Sign4Force;

/// MAKE SURE EVERYTHING IS THE SAME LENGTH
    assert(XPTS.size() == UPTS.size());
    assert(UPTS.size() == QPTS.size());
    assert(QPTS.size() == ElecPTS.size());
    assert(ElecPTS.size() == POTPTS.size());

    double x_scale = problem.characteristicLength;
    double t_scale = problem.characteristicTime;
    double den_scale = problem.characteristicDensity;
    double pot_scale = problem.thermalVoltage;
    double elec_scale = pot_scale/x_scale;
    double curr_scale = problem.electronCharge * den_scale * x_scale / t_scale;

    //cout << curr_scale << endl;
    // PRINTS THE LINE (SPACE DELIMTED)
    // pt uPTval qPTval ElectricFieldPTVal PotentialPTVal JPTval

    for(int i = 0; i < UPTS.size(); ++i)
    {
        ptr << x_scale * XPTS(i) /*+ problem.startEndPT */  << " "
            << den_scale * UPTS(i) << " "
            << elec_scale * ElecPTS(i) << " "
            << pot_scale * POTPTS(i) << " "
            << curr_scale * QPTS(i) << endl;
    }
    ptr.close();

    return 0;
}

int
ddpPrintState(ddpProblemInfo_type const & problem,
              ddpCarrier_type  & electrons,
              ddpCarrier_type  & holes,
              ddpPoisson_type const & poisson,
              int const & timeStamp,
              double const & tCurrent
             )
{
    // THIS WILL PRINT THE STATE OF THE SIMULATION AT TIME STAMPS
    // THE FILE WILL BE NAMED "State(TimeStamp#).dat" and will contain
    // Space delimited values:
    // pt uPTval qPTval ElectricFieldPTVal PotentialPTVal JPTval

    // this is to deal with the name of the file
    ofstream ptr;
    std::string prefix = "State";
    std::string extension = ".dat";
    std::string filename;
    std::stringstream ss;

    char snapShotNumberString[3];

    sprintf (snapShotNumberString, "%4.4d", timeStamp);

    // Use stringstream toconcatenate the strings using C++ including '\0'
    // character to signify end of string.
    ss << prefix << snapShotNumberString << extension << '\0';

    // convert the string stream into a C++ string and then int a C string.
    filename = ss.str ();
    char str[filename.size ()];

    // convert the C++ string to a C string by copying by value.
    for (unsigned int i = 0; i < filename.size (); ++i)
    {
        str[i] = filename[i];
    }

    // Now can finally use the string str = "uUvals(snapShotNumber).txt\0"
    ptr.open (str);

    double x_scale = problem.characteristicLength;
    double t_scale = problem.characteristicTime;
    double den_scale = problem.characteristicDensity;
    double pot_scale = problem.thermalVoltage;
    double elec_scale = pot_scale/x_scale;
    double curr_scale = problem.electronCharge * den_scale * x_scale / t_scale;

    ddpDenseVector_type
    XPTS = electrons.PTSDense,
    ElectronPTS = electrons.VandeMondeMatrices.globalVandeMondeDG
                  * electrons.carrierState.uDof,
                  HolePTS = holes.VandeMondeMatrices.globalVandeMondeDG
                            * holes.carrierState.uDof,
                            ElectronsQPTS = electrons.VandeMondeMatrices.globalVandeMondeDG
                                            * electrons.carrierState.qDof,
                                            HolesQPTS = holes.VandeMondeMatrices.globalVandeMondeDG
                                                    * holes.carrierState.qDof,
                                                    ElecPTS = electrons.VandeMondeMatrices.globalVandeMondeMX
                                                            * poisson.PoissonState.elecDof,
                                                            POTPTS = electrons.VandeMondeMatrices.globalVandeMondeDG
                                                                    * poisson.PoissonState.potDof;

    int N = XPTS.size();

    // Prints the carriers on the semiconductor side
    for(int i = 0; i < N; ++i)
    {
        ptr << x_scale * XPTS(i)  << " "
            << den_scale * ElectronPTS(i) << " "
            << den_scale * HolePTS(i) << " "
            << elec_scale * ElecPTS(i) << " "
            << pot_scale * POTPTS(i) << " "
            << curr_scale * (ElectronsQPTS(i) - HolesQPTS(i)) << " "
            << t_scale * tCurrent
            << endl;
    }
    ptr.close();

    return 0;
}





///////////////////////////////////////////////////////////////////////////////
// Semiconductor-Electrolyte Interface TEST Printing
///////////////////////////////////////////////////////////////////////////////

int
ddpPrintStateTest(ddpProblemInfo_type const & problem,
                  ddpCarrier_type  & carrier,
                  ddpCarrier_type  & carrier1,
                  ddpCarrier_type  & carrier2,
                  ddpPoisson_type const & poisson,
                  int const & timeStamp
                 )
{
    // THIS WILL PRINT THE STATE OF THE SIMULATION AT TIME STAMPS
    // THE FILE WILL BE NAMED "State(TimeStamp#).dat" and will contain
    // Space delimited values:
    // pt uPTval qPTval ElectricFieldPTVal PotentialPTVal JPTval

    // this is to deal with the name of the file
    ofstream ptr;
    std::string prefix = "State";
    std::string extension = ".dat";
    std::string filename;
    std::stringstream ss;

    char snapShotNumberString[3];

    sprintf (snapShotNumberString, "%4.4d", timeStamp);

    // Use stringstream toconcatenate the strings using C++ including '\0'
    // character to signify end of string.
    ss << prefix << snapShotNumberString << extension << '\0';

    // convert the string stream into a C++ string and then int a C string.
    filename = ss.str ();
    char str[filename.size ()];

    // convert the C++ string to a C string by copying by value.
    for (unsigned int i = 0; i < filename.size (); ++i)
    {
        str[i] = filename[i];
    }

    // Now can finally use the string str = "uUvals(snapShotNumber).txt\0"
    ptr.open (str);


    ddpDenseVector_type
    XPTS1 = carrier.PTSDense,
    XPTS2 = carrier2.PTSDense,
    UPTS = carrier.VandeMondeMatrices.globalVandeMondeDG
           * carrier.carrierState.uDof,
           QPTS = carrier.VandeMondeMatrices.globalVandeMondeDG
                  * carrier.carrierState.qDof,
                  U1PTS = carrier1.VandeMondeMatrices.globalVandeMondeDG
                          * carrier1.carrierState.uDof,
                          Q1PTS = carrier1.VandeMondeMatrices.globalVandeMondeDG
                                  * carrier1.carrierState.qDof,
                                  U2PTS = carrier2.VandeMondeMatrices.globalVandeMondeDG
                                          * carrier2.carrierState.uDof,
                                          Q2PTS = carrier2.VandeMondeMatrices.globalVandeMondeDG
                                                  * carrier2.carrierState.qDof;

    double Sign = carrier.carrierProps.Sign4Force;
    double Sign1 = carrier1.carrierProps.Sign4Force;
    double Sign2 = carrier2.carrierProps.Sign4Force;

    int N = XPTS1.size();

    // Prints the carriers on the semiconductor side
    for(int i = 0; i < N; ++i)
    {
        ptr << XPTS1(i)  << " "
            << UPTS(i) << " "
            << 0.0 << " "
            << 0.0 << " "
            << 0.0 << " "
            << QPTS(i) << " "
            << 0.0 << " "
            << 0.0
            << endl;
    }
    // Prints the carriers on the electrolyte side
    for(int i = 0; i < N; i++)
    {
        ptr << XPTS2(i) << " "
            << U1PTS(i) << " "
            << U2PTS(i) << " "
            << 0.0 << " "
            << 0.0 << " "
            << -Q2PTS(i) << " "
            << 0.0 << " "
            << 0.0
            << endl;
    }
    ptr.close();

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
// Semiconductor-Electrolyte Interface Model Printing
///////////////////////////////////////////////////////////////////////////////
int
ddpPrintState(ddpProblemInfo_type const & problem,
              ddpCarrier_type const & electrons,
              ddpCarrier_type const & holes,
              ddpCarrier_type const & reductants,
              ddpCarrier_type const & oxidants,
              ddpPoisson_type const & poisson,
              int const & timeStamp,
              double const & tCurrent)
{
    // THIS WILL PRINT THE STATE OF THE SIMULATION AT TIME STAMPS
    // THE FILE WILL BE NAMED "State(TimeStamp#).dat" and will contain
    // Space delimited values:
    // pt uPTval qPTval ElectricFieldPTVal PotentialPTVal JPTval

    // this is to deal with the name of the file
    ofstream ptr;
    std::string prefix = "State";
    std::string extension = ".dat";
    std::string filename;
    std::stringstream ss;

    char snapShotNumberString[3];
    int fileNumber;
    if(Continuation  == problem.IVP_Type)
        fileNumber = problem.NumFrames + timeStamp;
    else
        fileNumber = timeStamp;

    sprintf (snapShotNumberString, "%4.4d", fileNumber);


    // Use stringstream toconcatenate the strings using C++ including '\0'
    // character to signify end of string.
    ss << prefix << snapShotNumberString << extension << '\0';

    // convert the string stream into a C++ string and then int a C string.
    filename = ss.str ();
    char str[filename.size ()];

    // convert the C++ string to a C string by copying by value.
    for (unsigned int i = 0; i < filename.size (); ++i)
    {
        str[i] = filename[i];
    }

    // Now can finally use the string str = "uUvals(snapShotNumber).txt\0"
    ptr.open (str);

    ddpDenseVector_type SemiEFDof, ElectroEFDof, SemiPotDof, ElectroPotDof;
    poisson.getSemiconductorElecFieldDOFS(SemiEFDof);
    poisson.getElectrolyteElecFieldDOFS(ElectroEFDof);
    poisson.getSemiconductorPotentialDOFS(SemiPotDof);
    poisson.getElectrolytePotentialDOFS(ElectroPotDof);

    ddpDenseVector_type
    XPTS_semiconductor = electrons.PTSDense,
    XPTS_electrolyte = reductants.PTSDense,
    ElectronPTS = electrons.VandeMondeMatrices.globalVandeMondeDG
                  * electrons.carrierState.uDof,
                  HolePTS = holes.VandeMondeMatrices.globalVandeMondeDG
                            * holes.carrierState.uDof,
                            ReductantPTS =
                                reductants.VandeMondeMatrices.globalVandeMondeDG
                                * reductants.carrierState.uDof,
                                OxidantPTS = oxidants.VandeMondeMatrices.globalVandeMondeDG
                                        * oxidants.carrierState.uDof,
                                        ElectronsQPTS = 1.0
                                                * electrons.VandeMondeMatrices.globalVandeMondeDG
                                                * electrons.carrierState.qDof,
                                                HolesQPTS =  1.0
                                                        * holes.VandeMondeMatrices.globalVandeMondeDG
                                                        * holes.carrierState.qDof,
                                                        ReductantsQPTS = 1.0 //5e3
                                                                * reductants.VandeMondeMatrices.globalVandeMondeDG
                                                                * reductants.carrierState.qDof,
                                                                OxidantsQPTS = oxidants.VandeMondeMatrices.globalVandeMondeDG
                                                                        * oxidants.carrierState.qDof,
                                                                        SemiEFPTS = electrons.VandeMondeMatrices.globalVandeMondeMX
                                                                                * SemiEFDof,
                                                                                ElectroEFPTS =
                                                                                        reductants.VandeMondeMatrices.globalVandeMondeMX
                                                                                        * ElectroEFDof,
                                                                                        SemiPotPTS = electrons.VandeMondeMatrices.globalVandeMondeDG
                                                                                                * SemiPotDof,
                                                                                                ElectroPotPTS =
                                                                                                        reductants.VandeMondeMatrices.globalVandeMondeDG
                                                                                                        * ElectroPotDof;

    ddpDenseVector_type SemiGenPTS =
        electrons.VandeMondeMatrices.globalVandeMondeDG
        * electrons.carrierState.RHSFromGeneration;


    ddpDenseVector_type ElecGenPTS =
        reductants.VandeMondeMatrices.globalVandeMondeDG
        * reductants.carrierState.RHSFromGeneration;


    double x_scale = 1.0;// problem.characteristicLength;
    double t_scale = 1.0; //problem.characteristicTime;
    double den_scale = 1.0; //problem.characteristicDensity;
    double pot_scale = 1.0; //problem.thermalVoltage;
    double elec_scale = 1.0;//pot_scale/x_scale;
    double curr_scale = 1.0; //problem.electronCharge * den_scale * x_scale / t_scale;

    int N = XPTS_semiconductor.size();

    // Prints the carriers on the semiconductor side
    for(int i = 0; i < N; ++i)
    {
        ptr << x_scale * XPTS_semiconductor(i)  << " "
            << den_scale * ElectronPTS(i) << " "
            << den_scale * HolePTS(i) << " "
            << elec_scale * SemiEFPTS(i) << " "
            << pot_scale * SemiPotPTS(i) << " "
            << curr_scale * (ElectronsQPTS(i) - HolesQPTS(i)) << " "
            << (den_scale/t_scale) * SemiGenPTS(i) << " "
            << tCurrent * t_scale
            << endl;
    }
    // Prints the carriers on the electrolyte side
    for(int i = 0; i < N; i++)
    {
        ptr << x_scale * XPTS_electrolyte(i) << " "
            << den_scale * ReductantPTS(i) << " "
            << den_scale * OxidantPTS(i) << " "
            << elec_scale * ElectroEFPTS(i) << " "
            << pot_scale * ElectroPotPTS(i) << " "
            << curr_scale * ReductantsQPTS(i) << " " // - 0.5*OxidantsQPTS(i)) << " "
            << (den_scale/t_scale) * ElecGenPTS(i) << " "
            << tCurrent * t_scale
            << endl;
    }

    ptr.close();

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
// PRINTS OUT THE PROGRESS OF THE TIME STEPPING
///////////////////////////////////////////////////////////////////////////////
int
ddpPrintFinalDofs(ddpCarrier_type const & electrons,
                  ddpCarrier_type const & holes,
                  ddpPoisson_type const & poisson)
{
    // get size of DOFS
    int N = electrons.carrierState.uDof.size();
    ofstream ptr;


    // print out the electrons
    ptr.open("FINAL_ELECTRONS_DOF.dat");
    for(int i=0; i<N; i++)
    {
        ptr << electrons.carrierState.uDof(i)
            << " "
            << electrons.carrierState.qDof(i)
            << endl;
    }
    ptr.close();

    // print out the holes
    N = holes.carrierState.uDof.size();
    ptr.open("FINAL_HOLES_DOF.dat");
    for(int i=0; i<N; i++)
    {
        ptr << holes.carrierState.uDof(i)
            << " "
            << holes.carrierState.qDof(i)
            << endl;
    }
    ptr.close();


    int PotDof = poisson.PoissonState.potDof.size();
    N = poisson.PoissonState.elecDof.size();
    ptr.open("FINAL_POISSON_DOF.dat");
    for(int i=0; i<N; i++)
    {
        if(i < PotDof)
            ptr << poisson.PoissonState.potDof(i);
        else
            ptr << 0.0;

        ptr << " ";
        ptr << poisson.PoissonState.elecDof(i);
        ptr << endl;
    }

    ptr.close();

    return 0;
}

int
ddpPrintFinalDofs(ddpCarrier_type const & electrons,
                  ddpCarrier_type const & holes,
                  ddpCarrier_type const & reductants,
                  ddpCarrier_type const & oxidants,
                  ddpPoisson_type const & poisson)
{
    // get size of DOFS
    int N = electrons.carrierState.uDof.size();
    ofstream ptr;


    // print out the electrons
    ptr.open("FINAL_ELECTRONS_DOF.dat");
    for(int i=0; i<N; i++)
    {
        ptr << electrons.carrierState.uDof(i)
            << " "
            << electrons.carrierState.qDof(i)
            << endl;
    }
    ptr.close();

    // print out the holes
    N = holes.carrierState.uDof.size();
    ptr.open("FINAL_HOLES_DOF.dat");
    for(int i=0; i<N; i++)
    {
        ptr << holes.carrierState.uDof(i)
            << " "
            << holes.carrierState.qDof(i)
            << endl;
    }
    ptr.close();

    N = reductants.carrierState.uDof.size();
    // print out the reductants
    ptr.open("FINAL_REDUCTANTS_DOF.dat");
    for(int i=0; i<N; i++)
    {
        ptr << reductants.carrierState.uDof(i)
            << " "
            << reductants.carrierState.qDof(i)
            << endl;
    }
    ptr.close();

    N = oxidants.carrierState.uDof.size();
    // print out the holes
    ptr.open("FINAL_OXIDANTS_DOF.dat");
    for(int i=0; i<N; i++)
    {
        ptr << oxidants.carrierState.uDof(i)
            << " "
            << oxidants.carrierState.qDof(i)
            << endl;
    }
    ptr.close();


    int PotDof = poisson.PoissonState.potDof.size();
    N = poisson.PoissonState.elecDof.size();
    ptr.open("FINAL_POISSON_DOF.dat");
    for(int i=0; i<N; i++)
    {
        if(i < PotDof)
            ptr << poisson.PoissonState.potDof(i);
        else
            ptr << 0.0;

        ptr << " ";
        ptr << poisson.PoissonState.elecDof(i);
        ptr << endl;
    }

    ptr.close();

    return 0;
}

int
ddpReadInFinalStates(ddpCarrier_type & electrons,
                     ddpCarrier_type & holes,
                     ddpPoisson_type & poisson)
{
    int N = electrons.carrierState.uDof.size();
    int counter = 0;
    double dof1;
    double dof2;
    std::ifstream read;

    // READ IN ELECTONS DOF
    read.open("FINAL_ELECTRONS_DOF.dat");
    if(read)
    {
        while(!read.eof() && counter < N )
        {
            read >> dof1;
            read >> dof2;
            electrons.carrierState.uDof(counter) = dof1;
            electrons.carrierState.qDof(counter) = dof2;
            counter++;
        }
    }
    else
    {
        cout << "Couldnt Open FINAL_ELECTRONS_DOF.dat " << endl;
    }
    read.close();


    // READ IN HOLES DOF
    N = holes.carrierState.uDof.size();
    counter = 0;
    read.open("FINAL_HOLES_DOF.dat");
    if(read)
    {
        while(!read.eof() && counter < N )
        {
            read >> dof1;
            read >> dof2;
            holes.carrierState.uDof(counter) = dof1;
            holes.carrierState.qDof(counter) = dof2;
            counter++;
        }
    }
    else
    {
        cout << "Couldnt Open FINAL_HOLES_DOF.dat " << endl;
    }
    read.close();

    // READ IN ELECTRIC FIELD AND POTENTIAL
    int PotDof = poisson.PoissonState.potDof.size();
    N = poisson.PoissonState.elecDof.size();
    counter = 0;
    read.open("FINAL_POISSON_DOF.dat");
    if(read)
    {
        while(!read.eof() && counter < N)
        {
            read >> dof1;
            read >> dof2;

            if(counter < PotDof)
                poisson.PoissonState.potDof(counter) = dof1;

            poisson.PoissonState.elecDof(counter) = dof2;
            counter++;
        }
    }
    else
    {
        cout << "Couldnt Open FINAL_POISSON_DOF.dat " << endl;
    }
    read.close();

    return 0;
}

int ddpReadInElectronDOFS(ddpDenseVector_type & electrons_u,
                          ddpDenseVector_type & electrons_q)
{
    int N = electrons_u.size();
    int counter = 0;
    double dof1;
    double dof2;
    std::ifstream read;

    // READ IN ELECTONS DOF
    read.open("FINAL_ELECTRONS_DOF.dat");
    if(read)
    {
        while(!read.eof() && counter < N )
        {
            read >> dof1;
            read >> dof2;
            electrons_u(counter) = dof1;
            electrons_q(counter) = dof2;
            counter++;
        }
    }
    else
    {
        cout << "Couldnt Open FINAL_ELECTRONS_DOF.dat " << endl;
    }
    read.close();
}


int ddpReadInFinalStates(ddpCarrier_type & electrons,
                         ddpCarrier_type & holes,
                         ddpCarrier_type & reductants,
                         ddpCarrier_type & oxidants,
                         ddpPoisson_type & poisson)
{
    int N = electrons.carrierState.uDof.size();
    int counter = 0;
    double dof1;
    double dof2;
    std::ifstream read;

    // READ IN ELECTONS DOF
    read.open("FINAL_ELECTRONS_DOF.dat");
    if(read)
    {
        while(!read.eof() && counter < N )
        {
            read >> dof1;
            read >> dof2;
            electrons.carrierState.uDof(counter) = dof1;
            electrons.carrierState.qDof(counter) = dof2;
            counter++;
        }
    }
    else
    {
        cout << "Couldnt Open FINAL_ELECTRONS_DOF.dat " << endl;
    }
    read.close();


    // READ IN HOLES DOF
    N = holes.carrierState.uDof.size();
    counter = 0;
    read.open("FINAL_HOLES_DOF.dat");
    if(read)
    {
        while(!read.eof() && counter < N )
        {
            read >> dof1;
            read >> dof2;
            holes.carrierState.uDof(counter) = dof1;
            holes.carrierState.qDof(counter) = dof2;
            counter++;
        }
    }
    else
    {
        cout << "Couldnt Open FINAL_HOLES_DOF.dat " << endl;
    }
    read.close();

    // READ IN REDUCTANT DOF
    N = reductants.carrierState.uDof.size();
    counter = 0;
    read.open("FINAL_REDUCTANTS_DOF.dat");
    if(read)
    {
        while(!read.eof() && counter < N )
        {
            read >> dof1;
            read >> dof2;
            reductants.carrierState.uDof(counter) = dof1;
            reductants.carrierState.qDof(counter) = dof2;
            counter++;
        }
    }
    else
    {
        cout << "Couldnt Open FINAL_REDUCTANTS_DOF.dat " << endl;
    }
    read.close();

    // READ IN OXIDANT DOF
    N = oxidants.carrierState.uDof.size();
    counter = 0;
    read.open("FINAL_OXIDANTS_DOF.dat");
    if(read)
    {
        while(!read.eof() && counter < N )
        {
            read >> dof1;
            read >> dof2;
            oxidants.carrierState.uDof(counter) = dof1;
            oxidants.carrierState.qDof(counter) = dof2;
            counter++;
        }
    }
    else
    {
        cout << "Couldnt Open FINAL_OXIDANTS_DOF.dat " << endl;
    }
    read.close();

    // READ IN ELECTRIC FIELD AND POTENTIAL
    int PotDof = poisson.PoissonState.potDof.size();
    N = poisson.PoissonState.elecDof.size();
    counter = 0;
    read.open("FINAL_POISSON_DOF.dat");
    if(read)
    {
        while(!read.eof() && counter < N)
        {
            read >> dof1;
            read >> dof2;

            if(counter < PotDof)
                poisson.PoissonState.potDof(counter) = dof1;

            poisson.PoissonState.elecDof(counter) = dof2;
            counter++;
        }
    }
    else
    {
        cout << "Couldnt Open FINAL_POISSON_DOF.dat " << endl;
    }
    read.close();

    return 0;
}

int
progressBar(const int & current, const int & total,
            const int & width = 50)
{
    if(0 == current)
        cout << "\nPercentage of job complete:" << endl;

//	if( (current != total) && (current % (total/100) != 0) )
    //	return 0;

    double ratio = current/(double)total;
    int c = ratio * width;

    cout << std::setw(3) << (int) (ratio*100) << "% [";

    for(int i = 0; i < c; i++)
        cout << "=";

    for(int i = c; i < width; i++)
        cout << " ";

    cout << "]";

    printf("\r");
    fflush(stdout);

    return 0;
}
