#include "CCubicDomain.h" //ds domain structure
#include "Timer.h"        //ds time measurement
#include <iostream>       //ds cout



int main( int argc, char** argv )
{
    //ds check simple input arguments - CAUTION: the implementation expects real numbers, the simulation will be corrupted if invalid values are entered
    if( 4 != argc )
    {
        //ds inform
        std::cout << "usage: nbody_cpu_cl [Number of particles] [Number of time steps] [Target energy]" << std::endl;
        return 0;
    }

    //ds start timing
    Timer tmTimer; tmTimer.start( );

    //ds domain configuration
    const std::pair< double, double > pairBoundaries( -1.0, 1.0 );
    const double dDomainWidth( fabs( pairBoundaries.first ) + fabs( pairBoundaries.second ) );
    const unsigned int uNumberOfParticles( atoi( argv[1] ) );

    //ds current simulation configuration
    const double dTimeStepSize( 0.0001 );
    const unsigned int uNumberOfTimeSteps( atoi( argv[2] ) );
    const double dMinimumDistance( pow( 1.0/uNumberOfParticles, 1.0/3 ) );
    const double dPotentialDepth( 1.0 );

    //ds target kinetic energy
    const double dTargetKineticEnergy( atol( argv[3] ) );

    //ds cell list information
    const unsigned int uNumberOfCells1D( floor( dDomainWidth/( 2.5*dMinimumDistance ) ) );
    const unsigned int uMaximumCellIndex( uNumberOfCells1D + pow( uNumberOfCells1D, 2 ) + pow( uNumberOfCells1D, 3 ) + 1 );

    std::cout << "------- CPU SETUP -----------------------------------------------------------" << std::endl;
    std::cout << "  Number of particles: " << uNumberOfParticles << std::endl;
    std::cout << "        Boundary (3D): [" << pairBoundaries.first << ", " << pairBoundaries.second << "]" << std::endl;
    std::cout << "         Domain Width: " << dDomainWidth << std::endl;
    std::cout << "     Minimum distance: " << dMinimumDistance << std::endl;
    std::cout << "      Cutoff distance: " << 2.5*dMinimumDistance << std::endl;
    std::cout << "      Potential depth: " << dPotentialDepth << std::endl;
    std::cout << "Target kinetic energy: " << dTargetKineticEnergy << std::endl;
    std::cout << " Number of time steps: " << uNumberOfTimeSteps << std::endl;
    std::cout << "       Time step size: " << dTimeStepSize << std::endl;
    std::cout << "------- CELL LISTS ----------------------------------------------------------" << std::endl;
    std::cout << " Number of cells 1D M: " << uNumberOfCells1D << std::endl;
    std::cout << "   Maximum cell index: " << uMaximumCellIndex << std::endl;
    std::cout << "-----------------------------------------------------------------------------" << std::endl;

    //ds allocate a domain to work with specifying number of particles and timing
    NBody::CCubicDomain cDomain( pairBoundaries, uNumberOfParticles, dMinimumDistance, uNumberOfCells1D, uMaximumCellIndex );

    //ds create particles uniformly from a normal distribution
    cDomain.createParticlesUniformFromNormalDistribution( dTargetKineticEnergy );

    //ds information
    std::cout << "               Status:  0% done - current step: 0";

    //ds start simulation
    for( unsigned int uCurrentTimeStep = 1; uCurrentTimeStep < uNumberOfTimeSteps+1; ++uCurrentTimeStep )
    {
        //ds calculate percentage done
        const double dPercentageDone( 100.0*uCurrentTimeStep/uNumberOfTimeSteps );

        //ds get a formatted string -> 100% -> 3 digits
        char chBuffer[4];

        //ds fill the buffer
        std::snprintf( chBuffer, 4, "%3.0f", dPercentageDone );

        //ds print info
        std::cout << '\xd';
        std::cout << "               Status: " << chBuffer << "% done - current step: " << uCurrentTimeStep;

        //ds update particles
        cDomain.updateParticlesVelocityVerlet( dTimeStepSize, dMinimumDistance, dPotentialDepth );

        //ds record situation (we will write the stream to the file in one operation afterwards )
        cDomain.saveParticlesToStream( );
        cDomain.saveIntegralsToStream( dMinimumDistance, dPotentialDepth );
    }

    //ds save the streams to a file
    cDomain.writeParticlesToFile( "bin/simulation.txt", uNumberOfTimeSteps );
    cDomain.writeIntegralsToFile( "bin/integrals.txt", uNumberOfTimeSteps, dTimeStepSize );

    //ds stop timing
    const double dDurationSeconds( tmTimer.stop( ) );

    //ds cause an output ostream
    std::cout << std::endl;
    std::cout << "     Computation time: " << dDurationSeconds << std::endl;
    std::cout << "-----------------------------------------------------------------------------" << std::endl;

    return 0;
}
