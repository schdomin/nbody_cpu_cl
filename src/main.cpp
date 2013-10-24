#include "CCubicDomain.h" //ds domain structure
#include "Timer.h"        //ds time measurement
#include <iostream>       //ds cout



int main( int argc, char** argv )
{
    //ds start timing
    Timer tmTimer; tmTimer.start( );

    //ds domain configuration
    const std::pair< double, double > pairBoundaries( -1.0, 1.0 );
    const unsigned int uNumberOfParticles( 100 );

    //ds allocate a domain to work with specifying number of particles and timing
    NBody::CCubicDomain cDomain( pairBoundaries, uNumberOfParticles );

    //ds target kinetic energy
    const double dTargetKineticEnergy( 1000.0 );

    //ds create particles uniformly from a normal distribution
    cDomain.createParticlesUniformFromNormalDistribution( dTargetKineticEnergy );

    //ds current simulation configuration
    const double dTimeStepSize( 0.0001 );
    const unsigned int uNumberOfTimeSteps( 5000 );
    const double dMinimumDistance( 5/uNumberOfParticles ); //pow( 1.0/uNumberOfParticles, 1.0/3 ) ); <- causes massive accelerations
    const double dPotentialDepth( 0.01 );

    std::cout << "--------CPU SETUP------------------------------------------------------------" << std::endl;
    std::cout << "  Number of particles: " << uNumberOfParticles << std::endl;
    std::cout << "Target kinetic energy: " << dTargetKineticEnergy << std::endl;
    std::cout << "  Number of timesteps: " << uNumberOfTimeSteps << std::endl;

    //ds information
    std::cout << "               Status:  0% done - current step: 0";

    //ds start simulation
    for( unsigned int uCurrentTimeStep = 0; uCurrentTimeStep < uNumberOfTimeSteps; ++uCurrentTimeStep )
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
