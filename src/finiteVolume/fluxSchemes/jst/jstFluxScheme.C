/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 2021 Engys Ltd
    Copyright (C) 1991-2008 OpenCFD Ltd.

-------------------------------------------------------------------------------
License
    This file is part of HiSA.

    HiSA is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HiSA is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with HiSA.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "jstFluxScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "bound.H"
#include "fvcSurfaceReconstruct.H"
#include "orthogonalSnGrad.H"
#include "fvcMeshPhi.H"
#include "linear.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "thermodynamicConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(jstFluxScheme, 0);
addToRunTimeSelectionTable(fluxScheme, jstFluxScheme, dictionary);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

jstFluxScheme::jstFluxScheme
(
    const dictionary& dict,
    const fluidThermo& thermo,
    const volScalarField& rho,
    const volVectorField& U,
    const volVectorField& rhoU,
    const volScalarField& rhoE
)
:
    fluxScheme(typeName, dict, U.mesh()),
    mesh_(U.mesh()),
    thermo_(thermo),
    rho_(rho),
    U_(U),
    rhoU_(rhoU),
    rhoE_(rhoE),
    dict_(dict)
{
    //Info << mesh_.thisDb().classes() << endl;
    const IOdictionary& fvSchemes = mesh_.thisDb().lookupObject<IOdictionary>("fvSchemes");
    jst_k2_ = fvSchemes.getScalar("jst_k2");
    jst_k4_ = fvSchemes.getScalar("jst_k4");

    // get molWeight and Cp from thermophysicalProperties
    const IOdictionary& thermoDict = mesh_.thisDb().lookupObject<IOdictionary>("thermophysicalProperties");
    dictionary mixSubDict = thermoDict.subDict("mixture");
    dictionary specieSubDict = mixSubDict.subDict("specie");
    molWeight_ = specieSubDict.getScalar("molWeight");
    dictionary thermodynamicsSubDict = mixSubDict.subDict("thermodynamics");
    Cp_ = thermodynamicsSubDict.getScalar("Cp");
    
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

jstFluxScheme::~jstFluxScheme()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::jstFluxScheme::calcFlux(surfaceScalarField& phi, surfaceVectorField& phiUp, surfaceScalarField& phiEp, surfaceVectorField& Up)
{
    
    const volScalarField& p = thermo_.p();
    const fvMesh& mesh = mesh_;
    const volScalarField& T = mesh_.thisDb().lookupObjectRef<volScalarField>("T");

    // add the central term
    phi = linearInterpolate(rhoU_) & mesh.Sf();
    phiUp = linearInterpolate(rhoU_*U_ + p*tensor::I) & mesh.Sf();
    phiEp = linearInterpolate((rhoE_ + p)*U_) & mesh.Sf();

    // force to use oriented fluxes with signs
    phi.setOriented();
    phiUp.setOriented();
    phiEp.setOriented();

    // Wave speed
    volScalarField c(sqrt(thermo_.gamma()/thermo_.psi()));
    surfaceScalarField lambdaConv = (fvc::interpolate(c) + mag(fvc::interpolate(U_)&mesh_.Sf()/mesh_.magSf()));
    // lambdaConv should have no signs
    lambdaConv.setOriented(false);

    // JST pressure switch
    /*
    dimensionedScalar smallP
    (
        "smallP",                            
        dimensionSet(1, -1, -2, 0, 0, 0, 0), 
        1e-16                     
    );
    surfaceScalarField pDiff = mag(fv::orthogonalSnGrad<scalar>(mesh).snGrad(p)) / mesh.deltaCoeffs();
    surfaceScalarField pSum  = 2.0 * linearInterpolate(p);
    surfaceScalarField pSensor = pDiff / (pSum + smallP);
    pSensor.setOriented(false);
    pSensor = min(max(pSensor, scalar(0)), scalar(1));

    // JST artificial dissipation coefficients
    surfaceScalarField eps2 = jst_k2_ * pSensor * lambdaConv;
    eps2.setOriented(false);
    dimensionedScalar eps_zero
    (
        "jst_k4",                            
        lambdaConv.dimensions(), 
        0.0                   
    );
    surfaceScalarField eps4 = max(eps_zero, jst_k4_*lambdaConv-eps2);
    eps4.setOriented(false);
    */

   scalar RR = Foam::constant::thermodynamic::RR;

   scalar R = RR / molWeight_;
    
    dimensionedScalar TRef("TRef", dimTemperature, 300.0);
    dimensionedScalar pRef("pRef", dimPressure,   101325.0);

    volScalarField s = Cp_ * log(T/TRef) - R *log(p/pRef);
    dimensionedScalar smallS
    (
        "smallP",                            
        s.dimensions(), 
        1e-16                     
    );
    surfaceScalarField sDiff = mag(fv::orthogonalSnGrad<scalar>(mesh).snGrad(s)) / mesh.deltaCoeffs();
    surfaceScalarField sSum  = 2.0 * linearInterpolate(s);
    surfaceScalarField sSensor = sDiff / (sSum + smallS);
    sSensor.setOriented(false);
    sSensor = min(max(sSensor, scalar(0)), scalar(1));

    surfaceScalarField eps2 = jst_k2_ * sSensor * lambdaConv;
    eps2.setOriented(false);
    dimensionedScalar eps_zero
    (
        "jst_k4",                            
        lambdaConv.dimensions(), 
        0.0                   
    );
    surfaceScalarField eps4 = max(eps_zero, jst_k4_*lambdaConv-eps2);
    eps4.setOriented(false);

    // delta = snGrad * d. we have to force to use the ortho snGrad without any corrections
    surfaceScalarField deltaRhoPhi = fv::orthogonalSnGrad<scalar>(mesh).snGrad(rho_) / mesh.deltaCoeffs(); 
    surfaceVectorField deltaRhoUPhi = fv::orthogonalSnGrad<vector>(mesh).snGrad(rhoU_) / mesh.deltaCoeffs(); 
    surfaceScalarField deltaRhoEPhi = fv::orthogonalSnGrad<scalar>(mesh).snGrad(rhoE_) / mesh.deltaCoeffs(); 
    deltaRhoPhi.setOriented();
    deltaRhoUPhi.setOriented();
    deltaRhoEPhi.setOriented();

    /*
    volScalarField meshV
    (
        IOobject
        (
            "meshV",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimVolume, 0.0),
        "zeroGradient"
    );
    
    meshV.primitiveFieldRef() = mesh.V();
    meshV.correctBoundaryConditions();

    volScalarField LRef = pow(meshV, 1.0/3.0);
    // LRef / d
    surfaceScalarField stretch = pow(linearInterpolate(LRef) * mesh.deltaCoeffs(), -2.0); 
    stretch.setOriented(false);
    stretch = min(scalar(1.0), max(scalar(0.01), stretch));
    */
    
    /*
    // sum over all face to get delta-Rho for each cell
    // volScalarField D2Rho = fvc::div(deltaRhoPhi * mesh.magSf());
    // volVectorField D2RhoU = fvc::div(deltaRhoUPhi * mesh.magSf());
    // volScalarField D2RhoE = fvc::div(deltaRhoEPhi * mesh.magSf());
    */

    dimensionedScalar oneCoeff("oneCoeff", dimless, 1.0);
    volScalarField D2Rho = fvc::laplacian(oneCoeff, rho_, "laplacian(jst)");
    volVectorField D2RhoU = fvc::laplacian(oneCoeff, rhoU_, "laplacian(jst)");
    volScalarField D2RhoE = fvc::laplacian(oneCoeff, rhoE_, "laplacian(jst)");

    // apply the snGrad on the delta-Rho to get fourth order dissipation term
    surfaceScalarField deltaD2RhoPhi = fv::orthogonalSnGrad<scalar>(mesh).snGrad(D2Rho) / mesh.deltaCoeffs() / mesh.deltaCoeffs() / mesh.deltaCoeffs();    
    surfaceVectorField deltaD2RhoUPhi = fv::orthogonalSnGrad<vector>(mesh).snGrad(D2RhoU) / mesh.deltaCoeffs() / mesh.deltaCoeffs() / mesh.deltaCoeffs();    
    surfaceScalarField deltaD2RhoEPhi = fv::orthogonalSnGrad<scalar>(mesh).snGrad(D2RhoE) / mesh.deltaCoeffs() / mesh.deltaCoeffs() / mesh.deltaCoeffs();  
    deltaD2RhoPhi.setOriented();
    deltaD2RhoUPhi.setOriented();
    deltaD2RhoEPhi.setOriented();
    
    // add the artificial terms
    phi -= (eps2 * deltaRhoPhi - eps4 * deltaD2RhoPhi)*mesh.magSf();
    phiUp -= (eps2 * deltaRhoUPhi - eps4 * deltaD2RhoUPhi)*mesh.magSf();
    phiEp -= (eps2 * deltaRhoEPhi - eps4 * deltaD2RhoEPhi)*mesh.magSf();

    // Face velocity for sigmaDotU (turbulence term)
    Up = linearInterpolate(U_)*mesh_.magSf();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
