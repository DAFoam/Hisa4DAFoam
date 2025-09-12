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
    dict_(dict),
    TRef_("TRef", dimTemperature, 300.0),
    pRef_("pRef", dimPressure, 101325.0)
{
    // read in parameters
    const IOdictionary& fvSchemes = mesh_.thisDb().lookupObject<IOdictionary>("fvSchemes");
    jst_k2_ = fvSchemes.getScalar("jst_k2");
    jst_k4_ = fvSchemes.getScalar("jst_k4");
    sensorName_ = fvSchemes.lookupOrDefault<word>("sensorName", "pressure");

    // get molWeight and Cp from thermophysicalProperties
    const IOdictionary& thermoDict = mesh_.thisDb().lookupObject<IOdictionary>("thermophysicalProperties");
    dictionary mixSubDict = thermoDict.subDict("mixture");
    dictionary specieSubDict = mixSubDict.subDict("specie");
    scalar molWeight = specieSubDict.getScalar("molWeight");
    dictionary thermodynamicsSubDict = mixSubDict.subDict("thermodynamics");
    Cp_ = thermodynamicsSubDict.getScalar("Cp");
    scalar RR = Foam::constant::thermodynamic::RR;
    R_ = RR / molWeight;
    
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

    // local spectral radius at faces  c + |U| * n
    volScalarField c(sqrt(thermo_.gamma()/thermo_.psi()));
    surfaceScalarField specR = (fvc::interpolate(c) + mag(fvc::interpolate(U_)&mesh_.Sf()/mesh_.magSf()));
    // specR should have no signs
    specR.setOriented(false);

    // shock sensor field: it can be pressure or entropy
    volScalarField sField("sField", p);
    if (sensorName_ == "pressure")
    {
        // do nothing
    }
    else if (sensorName_ == "entropy")
    {
        sField = Cp_ * log(T/TRef_) - R_ *log(p/pRef_);
        sField.correctBoundaryConditions();
    }
    else
    {
        FatalErrorIn("") << "sensor " << sensorName_ << " not supported. Options are: pressure or entropy" << exit(FatalError);
    }
    
    // calculate the sensor = |s_N - s_P | / |s_N + s_P + eps| we use the first order difference
    dimensionedScalar smallS
    (
        "smallS",                            
        sField.dimensions(), 
        1e-16                     
    );
    // |s_N - s_P |
    surfaceScalarField sDiff = mag(fv::orthogonalSnGrad<scalar>(mesh).snGrad(sField)) / mesh.deltaCoeffs();
    // |s_N + s_P |
    surfaceScalarField sSum  = 2.0 * linearInterpolate(sField);
    surfaceScalarField sensor = sDiff / (sSum + smallS);
    sensor.setOriented(false);
    // bound it between 0 and 1
    sensor = min(max(sensor, scalar(0)), scalar(1));

    // JST artificial dissipation coefficients
    surfaceScalarField eps2 = jst_k2_ * sensor;
    eps2.setOriented(false);
    surfaceScalarField eps4 = max(scalar(0.0), jst_k4_-eps2);
    eps4.setOriented(false);
    
    // dPhi = snGrad * d. we have to force to use the ortho snGrad without any corrections
    surfaceScalarField dRho = fv::orthogonalSnGrad<scalar>(mesh).snGrad(rho_) / mesh.deltaCoeffs(); 
    surfaceVectorField dRhoU = fv::orthogonalSnGrad<vector>(mesh).snGrad(rhoU_) / mesh.deltaCoeffs(); 
    surfaceScalarField dRhoE = fv::orthogonalSnGrad<scalar>(mesh).snGrad(rhoE_) / mesh.deltaCoeffs(); 
    dRho.setOriented();
    dRhoU.setOriented();
    dRhoE.setOriented();
    
    // sum over all face to get d2Rho_dx for each cell
    // NOTE: the fvc::div operator sums over all the face (without multiplying the area) and then divide it by the volume
    // So we have to manually multiply the surface area
    volScalarField d2Rho_dx = fvc::div(dRho * mesh.magSf());
    volVectorField d2RhoU_dx = fvc::div(dRhoU * mesh.magSf());
    volScalarField d2RhoE_dx = fvc::div(dRhoE * mesh.magSf());

    // apply the snGrad on the d2Rho_dx to get 3rd order differences
    // Note we have to do two mesh.deltaCoeffs() here because the first one is for snGrad and the 2nd one
    // is for d2Rho_dx (we want difference d2Rho, not d2Rho/dx). This makes sure the unit for d3Rho is the same as rho
    surfaceScalarField d3Rho = fv::orthogonalSnGrad<scalar>(mesh).snGrad(d2Rho_dx) / mesh.deltaCoeffs() / mesh.deltaCoeffs();    
    surfaceVectorField d3RhoU = fv::orthogonalSnGrad<vector>(mesh).snGrad(d2RhoU_dx) / mesh.deltaCoeffs() / mesh.deltaCoeffs();    
    surfaceScalarField d3RhoE = fv::orthogonalSnGrad<scalar>(mesh).snGrad(d2RhoE_dx) / mesh.deltaCoeffs() / mesh.deltaCoeffs();  
    d3Rho.setOriented();
    d3RhoU.setOriented();
    d3RhoE.setOriented();
    
    // add the artificial fluxes:
    // Note when integrating the dRho fluxes over the control volume, we get d2Rho/dx2 (2nd order dissipation)
    // same applies to d3Rho, integrating it will get d4Rho/dx4 (fourth order dissipation)
    phi -= (eps2 * dRho - eps4 * d3Rho) * mesh.magSf() * specR;
    phiUp -= (eps2 * dRhoU - eps4 * d3RhoU) * mesh.magSf() * specR;
    phiEp -= (eps2 * dRhoE - eps4 * d3RhoE) * mesh.magSf() * specR;

    // Face velocity for sigmaDotU (turbulence term)
    Up = linearInterpolate(U_)*mesh_.magSf();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
