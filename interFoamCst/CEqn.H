/*
 *  NMTFD-2 Deliverable Task II: Continuous Species Transfer
 */

const dimensionedScalar deltaN1("deltaN1",1e-8/pow(average(alpha1.mesh().V()), 1.0/3.0));

volScalarField D_harmonic =  (D1*D2)/(alpha1*D2 + (1-alpha1)*D1)  ;

volScalarField D_arithmetic  = alpha1*D1 + (1-alpha1)*D2 ;

dimensionedScalar D1_D2 = D1-D2;



const surfaceVectorField& Sf = mesh.Sf();
const surfaceScalarField& magSf = mesh.magSf();
const surfaceVectorField N_hat_s(Sf/magSf);
const volVectorField gradAlpha1(fvc::grad(alpha1));
const surfaceVectorField gradAlphaf1(fvc::interpolate(gradAlpha1));
const surfaceVectorField n_hat_f(gradAlphaf1/(mag(gradAlphaf1) + deltaN1));
const surfaceScalarField gamma_f = N_hat_s & n_hat_f;		
const surfaceScalarField gamma_f_magnitude = mag(gamma_f);		

// Coefficients to simplify CST Equation

volScalarField beta = 1/(alpha1 + ((1-alpha1)/H)); 
volScalarField K4 = D_harmonic * beta * (1-1/H);
volScalarField K5 = D1_D2 * alpha1 * (beta-1);
volScalarField K6 = beta * ((1/H*beta*D1_D2) - (D1-D2/H));

// Interpolating coefficients for fvm::div() and fvm::laplacian

surfaceScalarField psi4 = gamma_f_magnitude * fvc::interpolate(K4) * fvc::snGrad(alpha1) * mesh.magSf();
surfaceScalarField psi5 = (1-gamma_f_magnitude) * fvc::interpolate(K5);
surfaceScalarField psi6 = (1-gamma_f_magnitude) * fvc::interpolate(K6) * fvc::snGrad(alpha1) * mesh.magSf();

// CST Equation

fvScalarMatrix CEqn
(
      fvm::ddt(C)
    + fvm::div(phi,C)
    - fvm::laplacian(fvc::interpolate(D_arithmetic), C)
    + fvm::div(psi4, C)
    - fvm::laplacian(psi5, C)
    - fvm::div(psi6, C)
);

Info<< "\nSolving for C" << endl;
CEqn.solve();

Info<< "Min C: " << min(C).value() << endl;
Info<< "Max C: " << max(C).value() << endl;
