
//Volume of the bubble : sum of volume of each CV containing gas i.e (1-alpha1)
//You can access the volume of a CV as :  mesh.V() Volume of one CV

//containing gas : (1-alpha1)*mesh.V()

//To calculate mass transfer coefficient  use the following formula.
   
   // kappa = (bubble_volume * Delta_C)/(Bubble_area * Delta_t)
//	Delta_C = C at (n+1) timestep – C at n timestep.

//magnitude of gradient of alpha * the volume of the control

double bubleVol = 0; 
double bubleVol2 = sum((1-alpha1)*mesh.V());
double surfaceA = 0;
double Ctotal =0;//sum((1-alpha1)*C);
const volVectorField gradAlphaaa(fvc::grad(alpha1));
forAll(mesh.cells(),cellI)
	{ 
	Ctotal +=  ((1-alpha1[cellI])*mesh.V()[cellI])*C[cellI];
        surfaceA += mag(gradAlphaaa[cellI])* mesh.V()[cellI];
        bubleVol=bubleVol+(1-alpha1[cellI])*mesh.V()[cellI];//
	}
Info << "Caverage: " << Ctotal/bubleVol<<endl;
Info << "volume of the domain: " << bubleVol << endl;
Info <<"bubleVol2: "<<bubleVol2<<endl;
Info <<"surfaceA: "<<surfaceA<<endl;
Info << endl;
