import gov.lanl.yadas.*;

public class AttritionAnalysis {

    /* 
       AttritionAnalysis: fit the attrition model to auto racing data.
       input file must be matrix of 
       raceid|finishid|driverid|trackid
    */

    public static void main (String[] args) {
	
	int B = 10650;
	String direc = "";
	String filename = "attw02.dat";
	String filename2 = "mssw02.dat";
	int numtracks = 23;
	
	try {
	    if (args.length > 0)
		B = Integer.parseInt(args[0]);
	    if (args.length > 1) 
		numtracks = Integer.parseInt(args[1]);
	    if (args.length > 2) 
		direc = args[2];
	    if (args.length > 3)
		filename = args[3];
	    if (args.length > 4) 
		filename2 = args[4];
	}
	catch (NumberFormatException e) {
	    System.out.println("Poor argument list!" + e);
	}
	
	DataFrame d = new DataFrame (direc + filename);
	DataFrame mss = new DataFrame (direc + filename2);
	int numdrivers = mss.length();
	
	// define MCMCParameters.  
	/*
	  0: theta
	  1: phi
	  2: atheta
	  3: btheta
	*/

	
	MCMCParameter paramarray[] = new MCMCParameter[4];
	MCMCParameter theta, phi, atheta, btheta;
	
	paramarray[0] = theta = new MCMCParameter(d.r(0.0, numdrivers), 
						  mss.r("mss"),
						  direc + "theta");
	paramarray[1] = phi = new MCMCParameter(d.r(1.0, numtracks), 
						d.r(0.0, numtracks), 
						direc + "phi");
	paramarray[2] = atheta = new MCMCParameter(d.r(0.0, 1), 
						   d.r(0.0, 1), 
						   direc + "atheta");
	paramarray[3] = btheta = new MCMCParameter(d.r(1.0, 1), 
						   d.r(0.25, 1), 
						   direc + "btheta");
	
	/*
	  define MCMCBonds
	  0: (finish position) - theta - phi
	  1: theta - atheta - btheta
	  2: phi - aphi - bphi
	*/
	
	MCMCBond bondarray[] = new MCMCBond[3];
	bondarray[0] = new BasicMCMCBond 
	    ( new MCMCParameter[] { theta, phi },
	      new ArgumentMaker[] { new GroupArgument (0, d.i("driverid")), 
				    new GroupArgument (1, d.i("trackid")) },
	      new AttritionLikelihood( d.i("raceid"), d.i("finishid") ) );
	bondarray[1] = new BasicMCMCBond 
	    ( new MCMCParameter[] {theta, atheta, btheta}, 
	      new ArgumentMaker[]{new IdentityArgument(0),
				  new GroupArgument(1, d.i(0, theta.length())),
				  new GroupArgument(2, d.i(0, theta.length()))},
	      new Gaussian() );
	bondarray[2] = new BasicMCMCBond 
	    ( new MCMCParameter[] {btheta},
	      new ArgumentMaker[] { new IdentityArgument(0),
				    new ConstantArgument(1.0),
				    new ConstantArgument(4.0) },
	      new Gamma());

	TunableMultipleParameterUpdate tmpu;
	MCMCUpdate[] updatearray = new MCMCUpdate[] 
	{ 
	    new UpdateTuner (theta), 
	    new UpdateTuner (btheta),
	    new UpdateTuner 
	    ( tmpu = new TunableMultipleParameterUpdate 
	      ( new MCMCParameter[] { theta, btheta },
		new ScalePerturber (0.20),
		direc + "multupdate") ),
	};

	for (int b = 0; b < B; b++) {
	    if (((int)(b/100)) == b/1000.0)
		System.out.println(b);
	    for (int i = 0; i < updatearray.length; i++) {
		updatearray[i].update ();
	    }
	    for (int i = 0; i < paramarray.length; i++) {
		paramarray[i].output ();
	    }
	}
      
	String acc;
	for (int iii = 0; iii < updatearray.length; iii++) {
	    acc = updatearray[iii].accepted();
	    System.out.println("Update " + iii + ": " + acc);
	}
	tmpu.finish();
	    
	for (int i = 0; i < paramarray.length; i++) {
	    paramarray[i].finish();
	}
    }
}

