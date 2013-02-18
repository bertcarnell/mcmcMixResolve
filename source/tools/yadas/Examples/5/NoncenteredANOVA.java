import gov.lanl.yadas.*;
import java.util.*;

public class NoncenteredANOVA {

    public static void main (String[] args) {

	int B = 5000;
	String direc = "";
	String filename = "data.dat";
	String filename2 = "alpha.dat";
	String shortfilename = "scalars.dat";
	
	try {
	    if (args.length > 0)
		B = Integer.parseInt(args[0]);
	    if (args.length > 1) 
		direc = args[1];
	    if (args.length > 2)
		filename = args[2];
	    if (args.length > 3) 
		filename2 = args[3];
	    if (args.length > 4) 
		shortfilename = args[4];
	}
	catch (NumberFormatException e) {
	    System.out.println("Poor argument list!" + e);
	}
	
	DataFrame d = new DataFrame (direc + filename);
	DataFrame d2 = new DataFrame (direc + filename2);
	DataFrame d0 = new DataFrame (direc + shortfilename);

	// define MCMCParameters.  

	MCMCParameter mu, alpha;
	
	MCMCParameter[] paramarray = new MCMCParameter[] 
	{ 
	    mu = new MCMCParameter (d0.r("mu"), d0.r("mumss"), direc+"muplus"),
	    alpha = new MCMCParameter (d2.r("alpha"), d2.r("alphamss"), 
				       direc + "alphaplus"),
	};

	MCMCBond databond, alphaprior;

	ArrayList bondlist = new ArrayList ();

	bondlist.add ( databond = new BasicMCMCBond 
		       ( new MCMCParameter[] { mu, alpha },
			 new ArgumentMaker[] 
			   { new ConstantArgument (d.r("y")),
			     new FunctionalArgument 
				 (d.length(), 2, new int[] {0, 1}, 
				  new int[][] { d.i(0), d.i("group") },
				  new Function () 
				      { public double f (double[] args) {
					  return args[0] + args[1]; }}),
			     new ConstantArgument (d.r(d0.r("sigma")[0])) },
			 new Gaussian () ));

	bondlist.add ( alphaprior = new BasicMCMCBond 
		       ( new MCMCParameter[] { alpha },
			 new ArgumentMaker[] 
			   { new IdentityArgument (0),
			     new ConstantArgument (d2.r(d0.r("aalpha")[0])),
			     new ConstantArgument (d2.r(d0.r("balpha")[0])) },
			 new Gaussian () ));

	/* 
	   define MCMCUpdates:
	*/

	MCMCUpdate[] updatearray = new MCMCUpdate[] 
	{ 
	    mu, alpha, 
	};

	for (int b = 0; b < B; b++) {
	    if ((b/1000.0 - (int)(b/1000))== 0) System.out.println(b);
	    for (int i = 0; i < updatearray.length; i++) {
		updatearray[i].update ();
	    }
	    for (int i = 0; i < paramarray.length; i++) {
		paramarray[i].output ();
	    }
	}

	updatearray = new MCMCUpdate[] 
	{ 
	    mu, alpha, 
	    new MultipleParameterUpdate 
		( new MCMCParameter[] { mu, alpha },
		  new NewOneUpOneDownPerturber 
		      ( new int[][] { d0.i(0), d2.i(0) }, d0.r("mamss")) ),
	};

	for (int b = 0; b < B; b++) {
	    if ((b/1000.0 - (int)(b/1000))== 0) System.out.println(b);
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

	for (int i = 0; i < paramarray.length; i++) {
	    paramarray[i].finish();
	}
    }
}

