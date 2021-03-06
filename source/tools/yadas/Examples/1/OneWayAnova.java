/* 
   OneWayAnova.java:  Example 1 on the YADAS tutorial.  
   Data y_{ij} ~ N ( \mu_{i}, \sigma^2 ),
   \mu_{i} ~ N ( \theta, \delta^2 ),
   \theta ~ flat prior on (-\infty, \infty),
   \sigma ~ Gamma (a_\sigma, b_\sigma),
   \delta ~ Gamma (a_\delta, b_\delta).  
*/

import java.util.*;
import gov.lanl.yadas.*;

public class OneWayAnova {

    public static void main (String[] args) {

	/* Comment 1.
	   Here we begin the initialization phase, where we define
	   the input files.  The data are in the file 'data.dat', 
	   the information concerning the main effects are in the 
	   file 'mu.dat', and the scalar parameters in the file 
	   'scalars.dat'.  B is the number of MCMC iterations.  
	   All of these file names, the directory name, and the 
	   number of iterations, can be changed using command line 
	   arguments.  Just before Comment 2, the data in these input 
	   files are inserted into DataFrames or ScalarFrames 
	   for easy access later.  In particular, the command 
	   d.r("y") extracts the real-valued variable named 'y' from 
	   the DataFrame 'd' (which here contains the data in the file 
	   Ex1data.dat).  
	*/

	int B = 10650;
	String direc = "";
	String filename = "Ex1data.dat";
	String filename2 = "Ex1mu.dat";
	String shortfilename = "Ex1scalars.dat";
	
	// The code below processes command-line arguments ("args") if any.  
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
	ScalarFrame d0 = new ScalarFrame (direc + shortfilename);

	/* Comment 2.  
	   At this stage we define the parameters.  In this application, 
	   the parameters to be updated are mu (the group means), 
	   theta (the overall mean), sigma (the data standard deviation),
	   and delta (the main effect standard deviation).  
	   The definition of a parameter consists of an array of initial 
	   values, an array of step sizes, and a file name where the 
	   samples from the distribution of that parameter will be sent. 
	   For example, there are three mu_i's, and d2.r("mu") is 
	   an array of length three (and so is d2.r("mss")).  
	   This analysis contains a single theta, sigma, and delta, 
	   so d0.r("theta") etc. are arrays of length one.  
	 */

	MCMCParameter mu, theta, sigma, delta;
	
	MCMCParameter[] paramarray = new MCMCParameter[] 
	{ 
	    mu = new MCMCParameter (d2.r("mu"), d2.r("mumss"), 
				       direc + "mu"),
	    theta = new MCMCParameter (d0.r("theta"), d0.r("thetamss"), 
				       direc + "theta"),
	    sigma = new MCMCParameter (d0.r("sigma"), d0.r("sigmamss"), 
				       direc + "sigma"),
	    delta = new MCMCParameter (d0.r("delta"), d0.r("deltamss"), 
				       direc + "delta"),
	};

	/* Comment 3.
	   At this stage we define the bonds, which combine to define the 
	   statistical model being analyzed.  The four bonds are: 
	   databond: y_{ij} ~ N ( \mu_i, \sigma^2 ),
	   muprior: \mu_i ~ N ( \theta, \delta^2 ),
	   sigmaprior: sigma ~ Gamma (a_\sigma, b_\sigma),
	   deltaprior: delta ~ Gamma (a_\delta, b_\delta).  
	   No bond is necessary to state that \theta has a flat prior on 
	   (-\infty, \infty).  
	   Note that the Gaussian distribution is parameterized by its 
	   standard deviation, not by its variance (and certainly not 
	   by its precision).  The Gamma distribution is parameterized 
	   by its shape and scale parameters, so that its mean is the 
	   product of its two parameters.  
	 */

	MCMCBond databond, muprior, sigmaprior, deltaprior;

	ArrayList bondlist = new ArrayList ();

	bondlist.add ( databond = new BasicMCMCBond 
		       ( new MCMCParameter[] { mu, sigma },
			 new ArgumentMaker[] 
			   { new ConstantArgument (d.r("y")),
			     new GroupArgument (0, d.i("group")),
			     new GroupArgument (1, d.i(0)) },
			 new Gaussian () ));

	bondlist.add ( muprior = new BasicMCMCBond 
		       ( new MCMCParameter[] { mu, theta, delta },
			 new ArgumentMaker[] 
			   { new IdentityArgument (0),
			     new GroupArgument (1, d2.i(0)),
			     new GroupArgument (2, d2.i(0)) },
			 new Gaussian () ));

	bondlist.add ( sigmaprior = new BasicMCMCBond 
		       ( new MCMCParameter[] { sigma },
			 new ArgumentMaker[] 
			   { new IdentityArgument (0),
			     new ConstantArgument (d0.r("asigma")),
			     new ConstantArgument (d0.r("bsigma")) },
			 new Gamma () ));

	bondlist.add ( deltaprior = new BasicMCMCBond 
		       ( new MCMCParameter[] { delta },
			 new ArgumentMaker[] 
			   { new IdentityArgument (0),
			     new ConstantArgument (d0.r("adelta")),
			     new ConstantArgument (d0.r("bdelta")) },
			 new Gamma () ));

	/* Comment 4. 
	   Now we define the MCMC algorithm by listing the update steps 
	   contained in it.  Listing individual parameters means that 
	   they will be updated by componentwise Metropolis steps, 
	   with step sizes as listed in the MCMCParameter definitions.  
	 */

	MCMCUpdate[] updatearray = new MCMCUpdate[] 
	{ 
	    new UpdateTuner (mu), 
	    new UpdateTuner (theta), 
	    new UpdateTuner (sigma), 
	    new UpdateTuner (delta),
	};

	/* Comment 5.  
	   The remainder of the code is the same for essentially all 
	   YADAS applications.  The first for() loop says: for each 
	   MCMC iteration, loop over the updates and update them, 
	   then loop over the parameters and output their current 
	   values.  The final two loops list acceptance rates for 
	   each of the update steps, then close all the output files.  
	   The 1000 means that every 1000 iterations, the code prints 
	   out how many iterations have been completed.  This number 
	   could be changed in slower applications.  
	 */

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

