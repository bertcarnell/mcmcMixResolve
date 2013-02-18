/* 
   AddCommonPerturber.java:  Example 4 on the YADAS tutorial.  
   Data y_{ij} ~ N ( \mu_{i}, \sigma^2 ),
   \mu_{i} ~ N ( \theta, \delta^2 ),
   \theta ~ flat prior on (-\infty, \infty),
   \sigma ~ Gamma (a_\sigma, b_\sigma),
   \delta ~ Gamma (a_\delta, b_\delta).  
*/

import java.util.*;
import gov.lanl.yadas.*;

public class AddCommonPerturberExample {

    public static void main (String[] args) {

	/* Comment 1.
	   This example is very similar to Example 1, OneWayAnova.java.  
	   This time we have priors and data that indicate that the 
	   data variance is much larger than the variance of the random 
	   effects (\sigma >> \delta).  In such a scenario, one would 
	   ordinarily reparameterize so that y_{ij} ~ N(\theta + \alpha_i, \sigma^2), 
	   \theta ~ flat prior, and \alpha_i ~ N(0, \delta^2) and this 
	   would lead to much better mixing by either the Gibbs sampler 
	   or a related Metropolis implementation.  
	   (See G. O. Roberts and S. K. Sahu, JRSSB 1997, 59:291-397.)
	   In this example we intentionally use the wrong parameterization 
	   and solve the problem using a MultipleParameterUpdate and an 
	   AddCommonPerturber.  The first half of the MCMC iterations are 
	   performed without this additional update, whereas the second 
	   half use it and mix more effectively.  
	*/

	int B = 5000;
	String direc = "";
	String filename = "data.dat";
	String filename2 = "mu.dat";
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
	ScalarFrame d0 = new ScalarFrame (direc + shortfilename);

	/* Comment 2.  
	   At this stage we define the parameters.  In this application, 
	   the parameters to be updated are mu (the group means), 
	   theta (the overall mean), sigma (the data standard deviation),
	   and delta (the main effect standard deviation).  
	   The definition of a parameter consists of an array of initial 
	   values, an array of step sizes, and a file name where the 
	   samples from the distribution of that parameter will be sent.  
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
	    mu, theta, sigma, delta, 
	    new MultipleParameterUpdate 
		( new MCMCParameter[] {mu, theta},
		  new NewAddCommonPerturber ( new int[][] { d2.i(0), d0.i(0) }, 
					      d0.r("mtmss") )), 
	};

	/* Comment 5.  
	   For illustrative purposes, we have changed the code that usually appears 
	   at this stage:  for the first B iterations, we run the algorithm without 
	   the MultipleParameterUpdate, and then for the next B, we run the 
	   algorithm with all the updates.  
	 */
	
	for (int b = 0; b < B; b++) {
	    if ((b/1000.0 - (int)(b/1000))== 0) System.out.println(b);
	    for (int i = 0; i < updatearray.length - 1; i++) {
		updatearray[i].update ();
	    }
	    for (int i = 0; i < paramarray.length; i++) {
		paramarray[i].output ();
	    }
	}
	    
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

