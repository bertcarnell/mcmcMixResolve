/* 
   LinearModelExample.java:  Example 3 on the YADAS tutorial.  
   Data y_{i} ~ N ( \mu_{i}, \sigma^2 ),
   \mu_{i} = \beta_0 + \beta_1 x_{1i} + \beta_2 x_{2i} + \gamma_{g_1(i)} + 
   \delta_{g_2(i)},
   \beta, \gamma, \delta have normal priors,
   \sigma has a Gamma prior.  
*/

import java.util.*;
import gov.lanl.yadas.*;

public class LinearModelExample {

    public static void main (String[] args) {

	/* Comment 1.
	   Here we begin the initialization phase, where we define
	   the input files.  Data, covariates, and group labels are in 
	   "data.dat", quantities (initial values, step sizes, and prior 
	   parameters) relevant to \beta are in "beta.dat", similarly 
	   for \gamma and \delta, and quantities relevant to \sigma are in 
	   "scalars.dat".  
	*/

	int B = 10650;
	String direc = "";
	String filename = "data.dat";
	String filenameb = "beta.dat";
	String filenameg = "gamma.dat";
	String filenamed = "delta.dat";
	String shortfilename = "scalars.dat";
	
	try {
	    if (args.length > 0)
		B = Integer.parseInt(args[0]);
	    if (args.length > 1) 
		direc = args[1];
	    if (args.length > 2)
		filename = args[2];
	    if (args.length > 3) 
		filenameb = args[3];
	    if (args.length > 4) 
		filenameg = args[4];
	    if (args.length > 5) 
		filenamed = args[6];
	    if (args.length > 6) 
		shortfilename = args[6];
	}
	catch (NumberFormatException e) {
	    System.out.println("Poor argument list!" + e);
	}
	
	DataFrame d = new DataFrame (direc + filename);
	DataFrame dbeta = new DataFrame (direc + filenameb);
	DataFrame dgamma = new DataFrame (direc + filenameg);
	DataFrame ddelta = new DataFrame (direc + filenamed);
	ScalarFrame d0 = new ScalarFrame (direc + shortfilename);

	/* Comment 2.  
	   At this stage we define the parameters.  In this application, 
	   the parameters to be updated are beta (the linear regression 
	   coefficients), gamma and delta (the random effects), and sigma 
	   (the data standard deviation).
	   The definition of a parameter consists of an array of initial 
	   values, an array of step sizes, and a file name where the 
	   samples from the distribution of that parameter will be sent.  
	 */

	MCMCParameter beta, gamma, delta, sigma;
	
	MCMCParameter[] paramarray = new MCMCParameter[] 
	{ 
	    beta = new MCMCParameter (dbeta.r("beta"), dbeta.r("betamss"), 
				       direc + "beta"),
	    gamma = new MCMCParameter (dgamma.r("gamma"), dgamma.r("gammamss"), 
				       direc + "gamma"),
	    delta = new MCMCParameter (ddelta.r("delta"), ddelta.r("deltamss"), 
				       direc + "delta"),
	    sigma = new MCMCParameter (d0.r("sigma"), d0.r("sigmamss"), 
				       direc + "sigma"),
	};

	/* Comment 3.
	   At this stage we define the bonds, which combine to define the 
	   statistical model being analyzed.  The four bonds are: 
	   The first bond is the bond of interest in this example as it 
	   introduces LinearModelArgument.  
	   The arguments to the LinearModelArgument constructor are a 
	   data frame (d); an indicator of whether the model contains 
	   an intercept (1); an integer indicating which of the parameters 
	   is the vector of (intercept and) slope coefficients (the zeroth 
	   parameter, beta); an array of integers mapping the elements of 
	   beta to real-valued variables in the data frame (the 1 and 2 
	   mean that the #1 and #2 variables in the data frame, 'x1' and 'x2' 
	   since 'y' is the #0 variable, map to the first two elements of 
	   beta after the intercept; and an array of integers mapping the 
	   random effect parameters to the integer-valued variables in the 
	   data frame: the 1 and 2 mean that the #1 and #2 parameters, 
	   gamma and delta, go with the first two integer-value variables, 
	   'group1' and 'group2'.  One could add a sixth argument, a Function, 
	   if a link function were desired.  
	 */

	MCMCBond databond, betaprior, gammaprior, deltaprior, sigmaprior;

	ArrayList bondlist = new ArrayList ();

	bondlist.add ( databond = new BasicMCMCBond 
		       ( new MCMCParameter[] { beta, gamma, delta, sigma },
			 new ArgumentMaker[] 
			   { new ConstantArgument (d.r("y")),
			     new LinearModelArgument (d, 1, 0, new int[] { 1, 2 },
						      new int[] { 1, 2 }),
			     new GroupArgument (3, d.i(0)),
			   },
			 new Gaussian () ));

	bondlist.add ( betaprior = new BasicMCMCBond 
		       ( new MCMCParameter[] { beta },
			 new ArgumentMaker[] 
			   { new IdentityArgument (0),
			     new ConstantArgument (dbeta.r("abeta")),
			     new ConstantArgument (dbeta.r("bbeta")) },
			 new Gaussian () ));

	bondlist.add ( gammaprior = new BasicMCMCBond 
		       ( new MCMCParameter[] { gamma },
			 new ArgumentMaker[] 
			   { new IdentityArgument (0),
			     new ConstantArgument (dgamma.r("agamma")),
			     new ConstantArgument (dgamma.r("bgamma")) },
			 new Gaussian () ));

	bondlist.add ( deltaprior = new BasicMCMCBond 
		       ( new MCMCParameter[] { delta },
			 new ArgumentMaker[] 
			   { new IdentityArgument (0),
			     new ConstantArgument (ddelta.r("adelta")),
			     new ConstantArgument (ddelta.r("bdelta")) },
			 new Gaussian () ));

	bondlist.add ( sigmaprior = new BasicMCMCBond 
		       ( new MCMCParameter[] { sigma },
			 new ArgumentMaker[] 
			   { new IdentityArgument (0),
			     new ConstantArgument (d0.r("asigma")),
			     new ConstantArgument (d0.r("bsigma")) },
			 new Gamma () ));

	/* Comment 4. 
	   Now we define the MCMC algorithm by listing the update steps 
	   contained in it.  Listing individual parameters means that 
	   they will be updated by componentwise Metropolis steps, 
	   with step sizes as listed in the MCMCParameter definitions.  
	   These update steps are not adequate by themselves to ensure 
	   good mixing of the MCMC.  Due to correlation between beta[0] 
	   and delta[i] + sigma[j], MultipleParameterUpdates are needed, 
	   and we are not ready to introduce them yet.  
	 */

	MCMCUpdate[] updatearray = new MCMCUpdate[] 
	{ 
	    new UpdateTuner (beta), 
	    new UpdateTuner (gamma), 
	    new UpdateTuner (delta), 
	    new UpdateTuner (sigma),
	};

	/* Comment 5.  
	   The remainder of the code is the same for essentially all 
	   YADAS applications.  The first for() loop says: for each 
	   MCMC iteration, loop over the updates and update them, 
	   then loop over the parameters and output their current 
	   values.  The final two loops list acceptance rates for 
	   each of the update steps, then close all the output files.  
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

