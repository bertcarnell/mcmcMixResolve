import gov.lanl.yadas.*;
import java.util.*;

/* 
   Example 7 in the YADAS web site/tutorial.  
   We have two binomial samples, and we want to test the hypothesis that 
   the two probabilities are equal.  We do this by putting a mixture 
   prior on the two probabilities $(p_1, p_2)$.  Their prior is a mixture 
   between two independent betas, and a common beta under which they are 
   both equal.  
*/

public class BinomialHypothesisTest {

    public static void main (String[] args) {

	int B = 10000;
	String direc = "";
	String filename = "p.dat";
	String shortfilename = "scalars.dat";
	
	try {
	    if (args.length > 0)
		B = Integer.parseInt(args[0]);
	    if (args.length > 1) 
		direc = args[1];
	    if (args.length > 2)
		filename = args[2];
	    if (args.length > 3)
		shortfilename = args[3];
	}
	catch (NumberFormatException e) {
	    System.out.println("Poor argument list!" + e);
	}
	
	DataFrame d = new DataFrame (direc + filename);
	DataFrame d0 = new DataFrame (direc + shortfilename);

	// define MCMCParameters.  
	// there is a single parameter p = (p_1, p_2)

	MCMCParameter p;
	
	MCMCParameter[] paramarray = new MCMCParameter[] 
	{ 
	    p = new MCMCParameter (d.r("p"), d.r(1.0), direc + "p"),
	};

	/*  
	    This model contains two bonds: 
	    the data x_i ~ Binomial (n_i, p_i)
	    and the mixture prior for (p_1, p_2).
	*/

	MCMCBond databond, mixtureprior;

	ArrayList bondlist = new ArrayList ();

	bondlist.add( databond = new BasicMCMCBond 
	    ( new MCMCParameter[] { p }, 
	      new ArgumentMaker[] {
		  new ConstantArgument (d.r("x")),
		  new ConstantArgument (d.r("n")),
		  new IdentityArgument (0) },
	      new Binomial () ));

	bondlist.add ( mixtureprior = new MixtureBond 
	    ( new MCMCParameter[][] { new MCMCParameter[] { p },
				      new MCMCParameter[] { p } },
	      new ArgumentMaker[][] 
		{ new ArgumentMaker[] { new IdentityArgument (0),
					new ConstantArgument (d.r("ap")),
					new ConstantArgument (d.r("bp")) },
		  new ArgumentMaker[] { new GroupArgument (0, new int[] {0}),
					new ConstantArgument (d0.r("ap")),
					new ConstantArgument (d0.r("bp")) } },
	      new Likelihood[] { new Beta(), new Beta() },
	      new AreTheyEqualArgument (0, d0.r("lambda")[0]) ) );

	databond.setName("databond");
	mixtureprior.setName("mixtureprior");

	/* 
	   define MCMCUpdates:
	*/

	MCMCUpdate[] updatearray = new MCMCUpdate[] 
	{ 
	    new ReversibleJumpUpdate 
		(new MCMCParameter[] { p }, 2, 0, 
		 new double[] { d0.r("a00")[0], d0.r("a01")[0], d0.r("a10")[0],
				d0.r("a11")[0] },
		 new JumpPerturber[] 
		    { new LogitPerturber (0, d.r("pmss")),
		      new EqualizingPerturber (0, d.r("n")),
		      new SplittingLogitPerturber (0, d.r("psplitmss")),
		      new CommonLogitPerturber (0, d0.r("pmss")[0])},
		 direc),
	};

	for (int b = 0; b < B; b++) {
	    if ((b/1000.0 - (int)(b/1000))== 0) System.out.println(b);
	    for (int i = 0; i < updatearray.length; i++) {
		updatearray[i].update ();
	    }
	    for (int i = 0; i < paramarray.length; i++) {
		if ((b/1.0 - (int)(b/1))== 0) paramarray[i].output ();
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

