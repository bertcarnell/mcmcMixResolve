import gov.lanl.yadas.*;
import java.util.*;

// adapted from NewAddCommonPerturber.  General solution would be good.  

public class NewOneUpOneDownPerturber implements Perturber {

    public NewOneUpOneDownPerturber (int[][] expandmat, double[] mss) {
	this.expandmat = expandmat;
	this.mss = mss;
    }

    public void perturb (double[][] candarray, int whoseTurn) {
	Integer wt = new Integer(whoseTurn);
	double temp = mss[whoseTurn] * rand.nextGaussian();
	int where = 0;
	int j = 0;
	for (int k = 0; k < candarray[j].length; k++) {
	    if (expandmat[j][k] == whoseTurn) {
		candarray[j][k] += temp;
	    }
	}
	j = 1;
	for (int k = 0; k < candarray[j].length; k++) {
	    if (expandmat[j][k] == whoseTurn) {
		candarray[j][k] -= temp;
	    }
	}
    }

    public int numTurns () {
	return mss.length;
    }

    public double jacobian () {
	return 1.0;
    }

    private int[][] expandmat;
    private double[] mss;

    static Random rand = new Random(System.currentTimeMillis());
}
