package clade;

import ml.CellNet;
import dna.Data;
import shared.Tools;
import structures.FloatList;
import tax.TaxTree;

/**
 * Estimates the probability that a QuickClade hit is taxonomically correct
 * at a given level, using neural networks with calibration via
 * K*sigmoid(a*logit(x)+b)^c.
 *
 * Hybrid model: all-length NNs for order-domain (5 models),
 * per-length-bin NNs for species-family (33 models, 3 levels × 11 bins).
 * Falls back to 5-parameter sigmoid model if .bbnets file is missing.
 *
 * @author Brian Bushnell, Ady
 */
public class CladeConfidence {

	public static float probCorrect(int length, float gcdif, float strdif,
			float hhdif, float cagadif,
			float k3dif, float k4dif, float k5dif, int taxLevel) {
		int idx = levelToIndex(taxLevel);
		if(idx<0){return -1;}
		if(USE_NN && loaded!=null && k5dif<1.0f){
			CellNet net;
			float[] cal;
			if(loaded.allLenNets[idx]!=null){
				net=loaded.allLenNets[idx];
				cal=loaded.allLenCal[idx];
			}else{
				int bin=binIndex(length);
				net=loaded.binNets[idx][bin];
				cal=loaded.binCal[idx][bin];
			}
			if(net!=null){
				return predictNN(net, cal, idx, length, gcdif, strdif, hhdif, cagadif, k3dif, k4dif, k5dif);
			}
		}
		return predictSigmoid(idx, length, k3dif, k4dif, k5dif);
	}

	/** Backwards-compatible call without gcdif/strdif/hhdif/cagadif; always uses sigmoid. */
	public static float probCorrect(int length, float k3dif, float k4dif, float k5dif, int taxLevel) {
		int idx = levelToIndex(taxLevel);
		if(idx<0){return -1;}
		return predictSigmoid(idx, length, k3dif, k4dif, k5dif);
	}

	/*---------- NN prediction ----------*/

	private static float predictNN(CellNet master, float[] cal, int idx, int length,
			float gcdif, float strdif, float hhdif, float cagadif,
			float k3dif, float k4dif, float k5dif) {
		CellNet net = getThreadNet(master, idx, length);
		FloatList fl = getThreadInput();
		fl.size=0;
		fl.add((float)(Math.log(Math.max(length, 1)) * INV_LN2));
		fl.add((float)(0.001 * Math.sqrt(Math.max(length, 1))));
		fl.add(gcdif);
		fl.add(strdif);
		fl.add(hhdif);
		fl.add(cagadif);
		fl.add(k3dif);
		fl.add(k4dif);
		fl.add(k5dif);
		fl.add((float)(k3dif / (k5dif + 0.01)));
		net.applyInput(fl);
		float raw = net.feedForward();
		return calibrate(raw, cal);
	}

	/** p = K * sigmoid(a * logit(x) + b) ^ c */
	private static float calibrate(float x, float[] p) {
		x=Tools.mid(0.0001f, x, 0.9999f);
		float K=p[0], a=p[1], b=p[2], c=p[3];
		double lx=Math.log(x / (1.0 - x));
		double s=1.0 / (1.0 + Math.exp(-(a * lx + b)));
		return (float)Math.min(1.0, K * Math.pow(s, c));
	}

	/*---------- thread-local nets ----------*/

	private static CellNet getThreadNet(CellNet master, int idx, int length) {
		CellNet[][] local = threadNets.get();
		if(local==null){
			local = new CellNet[NUM_LEVELS][NUM_BINS+1];
			threadNets.set(local);
		}
		int slot;
		if(loaded.allLenNets[idx]!=null){
			slot=NUM_BINS;
		}else{
			slot=binIndex(length);
		}
		if(local[idx][slot]==null){
			local[idx][slot] = master.copy(false);
		}
		return local[idx][slot];
	}

	private static FloatList getThreadInput() {
		FloatList fl = threadInput.get();
		if(fl==null){
			fl = new FloatList(NUM_FEATURES);
			threadInput.set(fl);
		}
		return fl;
	}

	private static final ThreadLocal<CellNet[][]> threadNets = new ThreadLocal<>();
	private static final ThreadLocal<FloatList> threadInput = new ThreadLocal<>();

	/*---------- bin lookup ----------*/

	static int binIndex(int length) {
		int raw = 63 - Long.numberOfLeadingZeros(Math.max(1, (long)length / 2500));
		return Tools.mid(0, raw, NUM_BINS-1);
	}

	/*---------- sigmoid fallback ----------*/

	private static float predictSigmoid(int idx, int length, float k3dif, float k4dif, float k5dif) {
		double[] params = (k5dif < 1.0f) ? PROB_K5[idx] : PROB_K4[idx];
		float kdif = (k5dif < 1.0f) ? k5dif : k4dif;
		double floor = params[0];
		double log2len = Math.log(Math.max(length, 1)) * INV_LN2 - LOG2_REF;
		double z = params[1] + params[2] * log2len + params[3] * k3dif + params[4] * kdif;
		double sigmoid = 1.0 / (1.0 + Math.exp(z));
		return (float)(floor + (1.0 - floor) * sigmoid);
	}

	/*---------- level mapping ----------*/

	private static int levelToIndex(int taxLevel) {
		switch (taxLevel) {
			case TaxTree.SPECIES: return 0;
			case TaxTree.GENUS:   return 1;
			case TaxTree.FAMILY:  return 2;
			case TaxTree.ORDER:   return 3;
			case TaxTree.CLASS:   return 4;
			case TaxTree.PHYLUM:  return 5;
			case TaxTree.KINGDOM: return 6;
			case TaxTree.DOMAIN:  return 7;
			default: return -1;
		}
	}

	/*---------- initialization ----------*/

	private static SerialNNLoader.LoadedNets loadNets() {
		String path = Data.findPath("?confidence.bbnets", true);
		if(path==null){return null;}
		return SerialNNLoader.load(path);
	}

	/*---------- constants ----------*/

	static final boolean USE_NN = true;
	private static final int NUM_LEVELS = 8;
	private static final int NUM_BINS = 11;
	private static final int NUM_FEATURES = 10;
	private static final double INV_LN2 = 1.0 / Math.log(2);
	private static final double LOG2_REF = Math.log(2500) * INV_LN2;

	private static final SerialNNLoader.LoadedNets loaded = loadNets();

	// Sigmoid fallback parameters
	static final double[][] PROB_K4 = {
		{0.0113446, -3.5630508, 0.6403488, 51.9226674, 13.6112786},
		{0.0795275, -5.0429731, 0.4555063, 59.6982223, 10.8967730},
		{0.2127537, -4.4862228, 0.2292117, 48.4028932,  7.7315173},
		{0.3256019, -5.2138624, 0.1925119, 49.9547608,  7.6156285},
		{0.5084922, -5.4920396, 0.0667606, 63.4318036,  0.5094376},
		{0.6071749, -6.5750378, 0.0950648, 57.9167172,  6.3128603},
		{0.6801771, -7.1822199, 0.1554127, 43.4878516, 15.1889807},
		{0.8182326, -8.3517371, 0.1448564, 32.9901037, 21.5209495},
	};

	static final double[][] PROB_K5 = {
		{0.0058409, -2.5413033, 0.4705238, 70.8753575, -0.7640568},
		{0.0556486, -3.6911274, 0.3030928, 71.3175058, -2.1744591},
		{0.1795990, -3.3621772, 0.1272263, 57.1529654, -2.5656541},
		{0.2788433, -4.2430087, 0.1167844, 57.7251810, -2.0374220},
		{0.4376369, -4.2088439,-0.0017810, 63.6425910, -5.5906412},
		{0.5352429, -5.0486971, 0.0025206, 64.2620784, -4.7840435},
		{0.6134413, -5.4231686, 0.0501707, 55.7766015, -2.2099322},
		{0.7916569, -7.1943760, 0.0995638, 49.9422788,  2.2346173},
	};
}
