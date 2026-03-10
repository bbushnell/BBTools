package cardinality;

import java.io.File;

import dna.Data;
import fileIO.ByteFile;
import fileIO.FileFormat;
import parse.LineParser1;
import parse.Parse;
import shared.Tools;
import structures.FloatList;

/**
 * Correction factors for estimators.
 * Loads a TSV file with per-occupancy correction factors and provides fast
 * lookup by estimator type and filled-bucket count, with quadratic interpolation
 * when the DDL bucket count differs from the loaded matrix size.
 * <p>
 * CF_MATRIX[type][filled_buckets] = correction factor to multiply rawEstimate by.
 *
 * @author Brian Bushnell, Chloe
 * @date March 5, 2026
 */
public class CorrectionFactor{

	/*--------------------------------------------------------------*/
	/*----------------            Main              ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args) {
		String fname=args[0];
		int buckets=(args.length>1 ? Parse.parseIntKMG(args[1]) : 2048);
		initialize(fname, buckets);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static synchronized void initialize(String fname, int buckets) {
		correctionFile=fname;
		correctionBuckets=buckets=CardinalityTracker.powerOf2AtLeast(buckets);
		CF_MATRIX=loadFile(fname, buckets);
	}

	public static synchronized float[][] loadFile(String fname, int buckets){
		String path=(fname.indexOf('?')>=0 ? Data.findPath(fname) : fname);
		File f=new File(path);
		if(!f.canRead() || !f.isFile()){
			System.err.println("Can't read "+fname);
			return null;
		}
		FloatList[] lists=null;
		FileFormat ff=FileFormat.testInput(path, null, false);
		ByteFile bf=ByteFile.makeByteFile(ff, 1);
		LineParser1 lp=new LineParser1('\t');
		lastCardMatrix=null;
		lastCardKeys=null;
		boolean inCardSection=false;
		FloatList[] cardLists=null;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			lp.set(line);
			if(lp.startsWith('#')){
				final String s=new String(line).trim();
				if(s.equals("#SECTION=cardinality")){inCardSection=true;}
				else if(s.equals("#SECTION=occupancy")){inCardSection=false;}
				continue;
			}
			if(!Tools.isNumeric(line[0])){continue;}// skip header/text lines
			if(!inCardSection){
				if(lists==null){
					lists=new FloatList[lp.terms()];
					for(int i=0; i<lists.length; i++){lists[i]=new FloatList();}
				}
				for(int i=0; i<lp.terms(); i++){lists[i].add(lp.parseFloat(i));}
			}else{
				if(cardLists==null){
					cardLists=new FloatList[lp.terms()];
					for(int i=0; i<cardLists.length; i++){cardLists[i]=new FloatList();}
				}
				for(int i=0; i<lp.terms(); i++){cardLists[i].add(lp.parseFloat(i));}
			}
		}
		bf.close();
		// Build cardinality-indexed table if a #SECTION=cardinality block was found.
		// File layout: col 0=MeanEst, col 1=Samples (skipped), col 2..N=CF values.
		// matrixCard layout: col 0=MeanEst keys, col 1=Mean_cf (type MEAN=1), col 2=HMean_cf, ...
		// We skip the Samples column when building matrixCard so type indices stay aligned.
		if(cardLists!=null && cardLists[0].size()>0){
			final int n=cardLists[0].size();
			// Detect Samples column: present when cardLists has more columns than the occupancy matrix.
			// Occupancy matrix has cols: Slot + CF columns (no Samples).
			// Cardinality matrix has cols: MeanEst + Samples + CF columns.
			// So skip col 1 (Samples) and map file col c to matrixCard col (c-1) for c>=2.
			final boolean hasSamples=(cardLists.length>2);
			final int firstCF=(hasSamples ? 2 : 1); // first CF column in file
			final int numCF=cardLists.length-firstCF; // number of CF columns
			final int matCols=numCF+1; // +1 for MeanEst keys at col 0
			final float[][] cardMat=new float[matCols][];
			// col 0 = MeanEst keys
			cardMat[0]=new float[n];
			for(int i=0; i<n; i++){cardMat[0][i]=cardLists[0].get(i);}
			// cols 1..numCF = CF values (skipping Samples column in file)
			for(int c=0; c<numCF; c++){
				cardMat[c+1]=new float[n];
				for(int i=0; i<n; i++){cardMat[c+1][i]=cardLists[firstCF+c].get(i);}
			}
			lastCardMatrix=cardMat;
			lastCardKeys=cardMat[0]; // col 0 = MeanEst key values
		}
		if(lists==null){return null;}
		// Fix point 0: quadratic extrapolation from points 1, 2, 3.
		// Slot 0 (0 filled buckets) has no meaningful estimate so computeCF returns 1.0 by convention,
		// which is wrong and distorts interpolations near low occupancy.
		if(lists[0].size()>=4){
			for(int col=1; col<lists.length; col++){
				lists[col].set(0,
					(float)quadraticPredict(1,lists[col].get(1), 2,lists[col].get(2), 3,lists[col].get(3), 0));
			}
		}
		float[][] m=interpolate(lists, buckets);
		// If MEAN99 column is absent (9-column legacy file), fall back to MEAN column.
		if(m!=null && m.length<=MEAN99){
			final float[][] m2=new float[MEAN99+1][];
			System.arraycopy(m, 0, m2, 0, m.length);
			m2[MEAN99]=m[MEAN];
			m=m2;
		}
		return m;
	}

	/**
	 * Expands or compresses loaded data into a matrix of exactly buckets+1 slots.
	 * Uses quadratic interpolation (average of two overlapping parabolas) when upsampling,
	 * nearest-slot subsampling when downsampling, and direct copy when sizes match.
	 */
	public static float[][] interpolate(FloatList[] lists, int buckets){
		final int cols=lists.length;
		final int slots=buckets+1;
		final float[][] matrix=new float[cols][slots];
		final int dataSize=lists[0].size();
		if(dataSize==slots){
			// Exact match: copy directly.
			for(int col=0; col<cols; col++){
				for(int s=0; s<slots; s++){matrix[col][s]=lists[col].get(s);}
			}
		}else if(dataSize>slots){
			// More data than slots: stride-subsample.
			final double stride=(double)(dataSize-1)/(slots-1);
			for(int col=0; col<cols; col++){
				for(int s=0; s<slots; s++){
					final int src=(int)Math.round(s*stride);
					matrix[col][s]=lists[col].get(Math.min(src, dataSize-1));
				}
			}
		}else if(dataSize>=3){
			// Fewer data points than slots: quadratic interpolation.
			// For each target slot, average two overlapping parabolas that bracket the position.
			final double invStride=(double)(dataSize-1)/(slots-1);
			for(int col=0; col<cols; col++){
				final FloatList fl=lists[col];
				for(int s=0; s<slots; s++){
					final double pos=s*invStride;// fractional data index
					final int lo=Math.min((int)pos, dataSize-2);
					final int a0=Math.max(0, Math.min(lo-1, dataSize-3));// left triple base
					final int b0=Math.max(0, Math.min(lo,   dataSize-3));// right triple base
					final double predA=quadraticPredict(a0,fl.get(a0), a0+1,fl.get(a0+1), a0+2,fl.get(a0+2), pos);
					final double predB=quadraticPredict(b0,fl.get(b0), b0+1,fl.get(b0+1), b0+2,fl.get(b0+2), pos);
					matrix[col][s]=(float)((predA+predB)*0.5);
				}
			}
		}else if(dataSize==2){
			// Only 2 data points: linear.
			final double invStride=(double)(dataSize-1)/(slots-1);
			for(int col=0; col<cols; col++){
				final FloatList fl=lists[col];
				for(int s=0; s<slots; s++){
					matrix[col][s]=(float)linearPredict(0,fl.get(0), 1,fl.get(1), s*invStride);
				}
			}
		}else{
			// Only 1 data point: constant.
			for(int col=0; col<cols; col++){
				final float v=lists[col].get(0);
				for(int s=0; s<slots; s++){matrix[col][s]=v;}
			}
		}
		return matrix;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Accessors           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Returns the correction factor for the given filled-bucket count and estimator type.
	 * Maps occupancy/buckets to a fractional position in the CF array and uses quadratic
	 * interpolation (average of two parabolas) for accuracy across any bucket count.
	 * @param occupancy Filled bucket count in the DDL
	 * @param buckets   Total bucket count in the DDL
	 * @param type      Estimator type constant (MEAN, HLL, etc.)
	 */
	/**
	 * Returns the correction factor for the given filled-bucket count and estimator type.
	 * Returns 1 if USE_CORRECTION is false, CF_MATRIX is null, or type is LINEAR.
	 * Uses quadratic interpolation (average of two parabolas) for any bucket count.
	 */
	public static float getCF(int occupancy, int buckets, int type){
		if(!USE_CORRECTION || CF_MATRIX==null || type==LINEAR){return 1;}
		final float[] array=CF_MATRIX[type];
		final int maxSlot=array.length-1;
		if(buckets==correctionBuckets){return array[Math.min(occupancy, maxSlot)];}
		final double pos=(double)occupancy*maxSlot/buckets;
		if(maxSlot<2){
			if(maxSlot==0){return array[0];}
			return (float)linearPredict(0,array[0], 1,array[1], pos);
		}
		final int lo=Math.min((int)pos, maxSlot-2);
		final int a0=Math.max(0, Math.min(lo-1, maxSlot-2));
		final int b0=Math.max(0, Math.min(lo,   maxSlot-2));
		final double predA=quadraticPredict(a0,array[a0], a0+1,array[a0+1], a0+2,array[a0+2], pos);
		final double predB=quadraticPredict(b0,array[b0], b0+1,array[b0+1], b0+2,array[b0+2], pos);
		return (float)((predA+predB)*0.5);
	}

	/**
	 * Returns the correction factor using an explicit matrix, for per-class CF support.
	 * Identical logic to getCF(int,int,int) but uses the provided matrix and corrBuckets
	 * instead of the global CF_MATRIX and correctionBuckets.
	 */
	public static float getCF(float[][] matrix, int corrBuckets, int occupancy, int buckets, int type){
		if(!USE_CORRECTION || matrix==null || type==LINEAR || type>=matrix.length){return 1;}
		final float[] array=matrix[type];
		final int maxSlot=array.length-1;
		if(buckets==corrBuckets){return array[Math.min(occupancy, maxSlot)];}
		final double pos=(double)occupancy*maxSlot/buckets;
		if(maxSlot<2){
			if(maxSlot==0){return array[0];}
			return (float)linearPredict(0,array[0], 1,array[1], pos);
		}
		final int lo=Math.min((int)pos, maxSlot-2);
		final int a0=Math.max(0, Math.min(lo-1, maxSlot-2));
		final int b0=Math.max(0, Math.min(lo,   maxSlot-2));
		final double predA=quadraticPredict(a0,array[a0], a0+1,array[a0+1], a0+2,array[a0+2], pos);
		final double predB=quadraticPredict(b0,array[b0], b0+1,array[b0+1], b0+2,array[b0+2], pos);
		return (float)((predA+predB)*0.5);
	}

	/** Returns the correction factor at full occupancy for the given estimator type. */
	public static float getFullCF(int type){
		if(!USE_CORRECTION || CF_MATRIX==null || type==LINEAR){return 1;}
		final float[] array=CF_MATRIX[type];
		return array[array.length-1];
	}

	/**
	 * Three-domain correction factor lookup using raw (uncorrected) mean estimate.
	 * <ul>
	 *   <li>loadFactor &le; 1: occupancy table only (count is still informative).</li>
	 *   <li>1 &lt; loadFactor &lt; 4: log-linear blend of both tables.</li>
	 *   <li>loadFactor &ge; 4: cardinality table only (occupancy saturates).</li>
	 * </ul>
	 * Falls back to occupancy-only if matrixCard is null.
	 * @param rawMeanEst  Raw (uncorrected) mean cardinality estimate
	 */
	public static float getCF(float[][] matrix, int corrBuckets,
			float[][] matrixCard, float[] cardKeys,
			int count, int buckets, double rawMeanEst, int type){
		if(!USE_CORRECTION || matrix==null || type==LINEAR){return 1;}
		final double loadFactor=rawMeanEst/buckets;
		if(matrixCard==null || loadFactor<=1.0){
			return getCF(matrix, corrBuckets, count, buckets, type);
		}else if(loadFactor>=4.0){
			return getCardCF(matrixCard, cardKeys, loadFactor, type);
		}else{
			final double t=Math.log(loadFactor)/Math.log(4.0); // 0 at 1×, 1 at 4×
			final float cfOcc=getCF(matrix, corrBuckets, count, buckets, type);
			final float cfCard=getCardCF(matrixCard, cardKeys, loadFactor, type);
			return (float)((1-t)*cfOcc+t*cfCard);
		}
	}

	/** Binary-search + linear interpolation into the cardinality-indexed CF table. */
	private static float getCardCF(float[][] matrixCard, float[] cardKeys,
			double estCard, int type){
		if(matrixCard==null || cardKeys==null || type<=0 || type>=matrixCard.length){return 1;}
		final float[] cfCol=matrixCard[type];
		final int n=cardKeys.length;
		if(n==0){return 1;}
		if(estCard<=cardKeys[0]){return cfCol[0];}
		if(estCard>=cardKeys[n-1]){return cfCol[n-1];}
		int lo=0, hi=n-1;
		while(lo<hi-1){
			final int mid=(lo+hi)>>>1;
			if(cardKeys[mid]<=(float)estCard){lo=mid;}else{hi=mid;}
		}
		final float k0=cardKeys[lo], k1=cardKeys[hi];
		if(k1==k0){return cfCol[lo];}
		final float frac=(float)((estCard-k0)/(k1-k0));
		return cfCol[lo]+frac*(cfCol[hi]-cfCol[lo]);
	}

	/*--------------------------------------------------------------*/
	/*----------------       Interpolation Math     ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Predicts Y4 from a quadratic (Lagrange) fit through (X1,Y1), (X2,Y2), (X3,Y3).
	 * X values need not be evenly spaced.
	 */
	public static double quadraticPredict(double x1, double y1,
			double x2, double y2, double x3, double y3, double x4){
		final double l1=(x4-x2)*(x4-x3)/((x1-x2)*(x1-x3));
		final double l2=(x4-x1)*(x4-x3)/((x2-x1)*(x2-x3));
		final double l3=(x4-x1)*(x4-x2)/((x3-x1)*(x3-x2));
		return y1*l1+y2*l2+y3*l3;
	}

	/**
	 * Predicts Y3 from a linear fit through (X1,Y1), (X2,Y2).
	 * X values need not be evenly spaced.
	 */
	public static double linearPredict(double x1, double y1,
			double x2, double y2, double x3){
		return y1+(y2-y1)*(x3-x1)/(x2-x1);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	public static String correctionFile="?cardinalityCorrection.tsv.gz";
	public static int correctionBuckets=2048;
	/** When false, getCF always returns 1 (raw estimates, no correction applied). */
	public static boolean USE_CORRECTION=true;

	/** Estimator type constants: index into CF_MATRIX rows.
	 * Matches CF file column order (Slot + CF columns; Hybrid at end as column 10). */
	public static final int OCCUPIED=0, MEAN=1, HMEAN=2, HMEANM=3, GMEAN=4, HLL=5,
		LINEAR=6, MWA=7, MEDCORR=8, MEAN99=9, HYBRID=10;

	/** Per-occupancy correction factor matrix: CF_MATRIX[type][filled_buckets]. Null until initialize() is called. */
	public static float[][] CF_MATRIX=null;

	/**
	 * Cardinality-indexed CF matrix populated by the most recent bipartite loadFile() call.
	 * lastCardMatrix[type][slot] = CF value; lastCardMatrix[0] = MeanEst key values.
	 * Null if the loaded file has no #SECTION=cardinality block.
	 */
	public static float[][] lastCardMatrix=null;
	/** MeanEst key array for binary-search into lastCardMatrix; same reference as lastCardMatrix[0]. */
	public static float[] lastCardKeys=null;

}
