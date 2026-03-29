package cardinality;

import java.io.File;
import java.util.Arrays;

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
		if(path==null){throw new RuntimeException("Cannot find CF file: "+fname);}
		File f=new File(path);
		if(!f.canRead() || !f.isFile()){
			throw new RuntimeException("Cannot read CF file: "+path);
		}
		// Peek at first line for version detection
		FileFormat ff0=FileFormat.testInput(path, null, false);
		ByteFile bf0=ByteFile.makeByteFile(ff0, 1);
		byte[] firstLine=bf0.nextLine();
		bf0.close();
		if(firstLine!=null){
			final String firstStr=new String(firstLine).trim();
			if(firstStr.startsWith("#VERSION=")){
				final int ver=Integer.parseInt(firstStr.substring(9));
				if(ver>=1){return loadFileV1(path, buckets, ver);}
			}else if(firstStr.startsWith("#Version\t")){
				final int ver=Integer.parseInt(firstStr.substring(9).trim());
				if(ver>=5){return loadFileV1(path, buckets, ver);}
			}
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
	/*----------------          v1 Loading          ----------------*/
	/*--------------------------------------------------------------*/

	/** Column name → type constant mapping for v1 header parsing. */
	private static int colNameToType(String name){
		if("Mean_cf".equals(name)){return MEAN;}
		if("HMean_cf".equals(name)){return HMEAN;}
		if("HMeanM_cf".equals(name)){return HMEANM;}
		if("GMean_cf".equals(name)){return GMEAN;}
		if("DLC3B_cf".equals(name)){return DLC3B;}
		if("HLL_cf".equals(name)) {return HLL;}
		if("Hybrid_cf".equals(name)){return HYBRID;}
		if("HybDLC_cf".equals(name)){return HYBDLC;}
		if("DLC_cf".equals(name)){return DLC;}
		if("HybDLC50_cf".equals(name)){return HYBDLC50;}
		if("DLCBest_cf".equals(name)){return DLCBEST;}
		if("DThHyb_cf".equals(name)){return DTHTHYB;}
		return -1;
	}

	/**
	 * Loads a v1 CF table: unified, cardinality-indexed by raw DLC3B estimate.
	 * Format: #VERSION=N, #ESTIMATOR=DLC3B, header line, then data rows.
	 * Col 0 = DLC3B_est (key), remaining cols = CF values per estimator.
	 * Returns null for CF_MATRIX (v1 has no occupancy-indexed table).
	 */
	private static float[][] loadFileV1(String path, int buckets, int ver){
		tableVersion=ver;
		v1Buckets=0; // reset; will be set from #Buckets header for v5+
		FileFormat ff=FileFormat.testInput(path, null, false);
		ByteFile bf=ByteFile.makeByteFile(ff, 1);
		LineParser1 lp=new LineParser1('\t');
		int[] colTypes=null; // colTypes[fileCol] → type constant; -1 for key column or unknown
		int numTypes=0;
		FloatList keyList=new FloatList();
		FloatList[] cfLists=null; // cfLists[typeIdx] parallel to colTypes

		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			lp.set(line);
			if(lp.startsWith('#')){
				final String s=new String(line).trim();
				if(s.startsWith("#VERSION=") || s.startsWith("#ESTIMATOR=")){continue;}
				// v5+ tab-delimited metadata lines: #Key\tValue
				// Distinguished from column header by absence of "_cf" suffix in content.
				if(ver>=5 && !s.contains("_cf")){
					final int tab=s.indexOf('\t');
					if(tab>0){
						final String key=s.substring(0,tab);
						final String val=s.substring(tab+1).trim();
						if("#Buckets".equals(key)){v1Buckets=Integer.parseInt(val);}
					}
					continue;
				}
				// Header line: #DLC3B_est\tMean_cf\tHMean_cf\t...
				if(colTypes==null){
					final String[] cols=s.split("\t");
					colTypes=new int[cols.length];
					numTypes=0;
					for(int i=0; i<cols.length; i++){
						colTypes[i]=colNameToType(cols[i].replace("#",""));
						if(colTypes[i]>=0){numTypes++;}
					}
					cfLists=new FloatList[cols.length];
					for(int i=0; i<cols.length; i++){cfLists[i]=new FloatList();}
				}
				continue;
			}
			if(!Tools.isNumeric(line[0])){continue;}
			if(colTypes==null){continue;} // no header yet
			keyList.add(lp.parseFloat(0));
			for(int i=1; i<lp.terms() && i<colTypes.length; i++){
				cfLists[i].add(lp.parseFloat(i));
			}
		}
		bf.close();

		if(keyList.size()==0 || colTypes==null){
			throw new RuntimeException("v1 CF table has no data rows: "+path);
		}

		// Build v1Matrix: v1Matrix[type][row] = CF value; v1Keys[row] = DLC3B key
		final int n=keyList.size();
		final int maxType=DTHTHYB+1; // need indices 0..DTHTHYB
		final float[][] mat=new float[maxType][n];
		v1Keys=new float[n];
		for(int i=0; i<n; i++){v1Keys[i]=keyList.get(i);}
		// Fill all types with 1.0 default
		for(int t=0; t<maxType; t++){Arrays.fill(mat[t], 1.0f);}
		// Populate from file columns
		for(int col=1; col<colTypes.length; col++){
			if(colTypes[col]<0 || colTypes[col]>=maxType){continue;}
			final float[] dest=mat[colTypes[col]];
			for(int i=0; i<n; i++){dest[i]=cfLists[col].get(i);}
		}
		v1Matrix=mat;

		// Extend the table to Long.MAX_VALUE by tiling the self-similar upper half.
		// DLL's tier structure repeats every cardinality doubling: CF(2C) ≈ CF(C).
		// Copy CF values from [maxKey/2, maxKey] into [maxKey, 2*maxKey], etc.
		extendV1Table();

		// v1 has no occupancy-indexed table
		lastCardMatrix=null;
		lastCardKeys=null;
		return null;
	}

	/**
	 * Extends v1Keys/v1Matrix to cover arbitrarily high cardinalities by tiling
	 * the self-similar upper half of the CF table. DLL's tier structure is self-similar
	 * across cardinality doublings: after a tier advance, the bucket distribution resets
	 * to the same shape as the previous era. So CF(2C) ≈ CF(C) at corresponding points.
	 * <p>
	 * Finds the upper half of the current key range [maxKey/2, maxKey], then appends
	 * copies with doubled keys for each octave until exceeding MAX_CF_KEY.
	 * At 1% spacing, each half-period has ~70 entries; ~40 doublings to reach Long.MAX_VALUE
	 * adds ~2800 entries. Negligible memory.
	 */
	private static void extendV1Table(){
		if(v1Keys==null || v1Keys.length<2 || v1Matrix==null){return;}
		final int origN=v1Keys.length;
		final float maxKey=v1Keys[origN-1];
		final float halfKey=maxKey*0.5f;

		// Find start of upper half: first key >= maxKey/2
		int halfIdx=origN-1;
		for(int i=0; i<origN; i++){
			if(v1Keys[i]>=halfKey){halfIdx=i; break;}
		}
		final int halfLen=origN-halfIdx; // number of entries in upper half
		if(halfLen<2){return;} // not enough data to tile

		// Count how many doublings we need
		int numDoublings=0;
		double testKey=maxKey;
		while(testKey<MAX_CF_KEY){testKey*=2; numDoublings++;}
		if(numDoublings==0){return;}

		// Build extended arrays
		final int extN=origN+numDoublings*halfLen;
		final float[] extKeys=new float[extN];
		System.arraycopy(v1Keys, 0, extKeys, 0, origN);
		final float[][] extMat=new float[v1Matrix.length][extN];
		for(int t=0; t<v1Matrix.length; t++){
			System.arraycopy(v1Matrix[t], 0, extMat[t], 0, origN);
		}

		// Tile: for each doubling, copy upper-half CF values with doubled keys
		int writeIdx=origN;
		double scale=2.0;
		for(int d=0; d<numDoublings; d++){
			for(int i=0; i<halfLen; i++){
				extKeys[writeIdx]=(float)(v1Keys[halfIdx+i]*scale);
				for(int t=0; t<v1Matrix.length; t++){
					extMat[t][writeIdx]=v1Matrix[t][halfIdx+i];
				}
				writeIdx++;
			}
			scale*=2.0;
		}

		v1Keys=extKeys;
		v1Matrix=extMat;
	}

	/**
	 * v1 correction factor lookup: binary search on raw DLC3B estimate with linear interpolation.
	 * Returns 1 if v1Matrix is null, USE_CORRECTION is false, or type is LINEAR.
	 */
	public static float getV1CF(double dlc3bEst, int type){
		if(!USE_CORRECTION || v1Matrix==null || type==LINEAR || type>=v1Matrix.length){return 1;}
		final float[] cfCol=v1Matrix[type];
		final int n=v1Keys.length;
		if(n==0){return 1;}
		final float key=(float)dlc3bEst;
		if(key<=v1Keys[0]){return cfCol[0];}
		if(key>=v1Keys[n-1]){return cfCol[n-1];}
		int lo=0, hi=n-1;
		while(lo<hi-1){
			final int mid=(lo+hi)>>>1;
			if(v1Keys[mid]<=key){lo=mid;}else{hi=mid;}
		}
		final float k0=v1Keys[lo], k1=v1Keys[hi];
		if(k1==k0){return cfCol[lo];}
		final float frac=(key-k0)/(k1-k0);
		return cfCol[lo]+frac*(cfCol[hi]-cfCol[lo]);
	}

	/*--------------------------------------------------------------*/
	/*----------------     Iterative CF Lookup      ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Binary search + linear interpolation into a single CF column.
	 * Clamps to edge values (no extrapolation beyond table bounds).
	 */
	public static double interpolateCF(double est, float[] cfCol, float[] keys){
		final int n=keys.length;
		if(n==0){return 1.0;}
		final float key=(float)est;
		if(key<=keys[0]){return cfCol[0];}
		if(key>=keys[n-1]){return cfCol[n-1];}
		int lo=0, hi=n-1;
		while(lo<hi-1){
			final int mid=(lo+hi)>>>1;
			if(keys[mid]<=key){lo=mid;}else{hi=mid;}
		}
		final float k0=keys[lo], k1=keys[hi];
		if(k1==k0){return cfCol[lo];}
		final float frac=(key-k0)/(k1-k0);
		return cfCol[lo]+frac*(cfCol[hi]-cfCol[lo]);
	}

	/**
	 * Fully generic iterative CF lookup.
	 * <p>
	 * Refines a cardinality estimate by iteratively looking up the CF for rawEst
	 * using the current best estimate as the lookup key:
	 * <pre>
	 *   bestEst = seedEst           // initial guess (e.g. raw DLC)
	 *   cf = cfCol[bestEst]         // look up CF at current best guess
	 *   bestEst = cf * rawEst       // refine: corrected estimate = CF × raw
	 *   repeat until converged or maxIters
	 * </pre>
	 * Flat CF curves converge in 1 iteration; divergent curves (DLL3) need 2-3.
	 *
	 * @param seedEst   initial cardinality estimate (e.g. raw DLC, CF near 1)
	 * @param rawEst    raw estimate from the target estimator
	 * @param cfCol     CF values for this estimator type
	 * @param keys      cardinality keys for binary search into cfCol
	 * @param maxIters  maximum refinement iterations (1 = single lookup)
	 * @param maxDif    convergence threshold: stop when |newCf - cf| &lt; maxDif
	 * @return correction factor to multiply rawEst by
	 */
	/**
	 * Scaled variant: multiplies lookup key by keyScale before interpolating.
	 * Use keyScale = tableBuckets / currentBuckets when the CF table was built
	 * at a different bucket count than the current instance (v5+ tables).
	 */
	public static double getCF(double seedEst, double rawEst,
			float[] cfCol, float[] keys, int maxIters, double maxDif, double keyScale){
		if(cfCol==null || keys==null || keys.length==0){return 1.0;}
		double bestEst=seedEst*keyScale;
		double cf=interpolateCF(bestEst, cfCol, keys);
		if(TRACE_CF){System.err.println("  getCF iter=0 bestEst="+String.format("%.2f",bestEst)
			+" cf="+String.format("%.8f",cf)+" rawEst="+String.format("%.2f",rawEst));}
		for(int i=1; i<maxIters; i++){
			bestEst=cf*rawEst*keyScale;
			if(bestEst<=0){break;}
			final double newCf=interpolateCF(bestEst, cfCol, keys);
			if(TRACE_CF){System.err.println("  getCF iter="+i+" bestEst="+String.format("%.2f",bestEst)
				+" cf="+String.format("%.8f",newCf)+" rawEst="+String.format("%.2f",rawEst));}
			if(Math.abs(newCf-cf)<maxDif){break;}
			cf=newCf;
		}
		return cf;
	}

	public static double getCF(double seedEst, double rawEst,
			float[] cfCol, float[] keys, int maxIters, double maxDif){
		if(cfCol==null || keys==null || keys.length==0){return 1.0;}
		double bestEst=seedEst;
		double cf=interpolateCF(bestEst, cfCol, keys);
		if(TRACE_CF){System.err.println("  getCF iter=0 bestEst="+String.format("%.2f",bestEst)
			+" cf="+String.format("%.8f",cf)+" rawEst="+String.format("%.2f",rawEst));}
		for(int i=1; i<maxIters; i++){
			bestEst=cf*rawEst;
			if(bestEst<=0){break;}
			final double newCf=interpolateCF(bestEst, cfCol, keys);
			if(TRACE_CF){System.err.println("  getCF iter="+i+" bestEst="+String.format("%.2f",bestEst)
				+" cf="+String.format("%.8f",newCf)+" rawEst="+String.format("%.2f",rawEst));}
			if(Math.abs(newCf-cf)<maxDif){break;}
			cf=newCf;
		}
		return cf;
	}

	/** Set to true to print iterative CF trace to stderr. */
	public static boolean TRACE_CF=false;

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
	 * Matches CF file column order (Slot + CF columns; Hybrid at end as column 10).
	 * v1 adds DLC3B and HYBDLC. */
	public static final int OCCUPIED=0, MEAN=1, HMEAN=2, HMEANM=3, GMEAN=4, HLL=5,
		LINEAR=6, MWA=7, MEDCORR=8, MEAN99=9, HYBRID=10, DLC3B=11, HYBDLC=12, DLC=13, HYBDLC50=14, DLCBEST=15, DTHTHYB=16;

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

	/** v1 table: cardinality-indexed CF keyed by raw DLC3B estimate. */
	public static float[][] v1Matrix=null;
	/** v1 table: DLC3B estimate keys for binary-search into v1Matrix. */
	public static float[] v1Keys=null;
	/** Table version: 0 = legacy bipartite, 1+ = unified DLC3B-indexed. */
	public static int tableVersion=0;
	/** Bucket count used when building v1Matrix; 0 = unknown (no scaling). Set from #Buckets header in v5+. */
	public static int v1Buckets=0;

	/** Maximum CF key value for table extension. Beyond this, CF is clamped to the last value.
	 *  Set to Long.MAX_VALUE cast to double; float precision limits to ~2^63 anyway. */
	public static final double MAX_CF_KEY=(double)Long.MAX_VALUE;

	/*--------------------------------------------------------------*/
	/*----------------   LC History CF Table        ----------------*/
	/*--------------------------------------------------------------*/

	/** Resource path for the LC history correction table. */
	public static String lcHistFile="?cardinalityCorrectionLC2BitHist.tsv.gz";

	/**
	 * LC history-aware correction table.
	 * LC_2BIT_CF_TABLE[filled][stateIdx] = expected distinct elements per bucket.
	 * filled: 0..lcHistBuckets (row 0 unused).
	 * stateIdx: 0..10 for the 11 valid states:
	 *   0=0.00, 1=1.00, 2=1.10, 3=2.00, 4=2.01, 5=2.10, 6=2.11,
	 *   7=3+.00, 8=3+.01, 9=3+.10, 10=3+.11
	 * Null until loadLCHistTable() is called.
	 */
	public static float[][] LC_2BIT_CF_TABLE=null;
	/** Bucket count the LC history table was built for. */
	public static int lcHistBuckets=0;

	/** Number of valid states in the loaded LC history table.
	 *  Set by the loader from the file's HistBits header.
	 *  Formula: 3 * (1 << hbits) - 1.  (5 for 1-bit, 11 for 2-bit, 23 for 3-bit.) */
	public static int LC_HIST_STATES=11;

	/** History bits the loaded LC history table was trained for. Set by loader. */
	public static int lcHistHbits=2;

	/** Returns the number of valid LC history states for the given history bit width.
	 *  hbits=1 → 5, hbits=2 → 11, hbits=3 → 23. */
	public static int lcHistNumStates(int hbits){
		return 3*(1<<hbits)-1;
	}

	/**
	 * Maps (nlzBin, histBits) to table column index for arbitrary history bit width.
	 * At NLZ bin k, the top min(k, hbits) history bits are valid; bottom bits must be 0.
	 * Returns -1 for structurally impossible states (nonzero invalid bits).
	 */
	public static int lcHistStateIndex(int nlzBin, int histBits, int hbits){
		final int bin=Math.min(nlzBin, hbits+1);
		final int validSlots=Math.min(bin, hbits);
		final int invalidBits=hbits-validSlots;
		// Bottom bits must be zero for the state to be valid
		if(invalidBits>0 && (histBits&((1<<invalidBits)-1))!=0){return -1;}
		// Offset for this bin: sum of 2^min(k, hbits) for k=0..bin-1 = (1<<min(bin, hbits+1))-1
		final int offset=(1<<Math.min(bin, hbits+1))-1;
		final int idx=(invalidBits>0) ? (histBits>>>invalidBits) : histBits;
		return offset+idx;
	}

	/** Convenience overload for 2-bit history (the common case). */
	public static int lcHistStateIndex(int nlzBin, int histBits){
		return lcHistStateIndex(nlzBin, histBits, 2);
	}

	/** Loads the LC history table from the resource file. */
	public static synchronized void loadLCHistTable(){
		String path=(lcHistFile.indexOf('?')>=0 ? Data.findPath(lcHistFile) : lcHistFile);
		if(path==null){return;} // silently skip if not found
		FileFormat ff=FileFormat.testInput(path, null, false);
		if(ff==null){return;}
		ByteFile bf=ByteFile.makeByteFile(ff, 1);
		if(bf==null){return;}

		int buckets=0;
		int hbits=2; // default; overridden by #HistBits header
		int numStates=-1; // determined from hbits or first data row
		FloatList[] lists=null; // lists[stateIdx] -> values per filled count
		int row=0;

		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			String s=new String(line).trim();
			if(s.startsWith("#Buckets")){
				String[] parts=s.split("\t");
				if(parts.length>=2){buckets=Integer.parseInt(parts[1].trim());}
				continue;
			}
			if(s.startsWith("#HistBits")){
				String[] parts=s.split("\t");
				if(parts.length>=2){hbits=Integer.parseInt(parts[1].trim());}
				continue;
			}
			if(s.startsWith("#")){continue;}
			// Data row: filled \t val0 \t val1 \t ...
			String[] parts=s.split("\t");
			if(numStates<0){numStates=parts.length-1;} // detect from first data row
			if(parts.length<numStates+1){continue;}
			if(lists==null){
				lists=new FloatList[numStates];
				for(int i=0; i<numStates; i++){lists[i]=new FloatList();}
			}
			for(int i=0; i<numStates; i++){
				lists[i].add(Float.parseFloat(parts[i+1]));
			}
			row++;
		}
		bf.close();

		if(lists==null || buckets==0){return;}

		// Set globals from loaded header
		lcHistHbits=hbits;
		LC_HIST_STATES=numStates;

		// Build table with era-averaging: table[f] = avg(raw[f], raw[f+1]).
		// Era 0 (f=0): unused. Era 1: raw[1] only. Era B: raw[B] only.
		// Eras 2..B-1: average of raw[f] and raw[f+1].
		// This corrects for mid-era state evolution (snapshots are at era starts).
		LC_2BIT_CF_TABLE=new float[buckets+1][numStates];
		for(int f=1; f<=row; f++){
			for(int si=0; si<numStates; si++){
				final float cur=lists[si].get(f-1); // raw[f] is at list index f-1
				if(f<row){
					final float next=lists[si].get(f); // raw[f+1]
					LC_2BIT_CF_TABLE[f][si]=(cur+next)*0.5f;
				}else{
					LC_2BIT_CF_TABLE[f][si]=cur; // last row: no next
				}
			}
		}
		lcHistBuckets=buckets;
	}

	/** Resource path for the LC history multiplier table. */
	public static String lcHistMultFile="?cardinalityCorrectionLC2BitHistMult.tsv.gz";

	/**
	 * LC history multiplier table.
	 * LC_2BIT_MULT_TABLE[filled][stateIdx] = correction factor (multiplier).
	 * At estimation time: estimate = LC * weightedAvg(CF across filled buckets).
	 * Same layout and loading as LC_2BIT_CF_TABLE.
	 */
	public static float[][] LC_2BIT_MULT_TABLE=null;
	public static int lcHistMultBuckets=0;

	/** Loads the LC history multiplier table. Same format as sum table.
	 *  Silently skips if the resource file does not exist. */
	public static synchronized void loadLCHistMultTable(){
		String path=(lcHistMultFile.indexOf('?')>=0 ? Data.findPath(lcHistMultFile, false) : lcHistMultFile);
		if(path==null){return;}
		FileFormat ff=FileFormat.testInput(path, null, false);
		if(ff==null){return;}
		ByteFile bf=ByteFile.makeByteFile(ff, 1);
		if(bf==null){return;}

		int buckets=0;
		int numStates=-1;
		FloatList[] lists=null;
		int row=0;

		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			String s=new String(line).trim();
			if(s.startsWith("#Buckets")){
				String[] parts=s.split("\t");
				if(parts.length>=2){buckets=Integer.parseInt(parts[1].trim());}
				continue;
			}
			if(s.startsWith("#")){continue;}
			String[] parts=s.split("\t");
			if(numStates<0){numStates=parts.length-1;}
			if(parts.length<numStates+1){continue;}
			if(lists==null){
				lists=new FloatList[numStates];
				for(int i=0; i<numStates; i++){lists[i]=new FloatList();}
			}
			for(int i=0; i<numStates; i++){
				lists[i].add(Float.parseFloat(parts[i+1]));
			}
			row++;
		}
		bf.close();

		if(lists==null || buckets==0){return;}
		LC_2BIT_MULT_TABLE=new float[buckets+1][numStates];
		for(int f=0; f<row; f++){
			for(int si=0; si<numStates; si++){
				LC_2BIT_MULT_TABLE[f+1][si]=lists[si].get(f);
			}
		}
		lcHistMultBuckets=buckets;
	}

}
