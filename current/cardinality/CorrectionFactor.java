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
	/*----------------   SBS CF Table        ----------------*/
	/*--------------------------------------------------------------*/

	/** Resource path for the SBS correction table. */
	public static String sbsFile="?cardinalityCorrectionLC2BitHist.tsv.gz";

	/**
	 * SBS-aware correction table.
	 * SBS_CF_TABLE[filled][stateIdx] = expected distinct elements per bucket.
	 * filled: 0..sbsBuckets (row 0 unused).
	 * stateIdx: 0..10 for the 11 valid states:
	 *   0=0.00, 1=1.00, 2=1.10, 3=2.00, 4=2.01, 5=2.10, 6=2.11,
	 *   7=3+.00, 8=3+.01, 9=3+.10, 10=3+.11
	 * Null until loadSbsTable() is called.
	 */
	public static float[][] SBS_CF_TABLE=null;
	/** Bucket count the SBS table was built for. */
	public static int sbsBuckets=0;

	/** Number of valid states in the loaded SBS table.
	 *  Set by the loader from the file's HistBits header.
	 *  Formula: 3 * (1 << hbits) - 1.  (5 for 1-bit, 11 for 2-bit, 23 for 3-bit.) */
	public static int SBS_STATES=11;

	/** History bits the loaded SBS table was trained for. Set by loader. */
	public static int sbsHbits=2;

	/** Returns the number of valid SBS states for the given history bit width.
	 *  hbits=1 → 5, hbits=2 → 11, hbits=3 → 23. */
	public static int sbsNumStates(int hbits){
		return 3*(1<<hbits)-1;
	}

	/*--------------------------------------------------------------*/
	/*---    SBS Formula Mode (replaces table with closed-form)  ---*/
	/*--------------------------------------------------------------*/

	/** When true, SBS uses closed-form formulas instead of the loaded table.
	 *  The formula y = base + a*L + b*L^2 + c*L^3 where L = ln(B/V)
	 *  is B-independent (verified on B=2048 and B=256) and eliminates
	 *  the need for resource files and bucket-count interpolation.
	 *  Coefficients fitted jointly on B=2048 (32M trials) and B=256 (300K trials).
	 *  Combined R^2 >= 0.999 for all 11 states.  @author Ady */
	public static boolean USE_SBS_FORMULA=false;

	/** Base values per state: minimum distinct count implied by the history pattern. */
	static final int[] SBS_FORMULA_BASE={1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3};

	/** Linear coefficient per state for y = base + a*L + b*L^2 + c*L^3. */
	static final double[] SBS_FORMULA_A={
		0.22105100, 0.12252876, 0.36952900, 0.06133797, 0.30751800,
		0.18526204, 0.43143072, 0.67720200, 0.70123350, 0.64688227, 0.73638900
	};

	/** Quadratic coefficient per state. */
	static final double[] SBS_FORMULA_B={
		0.03969600, 0.00756061, 0.03269600, 0.00236917, 0.02818300,
		0.00912048, 0.03450029, 0.03911000, 0.01746612, 0.01820189, 0.01219500
	};

	/** Cubic coefficient per state (always negative: corrects over-prediction near singularity). */
	static final double[] SBS_FORMULA_C={
		-0.00321100, -0.00065621, -0.00192200, -0.00029760, -0.00176600,
		-0.00068916, -0.00183925, -0.00248900, -0.00112447, -0.00079415, -0.00085800
	};

	/**
	 * Computes expected distinct elements for one bucket using the closed-form formula.
	 * y = base + a*L + b*L^2 + c*L^3 where L = ln(B / V), V = B - filled.
	 * B-independent: same coefficients work for any bucket count.
	 *
	 * @param stateIdx  SBS state index (0-10 for 2-bit history)
	 * @param filled    number of filled buckets
	 * @param B         total bucket count
	 * @return expected distinct elements for this bucket state at this occupancy
	 */
	public static double sbsFormula(int stateIdx, int filled, int B){
		final double V=B-filled;
		if(V<=0){return B;} // singularity cap
		final double L=Math.log((double)B/V);
		return SBS_FORMULA_BASE[stateIdx]
			+SBS_FORMULA_A[stateIdx]*L
			+SBS_FORMULA_B[stateIdx]*L*L
			+SBS_FORMULA_C[stateIdx]*L*L*L;
	}

	/*--------------------------------------------------------------*/
	/*-----------   Mean CF Formula (replaces table)   -------------*/
	/*--------------------------------------------------------------*/

	/** When true, Mean CF uses closed-form 3S2G formula instead of table lookup.
	 *  Formula: CF = a0 + a1*S(lc,c1,w1) + a2*S(lc,c2,w2) + a3*S(lc,c3,w3)
	 *               + g1*G(lc,gc1,gw1) + g2*G(lc,gc2,gw2)
	 *  where lc=log2(card), S(lc,c,w)=(1+tanh((lc-c)/w))/2, G(lc,c,w)=exp(-((lc-c)/w)^2).
	 *  Per-class coefficients fitted on B=2048 (2M-8M estimators).
	 *  Terminal value = a0+a1+a2+a3 (all S→1, all G→0 as card→∞).
	 *  Stable for arbitrarily large cardinalities (no divergence).  @author Eru */
	public static boolean USE_MEAN_CF_FORMULA=false;

	/** When true, enables all available formula modes (SBS, Mean CF, HC CF).
	 *  Individual flags can still override.  @author Eru */
	public static boolean USE_FORMULAS=false;

	/** Per-class 3S2G coefficient arrays.
	 *  Order: {a0, a1, c1, w1, a2, c2, w2, a3, c3, w3, g1, gc1, gw1, g2, gc2, gw2} */
	// DLL4: R²=0.999998, maxErr=0.00045, terminal=0.72096
	static final double[] MCF_DLL4={
		-1.015179396169067, 1.594540127002307, -0.519234303092515, 1.249974670798734,
		-0.179798404560318, 13.595342271669020, 1.208514276422193,
		0.321396896321896, 12.883995972095798, 1.342829006195819,
		0.074973659412974, 4.979261777628070, 5.136222405328894,
		-0.062220148621530, 11.180707724131452, 2.365584912098679};
	// LL6: R²=0.999998, maxErr=0.00055, terminal=0.72095
	static final double[] MCF_LL6={
		-1.708217544869827, 2.319470775861269, -0.790999171032084, 1.301305531370796,
		-0.093784622097852, 9.446809471143753, 1.845004109302430,
		0.203478241095712, 12.034559708649349, 0.843044293514485,
		0.043880560634460, 5.112264691560116, 4.200333607494791,
		-0.042569886954597, 12.461743367200754, 0.845306749606533};
	// DLL3: R²=0.999857, maxErr=0.00211, terminal=0.74492
	static final double[] MCF_DLL3={
		-0.284297841984091, 0.999031401260891, -0.228556730148101, 1.564116182580805,
		0.336114204102007, 12.470057109846664, 1.310101358590010,
		-0.305928673627107, 13.625570459214497, 2.731038958545785,
		-0.056686922608321, 3.023963219338324, 4.533572112498229,
		-0.160434763240275, 11.908935351173893, 4.338474558994728};
	// BDLL3 (cof): R²=0.999926, maxErr=0.00222, terminal=0.74380
	static final double[] MCF_BDLL3_COF={
		-1.126804203356419, 1.792187970103239, -0.701709745229426, 2.227707325270951,
		-0.426430920451125, 13.188527400138502, 2.198089902040401,
		0.504844645771836, 12.551416754053720, 1.473548488950153,
		0.223251242983487, -0.727524575094514, 2.775910307887784,
		-0.127223763102901, 12.095548363616729, 3.804649590273570};
	// BDLL3 (cot): R²=0.999798, maxErr=0.00297, terminal=0.71585
	static final double[] MCF_BDLL3_COT={
		-0.352221510448318, 1.007775775520781, -0.318488804091053, 1.581978485693328,
		0.243791330000102, 10.838990524950210, 1.414837159956376,
		-0.183498884559245, 13.774705127855546, 1.118901895913420,
		-0.221253722109374, 11.901037021849644, 1.853468035403414,
		-0.071580871833045, 10.101440974300839, 2.481107954178233};
	// UDLL6: R²=0.999973, maxErr=0.00071, terminal=0.72091
	static final double[] MCF_UDLL6={
		-0.619543782071058, 1.277187700611610, -0.230145324163336, 1.244575247161501,
		-0.151247923541826, 9.319562771741442, 2.204598330976844,
		0.214509902904495, 11.405161862170619, 0.878558848901602,
		-0.017905605237010, 1.943899332927732, 1.675455792371034,
		-0.097326679808376, 12.002456847327048, 1.105356472907433};

	/** Active Mean CF coefficients. Set by each class's initializeCF(). */
	public static double[] meanCfCoeffs=MCF_DLL4;

	/**
	 * Computes Mean correction factor from closed-form 3S2G formula.
	 * Uses the active coefficient array (set per class).
	 *
	 * @param card  estimated cardinality (typically dlcRawF)
	 * @return multiplicative correction factor for Mean estimate
	 */
	public static double meanCfFormula(double card){
		return meanCfFormula(card, meanCfCoeffs);
	}

	/**
	 * Computes Mean CF from 3S2G formula with explicit coefficient array.
	 * CF = a0 + a1*S(lc,c1,w1) + a2*S(lc,c2,w2) + a3*S(lc,c3,w3)
	 *      + g1*G(lc,gc1,gw1) + g2*G(lc,gc2,gw2)
	 * Terminal value = p[0]+p[1]+p[4]+p[7] (all sigmoids→1, Gaussians→0).
	 *
	 * @param card  estimated cardinality (typically dlcRawF)
	 * @param p     16-element coefficient array
	 * @return multiplicative correction factor
	 */
	public static double meanCfFormula(double card, double[] p){
		if(card<=0){return p[0]+p[1]+p[4]+p[7];} // terminal
		final double lc=Math.log(card)*INV_LOG2;
		final double s1=(1+Math.tanh((lc-p[2])/p[3]))*0.5;
		final double s2=(1+Math.tanh((lc-p[5])/p[6]))*0.5;
		final double s3=(1+Math.tanh((lc-p[8])/p[9]))*0.5;
		final double d1=(lc-p[11])/p[12]; final double g1=Math.exp(-d1*d1);
		final double d2=(lc-p[14])/p[15]; final double g2=Math.exp(-d2*d2);
		return p[0]+p[1]*s1+p[4]*s2+p[7]*s3+p[10]*g1+p[13]*g2;
	}

	private static final double INV_LOG2=1.0/Math.log(2.0);

	/*--------------------------------------------------------------*/
	/*----------------     SBS State Index             ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Maps (nlzBin, histBits) to table column index for arbitrary history bit width.
	 * At NLZ bin k, the top min(k, hbits) history bits are valid; bottom bits must be 0.
	 * Returns -1 for structurally impossible states (nonzero invalid bits).
	 */
	public static int sbsStateIndex(int nlzBin, int histBits, int hbits){
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
	public static int sbsStateIndex(int nlzBin, int histBits){
		return sbsStateIndex(nlzBin, histBits, 2);
	}

	/** Loads the SBS table from the resource file. */
	public static synchronized void loadSbsTable(){
		String path=(sbsFile.indexOf('?')>=0 ? Data.findPath(sbsFile) : sbsFile);
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
		sbsHbits=hbits;
		SBS_STATES=numStates;

		// Build table with era-averaging: table[f] = avg(raw[f], raw[f+1]).
		// Era 0 (f=0): unused. Era 1: raw[1] only. Era B: raw[B] only.
		// Eras 2..B-1: average of raw[f] and raw[f+1].
		// This corrects for mid-era state evolution (snapshots are at era starts).
		SBS_CF_TABLE=new float[buckets+1][numStates];
		for(int f=1; f<=row; f++){
			for(int si=0; si<numStates; si++){
				final float cur=lists[si].get(f-1); // raw[f] is at list index f-1
				if(f<row){
					final float next=lists[si].get(f); // raw[f+1]
					SBS_CF_TABLE[f][si]=(cur+next)*0.5f;
				}else{
					SBS_CF_TABLE[f][si]=cur; // last row: no next
				}
			}
		}
		sbsBuckets=buckets;
	}

	/** Resource path for the SBS multiplier table. */
	public static String sbsMultFile="?cardinalityCorrectionLC2BitHistMult.tsv.gz";

	/**
	 * SBS multiplier table.
	 * SBS_MULT_TABLE[filled][stateIdx] = correction factor (multiplier).
	 * At estimation time: estimate = LC * weightedAvg(CF across filled buckets).
	 * Same layout and loading as SBS_CF_TABLE.
	 */
	public static float[][] SBS_MULT_TABLE=null;
	public static int sbsMultBuckets=0;

	/** Loads the SBS multiplier table. Same format as sum table.
	 *  Silently skips if the resource file does not exist. */
	public static synchronized void loadSbsMultTable(){
		String path=(sbsMultFile.indexOf('?')>=0 ? Data.findPath(sbsMultFile, false) : sbsMultFile);
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
		SBS_MULT_TABLE=new float[buckets+1][numStates];
		for(int f=0; f<row; f++){
			for(int si=0; si<numStates; si++){
				SBS_MULT_TABLE[f+1][si]=lists[si].get(f);
			}
		}
		sbsMultBuckets=buckets;
	}

	/*--------------------------------------------------------------*/
	/*---       HC CF Formula (correction factor for HC)         ---*/
	/*--------------------------------------------------------------*/

	/** HC correction factor: exponential base + linear-amplitude sine/cosine.
	 *  CF(C) = t - a*exp(-b*lc) + (s0+s1*lc)*sin(2*pi*lc) + (c0+c1*lc)*cos(2*pi*lc)
	 *  where lc = log2(C).
	 *  The sine frequency is EXACTLY 1 cycle per octave (structural, not fitted).
	 *  R²=0.99995, maxRes=0.00076, 0 strictly-better violations (B=2048, 2M estimators).
	 *  Returns 1.0 when formulas are disabled.  @author Eru, Ady */
	public static boolean USE_HC_CF_FORMULA=false;

	private static final double HC_CF_T =   0.998484777537405;
	private static final double HC_CF_A = 562.709161557896778;
	private static final double HC_CF_B =   1.011821147415617;
	private static final double HC_CF_S0=   0.000344111534286;
	private static final double HC_CF_S1=  -0.000037079356918;
	private static final double HC_CF_C0=   0.000424944596363;
	private static final double HC_CF_C1=  -0.000035494763272;

	/** Returns the HC correction factor for a given cardinality estimate.
	 *  When USE_HC_CF_FORMULA is false, returns 1.0 (no correction). */
	public static double hcCfFormula(double card){
		if(!USE_HC_CF_FORMULA && !USE_FORMULAS){return 1.0;}
		if(card<=0){return HC_CF_T;} // terminal value
		final double lc=Math.log(card)*INV_LOG2;
		final double base=HC_CF_T-HC_CF_A*Math.exp(-HC_CF_B*lc);
		final double angle=TWO_PI*lc;
		final double s=Math.sin(angle), c=Math.cos(angle);
		return base+(HC_CF_S0+HC_CF_S1*lc)*s+(HC_CF_C0+HC_CF_C1*lc)*c;
	}

	private static final double TWO_PI=2.0*Math.PI;

}
