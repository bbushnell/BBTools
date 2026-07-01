package var2;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntList;

/**
 * Balances a labeled training-vector TSV (as produced by {@link VcfToTrainingVectors})
 * into train/validation splits for neural-network training, using stratified negative
 * sampling so that rare-but-hard false positives are not drowned out by common easy ones.
 *
 * <p>Strategy (deterministic given a seed):
 * <ol>
 * <li>Keep ALL positives (label &gt; 0.5).
 * <li>Choose a target negative count so that positives are {@code posfraction} of the output
 *     (default 0.3 -&gt; 30% positive / 70% negative).
 * <li>ENRICHED negatives: sample a fixed per-category quota across category axes built from
 *     real {@link VectorUMP45} columns -- variant type x {depth, score, event-length}, depth
 *     alone, and (optionally) the artifact axes where false positives concentrate:
 *     type x homopolymer (slippage), allele fraction (low-VAF errors), strand ratio, and mapq.
 *     New (artifact) axes are sampled at {@code newweight} times the base quota -- a smaller
 *     fraction, not an equal split.
 * <li>REPRESENTATIVE negatives: fill the remainder up to the target with a random sample of
 *     nontrivial negatives (composite score &gt; 0.01) not already chosen.  This makes the
 *     enriched set a minority and pins the final ratio to {@code posfraction}.
 * <li>Shuffle, split off {@code valfraction} as validation, and write both.
 * </ol>
 *
 * <p>With {@code noscore} (default true) the composite-score column (29) is zeroed on output,
 * forcing the network to learn its own scoring rather than copy the existing score.
 *
 * <p>With {@code depth1count}&gt;0, depth-1 variants (allelic depth AD==1, recovered from vector
 * column 9) are handled as a SEPARATE category: they are pulled out of the standard pos/neg pools
 * (so the stratified balancing above applies only to AD&ge;2), then sampled at their natural
 * true/false frequency up to {@code depth1count} rows -- enriched to a {@code depth1floor} positive
 * fraction (default 1%) only when the natural positive rate is below it -- and appended to the
 * balanced AD&ge;2 set before the train/val split.  {@code depth1count=0} (default) restores the
 * original behavior (AD==1 balanced together with everything else).
 *
 * <p>This is the Java, reproducible counterpart of the Python prototype balancer; the column
 * indices below are bound to {@link VectorUMP45}'s 33-feature layout.
 *
 * @author UMP45
 * @date June 7, 2026
 */
public class BalanceVectors {

	public static void main(String[] args){
		Timer t=new Timer();
		BalanceVectors x=new BalanceVectors(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	/**
	 * Parses command-line arguments and validates I/O.
	 * @param args Command-line arguments
	 */
	public BalanceVectors(String[] args){
		{
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=", 2);
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("in")){in=b;}
			else if(a.equals("outtrain") || a.equals("train")){outTrain=b;}
			else if(a.equals("outval") || a.equals("val") || a.equals("validate")){outVal=b;}
			else if(a.equals("posfraction") || a.equals("ratio")){posFraction=Double.parseDouble(b);}
			else if(a.equals("valfraction")){valFraction=Double.parseDouble(b);}
			else if(a.equals("newweight")){newWeight=Double.parseDouble(b);}
			else if(a.equals("enrich")){enrich=Parse.parseBoolean(b);}
			else if(a.equals("newcats") || a.equals("artifactcats")){newCats=Parse.parseBoolean(b);}
			else if(a.equals("noscore")){noScore=Parse.parseBoolean(b);}
			else if(a.equals("seed")){seed=Long.parseLong(b);}
			else if(a.equals("quota")){quotaOverride=Integer.parseInt(b);}
			else if(a.equals("depth1count") || a.equals("d1count")){depth1Count=Integer.parseInt(b);}
			else if(a.equals("depth1floor") || a.equals("d1floor")){depth1Floor=Double.parseDouble(b);}
			else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		if(parser.in1!=null && in==null){in=parser.in1;}
		assert(in!=null) : "Missing in= (input vector TSV)";
		assert(outTrain!=null) : "Missing outtrain=";
		assert(outVal!=null) : "Missing outval=";
		assert(posFraction>0 && posFraction<1) : "posfraction must be in (0,1): "+posFraction;
		assert(valFraction>=0 && valFraction<1) : "valfraction must be in [0,1): "+valFraction;

		if(!ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}
		ffin=FileFormat.testInput(in, FileFormat.TXT, null, true, true);
		ffoutTrain=FileFormat.testOutput(outTrain, FileFormat.TXT, null, true, true, false, false);
		ffoutVal=FileFormat.testOutput(outVal, FileFormat.TXT, null, true, true, false, false);
	}

	/**
	 * Loads vectors, performs stratified balancing, and writes train/val splits.
	 * @param t Timer for elapsed-time reporting
	 */
	void process(Timer t){
		final Random randy=new Random(seed);

		//Load: separate positives (kept whole) from negatives (sampled)
		final ArrayList<byte[]> pos=new ArrayList<byte[]>();
		final ArrayList<byte[]> neg=new ArrayList<byte[]>();

		//Depth-1 (AD==1) variants, gathered separately when depth1Count>0 (sampled, not balanced)
		final ArrayList<byte[]> d1pos=new ArrayList<byte[]>();
		final ArrayList<byte[]> d1neg=new ArrayList<byte[]>();

		//Category buckets: key -> list of negative indices
		final HashMap<String, IntList> existing=new HashMap<String, IntList>();
		final HashMap<String, IntList> newc=new HashMap<String, IntList>();

		byte[] header=null;
		ByteFile bf=ByteFile.makeByteFile(ffin);
		final float[] f=new float[VectorUMP45.DIMS+1];
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			if(line[0]=='#'){header=line; continue;}//#dims header
			linesIn++;
			parseFields(line, f);
			final float label=f[f.length-1];
			if(depth1Count>0 && recoverAD(f)==1){//Route depth-1 aside; not run through standard balancing
				if(label>0.5f){d1pos.add(line);}else{d1neg.add(line);}
				continue;
			}
			if(label>0.5f){
				pos.add(line);
			}else{
				final int idx=neg.size();
				neg.add(line);
				if(enrich){addCategories(f, idx, existing, newc);}
			}
		}
		bf.close();
		outstream.println("Loaded "+linesIn+" vectors: "+pos.size()+" positive, "+neg.size()+" negative.");

		//Target negatives so that positives are posFraction of the total
		final long targetNeg=Math.round(pos.size()*(1-posFraction)/posFraction);

		//Per-category base quota: spread half the target across the existing categories
		//(matches the prototype's scale); new (artifact) categories get newWeight x this.
		final int baseQuota=quotaOverride>0 ? quotaOverride :
			(int)Math.max(1, (targetNeg/2)/Math.max(1, existing.size()));
		final int newQuota=(int)Math.max(1, Math.round(baseQuota*newWeight));

		//ENRICHED: sample distinct negative indices across all category buckets
		final boolean[] chosen=new boolean[neg.size()];
		long enrichedCount=0;
		enrichedCount+=sampleBuckets(existing, baseQuota, chosen, randy);
		if(newCats){enrichedCount+=sampleBuckets(newc, newQuota, chosen, randy);}
		outstream.println("Categories: "+existing.size()+" existing"+
			(newCats ? " + "+newc.size()+" artifact" : "")+
			"; quota existing="+baseQuota+(newCats ? ", artifact="+newQuota : "")+
			"; enriched negatives="+enrichedCount);

		//REPRESENTATIVE: fill the remainder with random negatives not already chosen by enrichment
		long fill=targetNeg-enrichedCount;
		if(fill>0){
			final IntList pool=new IntList();
			//TODO: Possible bug [var2/BalanceVectors#001] (DOC/QUESTION->UMP45) - class javadoc says the
			//representative fill samples "nontrivial negatives (composite score > 0.01)", but there is NO such
			//filter here: pool = ALL unchosen negatives. Either restore the filter (only meaningful when
			//vectors were built with includescore=t, else score col is 0 and it'd exclude everything) or fix
			//the javadoc. NN-pipeline-owned: flag, don't change balancing behavior unilaterally.
			for(int i=0; i<neg.size(); i++){
				if(!chosen[i]){pool.add(i);}
			}
			final int take=(int)Math.min(fill, pool.size());
			final int[] arr=pool.toArray();
			partialShuffle(arr, take, randy);
			for(int i=0; i<take; i++){chosen[arr[i]]=true;}
			outstream.println("Representative fill: "+take+" (of "+fill+" requested, pool="+pool.size()+").");
		}

		//Combine kept positives + chosen negatives (the depth-2+ balanced portion)
		final ArrayList<byte[]> combined=new ArrayList<byte[]>(pos.size()+(int)targetNeg+depth1Count);
		combined.addAll(pos);
		long negKept=0;
		for(int i=0; i<neg.size(); i++){
			if(chosen[i]){combined.add(neg.get(i)); negKept++;}
		}
		outstream.println("Depth-2+ balanced: "+pos.size()+" pos + "+negKept+" neg = "+combined.size()+
			" ("+String.format("%.3f", pos.size()/(double)Math.max(1, combined.size()))+" positive)"+
			(noScore ? "  [score col zeroed]" : ""));

		//Depth-1 (AD==1): sampled independently at natural pos/neg frequency, enriched to the depth1Floor
		//positive fraction ONLY if the natural positive rate is below it; NOT run through standard balancing.
		if(depth1Count>0){
			final long avail=d1pos.size()+d1neg.size();
			final double natFrac=(avail>0 ? d1pos.size()/(double)avail : 0);
			long npos=Math.round(depth1Count*natFrac);
			if(natFrac<depth1Floor){npos=Math.round(depth1Floor*depth1Count);}//enrich only when natural < floor
			npos=Math.min(npos, d1pos.size());
			final long baseNneg=depth1Count-npos;//baseline depth-1 negatives for this dataset (the floor)
			//Backfill: when the AD>=2 negative pool fell short of the 30/70 target (low depth), make up the
			//shortfall with EXTRA depth-1 negatives so positivity stays constant across depths.  Only negatives
			//backfill -- never positives -- so the net does not learn that low depth is more likely positive.
			final long deficit=Math.max(0, targetNeg-negKept);
			long nneg=Math.min(baseNneg+deficit, d1neg.size());
			sampleInto(combined, d1pos, (int)npos, randy);
			sampleInto(combined, d1neg, (int)nneg, randy);
			outstream.println("Depth-1 (AD==1): "+d1pos.size()+" pos + "+d1neg.size()+" neg avail; sampled "+
				npos+" pos + "+nneg+" neg (floor "+baseNneg+" + AD>=2 deficit "+deficit+
				"; natFrac="+String.format("%.4f", natFrac)+
				(natFrac<depth1Floor ? " -> pos enriched to "+String.format("%.0f%%", depth1Floor*100)+" floor" : "")+").");
		}

		//Shuffle the combined set (depth-2+ balanced + depth-1 sampled) and split off validation
		final byte[][] arr=combined.toArray(new byte[0][]);
		fullShuffle(arr, randy);

		final long total=arr.length;
		final long valN=Math.round(total*valFraction);
		final long trainN=total-valN;
		outstream.println("Total training set: "+total+" rows; train="+trainN+", val="+valN+".");

		//Write train then val
		writeSplit(ffoutTrain, header, arr, 0, (int)trainN);
		writeSplit(ffoutVal, header, arr, (int)trainN, arr.length);

		t.stop();
		outstream.println(Tools.timeLinesBytesProcessed(t, linesIn, 0, 8));
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Category Logic        ----------------*/
	/*--------------------------------------------------------------*/

	/** Adds negative index idx to every category bucket it belongs to. */
	private void addCategories(float[] f, int idx, HashMap<String, IntList> existing,
			HashMap<String, IntList> newc){
		final String vt=typeStr(f);
		final String db=depthBin(f);
		bucket(existing, vt+"_"+db, idx);
		bucket(existing, vt+"_"+scoreBin(f), idx);
		bucket(existing, vt+"_"+eventLenBin(f), idx);
		bucket(existing, db, idx);
		if(newCats){
			bucket(newc, vt+"_"+hpBin(f), idx);//type x homopolymer (slippage)
			bucket(newc, afBin(f), idx);       //allele fraction (low-VAF errors)
			bucket(newc, srBin(f), idx);       //strand ratio
			bucket(newc, mqBin(f), idx);       //mapq
		}
	}

	private static void bucket(HashMap<String, IntList> map, String key, int idx){
		IntList list=map.get(key);
		if(list==null){list=new IntList(); map.put(key, list);}
		list.add(idx);
	}

	private static String typeStr(float[] f){
		if(f[IDX_TYPE_SUB]>0.5f){return "SUB";}
		if(f[IDX_TYPE_INS]>0.5f){return "INS";}
		if(f[IDX_TYPE_DEL]>0.5f){return "DEL";}
		return "OTHER";
	}
	private static String depthBin(float[] f){
		int d=Math.max(1, (int)(Math.pow(2, f[IDX_DEPTH]*8)-1));
		if(d<10){return "d"+d;}
		if(d<100){return "d"+((d/10)*10)+"-"+((d/10)*10+9);}
		return "d100+";
	}
	private static String scoreBin(float[] f){
		return "s"+Math.min(19, (int)(f[IDX_SCORE]*100));
	}
	private static String eventLenBin(float[] f){
		int e=Math.max(0, (int)(Math.pow(2, f[IDX_EVENT_LEN]*8)-1));
		if(e<=1){return "len1";}
		if(e<=5){return "len2-5";}
		if(e<=20){return "len6-20";}
		return "len20+";
	}
	private static String hpBin(float[] f){
		float c=f[IDX_HP];
		int hpc=c>0 ? Math.round(1f/c-1) : 0;//invert 1/(hpc+1)
		if(hpc<=0){return "hp0";}
		if(hpc<=2){return "hp1-2";}
		if(hpc<=5){return "hp3-5";}
		return "hp6+";
	}
	private static String afBin(float[] f){
		float af=f[IDX_AF];
		if(af<0.1f){return "af00-10";}
		if(af<0.3f){return "af10-30";}
		if(af<0.6f){return "af30-60";}
		return "af60+";
	}
	private static String srBin(float[] f){
		float sr=f[IDX_STRAND];
		if(sr<0.1f){return "sr00-10";}
		if(sr<0.3f){return "sr10-30";}
		return "sr30+";
	}
	private static String mqBin(float[] f){
		float mq=f[IDX_MAPQ]*40;
		if(mq<10){return "mq00-10";}
		if(mq<20){return "mq10-20";}
		if(mq<40){return "mq20-40";}
		return "mq40+";
	}

	/*--------------------------------------------------------------*/
	/*----------------        Sampling Helpers      ----------------*/
	/*--------------------------------------------------------------*/

	/** Samples up to quota distinct indices from each bucket; marks them chosen. Returns newly chosen count. */
	private static long sampleBuckets(HashMap<String, IntList> map, int quota, boolean[] chosen, Random randy){
		long added=0;
		for(IntList list : map.values()){
			final int[] a=list.toArray();
			final int take=Math.min(quota, a.length);
			partialShuffle(a, take, randy);
			for(int i=0; i<take; i++){
				if(!chosen[a[i]]){chosen[a[i]]=true; added++;}
			}
		}
		return added;
	}

	/** Partial Fisher-Yates: randomizes the first k elements of a in place (distinct draw). */
	private static void partialShuffle(int[] a, int k, Random randy){
		final int n=a.length;
		for(int i=0; i<k && i<n; i++){
			int j=i+randy.nextInt(n-i);
			int tmp=a[i]; a[i]=a[j]; a[j]=tmp;
		}
	}

	/** Full Fisher-Yates shuffle of a byte[][] in place. */
	private static void fullShuffle(byte[][] a, Random randy){
		for(int i=a.length-1; i>0; i--){
			int j=randy.nextInt(i+1);
			byte[] tmp=a[i]; a[i]=a[j]; a[j]=tmp;
		}
	}

	/** Recovers integer allelic depth (AD) from the encoded VectorUMP45 column 9 (= log2(AD+1)/8). */
	private static int recoverAD(float[] f){
		return (int)Math.round(Math.pow(2, f[IDX_AD]*8.0)-1);
	}

	/** Randomly draws k distinct rows from src and appends them to dest (all of src when k>=src.size()). */
	private static void sampleInto(ArrayList<byte[]> dest, ArrayList<byte[]> src, int k, Random randy){
		final int n=src.size();
		if(k>=n){dest.addAll(src); return;}
		final int[] idx=new int[n];
		for(int i=0; i<n; i++){idx[i]=i;}
		partialShuffle(idx, k, randy);
		for(int i=0; i<k; i++){dest.add(src.get(idx[i]));}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Output             ----------------*/
	/*--------------------------------------------------------------*/

	/** Writes header + rows [from,to) to ff, zeroing the score column when noScore is set. */
	private void writeSplit(FileFormat ff, byte[] header, byte[][] rows, int from, int to){
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		ByteBuilder bb=new ByteBuilder(1<<16);
		if(header!=null){bb.append(header).nl();}
		for(int i=from; i<to; i++){
			if(noScore){appendZeroedScore(bb, rows[i]);}
			else{bb.append(rows[i]);}
			bb.nl();
			if(bb.length>=(1<<15)){bsw.print(bb); bb.clear();}
		}
		if(bb.length>0){bsw.print(bb);}
		errorState|=bsw.poisonAndWait();
	}

	/** Appends a vector line to bb, replacing the IDX_SCORE-th tab field with "0". */
	private static void appendZeroedScore(ByteBuilder bb, byte[] line){
		int field=0, i=0;
		while(i<line.length){
			if(field==IDX_SCORE){
				bb.append((byte)'0');//Substitute zeroed score
				while(i<line.length && line[i]!='\t'){i++;}//Skip original value
				if(i<line.length){bb.append((byte)'\t'); i++; field++;}
			}else{
				bb.append(line[i]);
				if(line[i]=='\t'){field++;}
				i++;
			}
		}
	}

	/** Parses just the value of the col-th tab field of a vector line. */
	private static float fieldAt(byte[] line, int col){
		int a=0, field=0;
		for(int i=0; i<=line.length; i++){
			if(i==line.length || line[i]=='\t'){
				if(field==col){return a<i ? Parse.parseFloat(line, a, i) : 0f;}
				field++; a=i+1;
			}
		}
		return 0f;
	}

	/** Parses a tab-delimited vector line into dest[0..DIMS] (DIMS features + label). */
	private static void parseFields(byte[] line, float[] dest){
		int a=0, field=0;
		for(int i=0; i<=line.length; i++){
			if(i==line.length || line[i]=='\t'){
				if(field<dest.length){dest[field]=(a<i ? Parse.parseFloat(line, a, i) : 0f);}
				field++;
				a=i+1;
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Column Indices       ----------------*/
	/*--------------------------------------------------------------*/

	/** VectorUMP45 column indices used for stratification (see VectorUMP45.makeVector). */
	private static final int IDX_TYPE_SUB=1, IDX_TYPE_INS=2, IDX_TYPE_DEL=3;
	private static final int IDX_DEPTH=8, IDX_AD=9, IDX_AF=10, IDX_MAPQ=12, IDX_EVENT_LEN=21;
	private static final int IDX_STRAND=22, IDX_HP=28, IDX_SCORE=29;

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in=null, outTrain=null, outVal=null;
	private final FileFormat ffin, ffoutTrain, ffoutVal;

	private double posFraction=0.3;//Target positive fraction of the output
	private double valFraction=0.1;//Fraction held out for validation
	private double newWeight=0.3;  //Artifact-category quota relative to existing categories
	private boolean enrich=true;   //Use enriched category sampling at all
	private boolean newCats=true;  //Include the artifact axes (homopolymer/AF/strand/mapq)
	private boolean noScore=true;  //Zero the composite-score column on output
	private long seed=42;
	private int quotaOverride=-1;  //If >0, overrides the computed per-existing-category quota
	private int depth1Count=0;     //If >0, gather AD==1 variants separately and sample this many (NOT balanced)
	private double depth1Floor=0.01;//Enrich depth-1 positives to this fraction only if the natural rate is below it

	private long linesIn=0;
	private PrintStream outstream=System.err;
	private boolean errorState=false;
}
