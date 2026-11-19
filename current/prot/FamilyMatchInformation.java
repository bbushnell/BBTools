package prot;

import java.nio.charset.StandardCharsets;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import dna.AminoAcid;
import fileIO.ByteFile;
import parse.LineParser1;
import structures.ByteBuilder;
import shared.Tools;

/** Frozen development background and view-specific scales for exact-match information.
 * Probabilities are indexed by AminoAcid.acidToNumberExtended, canonical codes0..19.
 * This is caller evidence, never an alignment score or a classification gate.
 * @author Sayu
 */
public final class FamilyMatchInformation{
	/** Protect the frozen calibration resource before an ordinary caller opens outputs. */
	public static void protectOutputs(String resource,String[] outputs)throws IOException{
		require(resource!=null && outputs!=null,"Information resource and output roster are required");
		final Path input=Paths.get(resource).toRealPath();
		for(String name:outputs){
			if(!Tools.isOutputFileName(name)){continue;}
			final Path out=Paths.get(name).toFile().getCanonicalFile().toPath();
			if(out.equals(input) || Files.exists(out) && Files.isSameFile(input,out)){
				throw new IllegalArgumentException("Output aliases frozen information resource: "+name);
			}
		}
	}
	/** Small frozen resource. Required rows: cohort, sampling, alpha=1, long_i0,
	 * short_i0, then exactly20 count/uppercase-residue/nonnegative-integer rows.
	 * Cohort and sampling are provenance assertions by the producer, not a validation
	 * that the referenced development population is independent of evaluation data.
	 */
	public static FamilyMatchInformation read(String path){
		require(path!=null && !path.isEmpty(),"Nonzero information mode requires an explicit frozen resource");
		final ByteFile input=ByteFile.makeByteFile1(path,false);
		final LineParser1 lp=new LineParser1('\t');
		final long[] counts=new long[20];final boolean[] seen=new boolean[20];
		double alpha=0,longI0=0,shortI0=0;boolean cohort=false,sampling=false;int n=0;
		try{
			require(Arrays.equals(input.nextLine(),"#bbgap_match_information\t1".getBytes(StandardCharsets.US_ASCII)),
				"Unexpected information resource schema");
			for(byte[] line=input.nextLine();line!=null;line=input.nextLine()){
				lp.set(line);require(lp.terms()>=2,"Information resource has an incomplete row");
				if(lp.termEquals("count",0)){
					require(lp.terms()==3,"A reference count row needs exactly a residue and count");
					final String letter=lp.parseString(1);
					require(letter.length()==1 && letter.charAt(0)>='A' && letter.charAt(0)<='Z',"Count residue must be one uppercase canonical amino acid");
					final int code=AminoAcid.acidToNumberExtended[letter.charAt(0)];
					require(code>=0 && code<20 && !seen[code],"Background needs every canonical residue exactly once");
					// These20 boundary values need overflow checking; Parse's byte loop may wrap.
					counts[code]=Long.parseLong(lp.parseString(2));require(counts[code]>=0,"Negative reference count");seen[code]=true;n++;
				}else{
					require(lp.terms()==2,"Information metadata rows need exactly two fields");
					if(lp.termEquals("cohort",0)){require(!cohort && !lp.parseString(1).trim().isEmpty(),"Missing or repeated development cohort");cohort=true;}
					else if(lp.termEquals("sampling",0)){require(!sampling && !lp.parseString(1).trim().isEmpty(),"Missing or repeated development sampling rule");sampling=true;}
					else if(lp.termEquals("alpha",0)){require(alpha==0,"Repeated pseudocount");alpha=Double.parseDouble(lp.parseString(1));require(alpha==1,"First information policy requires alpha=1");}
					else if(lp.termEquals("long_i0",0)){require(longI0==0,"Repeated LONG normalizer");longI0=positive(lp.parseString(1));}
					else if(lp.termEquals("short_i0",0)){require(shortI0==0,"Repeated SHORT normalizer");shortI0=positive(lp.parseString(1));}
					else{throw new IllegalArgumentException("Unknown information resource field: "+lp.parseString(0));}
				}
			}
		}finally{require(!input.close(),"Information resource read failed");}
		require(n==20 && cohort && sampling && alpha==1,"Incomplete information resource/background provenance");
		return fromCounts(counts,alpha,longI0,shortI0);
	}
	private static double positive(String text){
		final double value=Double.parseDouble(text);
		require(Double.isFinite(value) && value>0,"Information normalizers must be positive and finite");return value;
	}

	public FamilyMatchInformation(double[] probabilities,double longNormalizer,double shortNormalizer){
		require(probabilities!=null && probabilities.length==20,"Information background requires exactly20 canonical probabilities");
		double sum=0;
		this.probabilities=probabilities.clone();
		weights=new double[20];
		for(int i=0;i<20;i++){
			final double p=this.probabilities[i];
			require(Double.isFinite(p) && p>0 && p<1,"Every information background probability must be finite and strictly between0 and1");
			sum+=p;weights[i]=-Math.log(p)/Math.log(2);
		}
		require(Math.abs(sum-1)<=1e-12,"Information background probabilities must sum to1; no silent renormalization");
		longView=new Binding(this,View.LONG,longNormalizer);
		shortView=new Binding(this,View.SHORT,shortNormalizer);
	}

	/** Freeze alpha-smoothed reference-side development counts, not exact-match counts. */
	public static FamilyMatchInformation fromCounts(long[] counts,double alpha,double longNormalizer,double shortNormalizer){
		require(counts!=null && counts.length==20 && Double.isFinite(alpha) && alpha>0,
			"Information background requires20 nonnegative counts and a positive finite pseudocount");
		long total=0;
		for(long n:counts){require(n>=0,"Reference background counts must be nonnegative");total=Math.addExact(total,n);}
		require(total>0,"Empty development reference population cannot define an information background");
		final double denominator=total+20*alpha;
		require(Double.isFinite(denominator),"Smoothed background count overflow");
		final double[] p=new double[20];
		for(int i=0;i<20;i++){p[i]=(counts[i]+alpha)/denominator;}
		return new FamilyMatchInformation(p,longNormalizer,shortNormalizer);
	}

	/** Only call for canonical exact matches, inside the existing traceback walk. */
	public double weight(int code){
		assert(code>=0 && code<20) : "Only canonical exact matches contribute information";
		return weights[code];
	}
	public Binding bind(View view){
		require(view!=null,"Information normalizer must name its classifier view");
		return view==View.LONG ? longView : shortView;
	}
	/** Small reproducibility record in the actual residue encoding, never a fitted estimate. */
	public ByteBuilder recipe(){
		final ByteBuilder out=new ByteBuilder("view\tI0\nLONG\t").append(Double.toString(longView.normalizer))
			.append("\nSHORT\t").append(Double.toString(shortView.normalizer)).append("\ncode\tresidue\tprobability\tbits\n");
		for(char letter:"ACDEFGHIKLMNPQRSTVWY".toCharArray()){
			final int code=AminoAcid.acidToNumberExtended[letter];
			out.append(code).tab().append(letter).tab().append(Double.toString(probabilities[code])).tab().append(Double.toString(weights[code])).nl();
		}
		return out;
	}

	public static final class Binding{
		private Binding(FamilyMatchInformation background_,View view_,double normalizer_){
			require(Double.isFinite(normalizer_) && normalizer_>0,"Information normalizer must be a positive finite frozen development median");
			background=background_;view=view_;normalizer=normalizer_;
		}
		public double ratio(double information){
			require(Double.isFinite(information) && information>=0,"Measured information must be finite and nonnegative");
			final double value=information/normalizer;
			require(Double.isFinite(value),"Information divided by its normalizer overflowed");
			return value;
		}
		public final FamilyMatchInformation background;
		public final View view;
		public final double normalizer;
	}
	public enum View{LONG,SHORT}
	private final double[] weights;
	private final double[] probabilities;
	private final Binding longView,shortView;
	private static void require(boolean condition,String message){if(!condition){throw new IllegalArgumentException(message);}}
}
