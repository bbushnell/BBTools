package prot;

import java.nio.charset.StandardCharsets;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import fileIO.ByteFile;
import parse.LineParser1;
import structures.ByteBuilder;
import shared.Tools;

/** One frozen coefficient for the eight-class LONG/SHORT16K5 calibration.
 * Provenance identities are inherited from the reviewed execution packet, not
 * recomputed here. This artifact requires deployed negative replay before use
 * as an accepted policy; it is not itself evidence of replay or positive accuracy.
 * @author Keqing
 */
public final class FamilyLengthPolicy{

	/** Reject an output alias before a caller opens its writers, including existing hard links. */
	public static void protectOutputs(Path resource,String[] outputs)throws IOException{
		require(resource!=null && outputs!=null,"Length policy and output roster are required");
		final Path input=resource.toRealPath();
		for(String name:outputs){
			if(!Tools.isOutputFileName(name)){continue;}
			final Path out=Paths.get(name).toFile().getCanonicalFile().toPath();
			if(out.equals(input) || Files.exists(out) && Files.isSameFile(input,out)){
				throw new IllegalArgumentException("Output aliases frozen length policy: "+name);
			}
		}
	}

	public FamilyLengthPolicy(String reference_,String referencePin_,String source_,Requirement[] classes){
		require(FamilyLengthAcceptance.MIN_OVERLAP==28 && FamilyLengthAcceptance.MAX_OVERLAP==100 &&
			FamilyLengthAcceptance.BLOSUM_UNIT==.005,"Changed arithmetic constants require a new policy schema");
		path(reference_);referenceIdentity(referencePin_);
		require(source_!=null && source_.matches("[0-9a-f]{40}"),"Policy needs its frozen software source commit");
		require(classes!=null && classes.length==8,"One shared policy requires exactly eight complete classes");
		final Requirement[] ordered=new Requirement[8];int bits=0;
		for(Requirement row:classes){
			require(row!=null && ordered[row.index]==null,"Every view/class must occur exactly once");
			ordered[row.index]=row;bits=Math.max(bits,row.amplitudeBits);
		}
		for(Requirement row:ordered){require(row!=null,"A complete calibration class is missing");}
		reference=reference_;referencePin=referencePin_;source=source_;requirements=ordered;amplitudeBits=bits;
	}

	/** Reload with an explicit caller-side reference binding; no basename fallback. */
	public static FamilyLengthPolicy read(Path file,String expectedReference,String expectedPin){
		require(file!=null,"A frozen policy path is required");
		final ByteFile input=ByteFile.makeByteFile1(file.toString(),false);
		final LineParser1 lp=new LineParser1('\t');
		final Map<String,String> metadata=new HashMap<String,String>();
		final Requirement[] rows=new Requirement[8];int count=0;
		try{
			require(Arrays.equals(input.nextLine(),HEADER.getBytes(StandardCharsets.US_ASCII)),"Unexpected length-policy schema");
			for(byte[] line=input.nextLine();line!=null;line=input.nextLine()){
				lp.set(line);require(lp.terms()>=2,"Incomplete length-policy row");
				if(lp.termEquals("class",0)){
					require(lp.terms()==13 && count<8,"Policy class rows need exactly twelve values and eight rows");
					rows[count++]=new Requirement(lp.parseString(1),lp.parseString(2),
						Integer.parseInt(lp.parseString(3)),Integer.parseInt(lp.parseString(4)),Integer.parseInt(lp.parseString(5)),
						Integer.parseInt(lp.parseString(6)),Integer.parseInt(lp.parseString(7)),Integer.parseInt(lp.parseString(8)),
						Integer.parseInt(lp.parseString(9)),lp.parseString(10),lp.parseString(11),lp.parseString(12));
				}else{
					require(lp.terms()==2 && metadata.put(lp.parseString(0),lp.parseString(1))==null,
						"Policy metadata must have unique two-column rows");
				}
			}
		}finally{require(!input.close(),"Length-policy read failed");}
		require(count==8 && metadata.size()==11,"Incomplete or unknown length-policy metadata/classes");
		equal(metadata,"formula","inverse_overlap_v1");equal(metadata,"min_overlap","28");equal(metadata,"anchor_overlap","100");
		equal(metadata,"blosum_unit","0.005");equal(metadata,"verification",VERIFICATION);
		equal(metadata,"purpose","frozen_calibration_candidate_requires_deployed_replay");
		final FamilyLengthPolicy result=new FamilyLengthPolicy(value(metadata,"tailreference"),
			value(metadata,"tailreferencesha80"),value(metadata,"software_source"),rows);
		result.requireReference(expectedReference,expectedPin);
		require(Integer.parseInt(value(metadata,"amplitude_bits"))==result.amplitudeBits,
			"Serialized coefficient must equal the maximum of all eight class requirements");
		require(Float.floatToRawIntBits(Float.parseFloat(value(metadata,"amplitude")))==result.amplitudeBits,
			"Decimal coefficient must reload to the exact stored float bits");
		return result;
	}

	/** Final small artifact; all class minima and provenance remain available for audit. */
	public ByteBuilder encode(){
		final ByteBuilder out=new ByteBuilder(HEADER).nl().append("formula\tinverse_overlap_v1\n")
			.append("min_overlap\t28\nanchor_overlap\t100\nblosum_unit\t0.005\n")
			.append("tailreference\t").append(reference).nl().append("tailreferencesha80\t").append(referencePin).nl()
			.append("software_source\t").append(source).nl().append("verification\t").append(VERIFICATION).nl()
			.append("purpose\tfrozen_calibration_candidate_requires_deployed_replay\n")
			.append("amplitude_bits\t").append(amplitudeBits).nl()
			.append("amplitude\t").append(Float.toString(amplitude())).nl();
		for(Requirement r:requirements){
			out.append("class\t").append(view(r.index)).tab().append(cohort(r.index)).tab()
				.append(r.queries).tab().append(r.baseline).tab().append(r.budget).tab().append(r.immutable).tab()
				.append(r.amplitudeBits).tab().append(r.accepted).tab().append(r.precedingAccepted).tab()
				.append(r.input).tab().append(r.inputPin).tab().append(r.reportDirectory).nl();
		}
		return out;
	}
	public float amplitude(){return FamilyLengthAcceptance.amplitude(amplitudeBits);}
	public Requirement requirement(int index){checkIndex(index);return requirements[index];}
	public void requireReference(String expectedReference,String expectedPin){
		require(reference.equals(expectedReference) && referencePin.equals(expectedPin),
			"Length policy belongs to a different frozen joint reference binding");
	}

	/** Counts below are the predeclared complete populations, not fitted length bins. */
	public static int expectedQueries(int index){checkIndex(index);return QUERIES[index];}
	public static int expectedBaseline(int index){checkIndex(index);return BASELINES[index];}
	public static String view(int index){checkIndex(index);return index<4 ? "long" : "short";}
	public static String cohort(int index){checkIndex(index);return COHORTS[index%4];}
	public static int index(String view,String cohort){
		final int offset="long".equals(view) ? 0 : "short".equals(view) ? 4 : -1;
		require(offset>=0,"Policy view must be long or short");
		for(int i=0;i<4;i++){if(COHORTS[i].equals(cohort)){return offset+i;}}
		throw new IllegalArgumentException("Policy cohort must be k5, shadow, reversed or shuffled");
	}

	/** Validated local class minimum, with inherited input and software-packet provenance.
	 * Counts at this local amplitude are not counts at the final shared amplitude.
	 */
	public static final class Requirement{
		public Requirement(String view,String cohort,int queries_,int baseline_,int budget_,int immutable_,
			int bits_,int accepted_,int preceding_,String input_,String inputPin_,String reports_){
			index=FamilyLengthPolicy.index(view,cohort);
			require(queries_==QUERIES[index] && baseline_==BASELINES[index],
				"Class population/original baseline differs from the frozen study");
			require(budget_==FamilyLengthFit.allowed100ppm(queries_),"Class budget must be exactly100 accepts per million rounded down");
			require(bits_>=0 && bits_<=0x7f7fffff,"A class requirement must be a finite nonnegative float encoding");
			require(immutable_>=0 && immutable_<=accepted_ && accepted_<=budget_,
				"Immutable and total local accepts must fit the exact class budget");
			require(bits_==0 ? preceding_==accepted_ : preceding_>budget_ && preceding_<=queries_,
				"The preceding float must fail the class budget; a zero coefficient has no predecessor");
			path(input_);path(reports_);identity(inputPin_);
			queries=queries_;baseline=baseline_;budget=budget_;immutable=immutable_;amplitudeBits=bits_;
			accepted=accepted_;precedingAccepted=preceding_;input=input_;inputPin=inputPin_;reportDirectory=reports_;
		}
		public final int index,queries,baseline,budget,immutable,amplitudeBits,accepted,precedingAccepted;
		public final String input,inputPin,reportDirectory;
	}
	private static void checkIndex(int index){require(index>=0 && index<8,"Class index must be0..7");}
	private static void path(String value){
		require(value!=null && !value.isEmpty() && value.indexOf('\t')<0 && value.indexOf('\n')<0 && value.indexOf('\r')<0 &&
			Paths.get(value).isAbsolute(),"Policy provenance needs absolute single-field paths");
	}
	/** Structural joint manifests may explicitly lack a measured identity. Preserve
	 * that marker and exact path binding; never manufacture a digest for the caller. */
	private static void referenceIdentity(String value){
		require("unmeasured".equals(value) || value!=null && value.matches("[0-9a-f]{20}"),
			"Joint reference identity requires an inherited sha80 or explicit unmeasured marker");
	}
	private static void identity(String value){require(value!=null && value.matches("[0-9a-f]{20}"),"Policy provenance requires an inherited sha80");}
	private static String value(Map<String,String> map,String key){final String v=map.get(key);require(v!=null,"Missing policy field: "+key);return v;}
	private static void equal(Map<String,String> map,String key,String expected){require(expected.equals(value(map,key)),"Unexpected policy field: "+key);}
	private static void require(boolean condition,String message){if(!condition){throw new IllegalArgumentException(message);}}
	public final String reference,referencePin,source;
	public final int amplitudeBits;
	private final Requirement[] requirements;
	private static final String HEADER="#bbgap_length_policy\t1";
	private static final String VERIFICATION="inherited_identity_and_structural_stability";
	private static final String[] COHORTS={"k5","shadow","reversed","shuffled"};
	private static final int[] QUERIES={1000000,1000000,1000000,1000000,1246631,1000000,567080,567080};
	private static final int[] BASELINES={0,57,11,1,0,10,27,1};
}
