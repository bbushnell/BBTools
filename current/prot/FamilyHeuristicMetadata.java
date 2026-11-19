package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import ddl.DDLFamilySearch.StructuralInputs;
import fileIO.ByteFile;
import parse.LineParser1;
import structures.IntHashMap;

/** Immutable roster-bound member medians and hydrophobicity means for heuristic ranking.
 * @author Yelan
 */
public final class FamilyHeuristicMetadata{

	public FamilyHeuristicMetadata(FamilyTailCascade.SearchConfig config,Path medians_,String medianPin_,
			Path familyHydro_,String hydroPin_){
		this(ids(config), ranks(config), members(config), medians_, medianPin_, familyHydro_, hydroPin_);
	}

	/** Explicit inherited-identity mode retains the same live-roster and numeric checks. */
	public FamilyHeuristicMetadata(FamilyTailCascade.SearchConfig config,Path medians_,String medianPin_,
			Path familyHydro_,String hydroPin_,FamilyCalibration.VerificationMode verificationMode_){
		this(ids(config), ranks(config), members(config), medians_, medianPin_, familyHydro_, hydroPin_, verificationMode_);
	}

	public FamilyHeuristicMetadata(FamilyMemberPreparation.Prepared[] prepared,Path medians_,String medianPin_,
			Path familyHydro_,String hydroPin_){
		this(ids(prepared), ranks(prepared), members(prepared), medians_, medianPin_, familyHydro_, hydroPin_);
	}

	/** Resident preparation propagates its own policy; no process-global verification switch. */
	public FamilyHeuristicMetadata(FamilyMemberPreparation.Prepared[] prepared,Path medians_,String medianPin_,
			Path familyHydro_,String hydroPin_,FamilyCalibration.VerificationMode verificationMode_){
		this(ids(prepared), ranks(prepared), members(prepared), medians_, medianPin_, familyHydro_,
			hydroPin_, verificationMode_);
	}

	public FamilyHeuristicMetadata(String[] ids_,int[] ranks_,long[] members_,Path medians_,String medianPin_,
			Path familyHydro_,String hydroPin_){
		this(ids_, ranks_, members_, medians_, medianPin_, familyHydro_, hydroPin_, FamilyCalibration.VerificationMode.SHA80);
	}

	/** Compact family records already carry the validated live-order median length and mean KD. */
	public FamilyHeuristicMetadata(FamilyTailCascade.SearchConfig config,double[] lengths_,double[] hydrophobicities_){
		this(ids(config), ranks(config), members(config), lengths_, hydrophobicities_);
	}

	/** Structural checks bind measured metadata to live families, not newly measured byte identity. */
	public FamilyHeuristicMetadata(String[] ids_,int[] ranks_,long[] members_,Path medians_,String medianPin_,
			Path familyHydro_,String hydroPin_,FamilyCalibration.VerificationMode verificationMode_){
		require(verificationMode_!=null, "Member metadata needs an explicit verification policy");
		verificationMode=verificationMode_;
		structural=verificationMode==FamilyCalibration.VerificationMode.STRUCTURAL ? new StructuralInputs() : null;
		require(ids_!=null && ranks_!=null && members_!=null && medians_!=null,
			"Member metadata binding requires live ids, ranks, members and medians");
		require((familyHydro_==null)==(hydroPin_==null), "Member KD binding must name both the table and pin or neither");
		require(ids_.length>0 && ranks_.length==ids_.length && members_.length==ids_.length,
			"Member metadata arrays must name the same nonempty live roster");
		ids=ids_.clone(); ranks=ranks_.clone(); members=members_.clone();
		validateRoster();
		medians=medians_.toAbsolutePath().normalize(); medianPin=medianPin_;
		familyHydro=familyHydro_==null ? null : familyHydro_.toAbsolutePath().normalize(); hydroPin=hydroPin_;
		verify(); lengths=readMedians(); hydrophobicities=familyHydro==null ? null : readHydrophobicities(); verify();
	}

	FamilyHeuristicMetadata(String[] ids_,int[] ranks_,long[] members_,
			double[] lengths_,double[] hydrophobicities_){
		require(ids_!=null && ranks_!=null && members_!=null && lengths_!=null,
			"In-memory member metadata binding requires live ids, ranks, members and medians");
		require(ids_.length>0 && ranks_.length==ids_.length && members_.length==ids_.length &&
			lengths_.length==ids_.length && (hydrophobicities_==null || hydrophobicities_.length==ids_.length),
			"In-memory member metadata arrays must name the same nonempty live roster");
		verificationMode=null; structural=null; medians=familyHydro=null; medianPin=hydroPin=null;
		ids=ids_.clone(); ranks=ranks_.clone(); members=members_.clone();
		validateRoster();
		lengths=lengths_.clone(); hydrophobicities=hydrophobicities_==null ? null : hydrophobicities_.clone();
		for(double value : lengths){
			require(Double.isFinite(value) && value>=40 && value<=200 && value*2==Math.rint(value*2),
				"Member median must be a finite integer or half integer in40..200; no imputation");
		}
		if(hydrophobicities!=null){
			for(double value : hydrophobicities){FamilyHeuristicRanking.hydrophobicitySimilarity(value, value);}
		}
	}

	public int families(){return ids.length;}
	public int rank(int ordinal){
		require(ordinal>=0 && ordinal<ranks.length, "Family rank needs a live family ordinal");
		return ranks[ordinal];
	}
	public long members(int ordinal){
		require(ordinal>=0 && ordinal<members.length, "Family member count needs a live family ordinal");
		return members[ordinal];
	}
	public String id(int ordinal){
		require(ordinal>=0 && ordinal<ids.length, "Family ID needs a live family ordinal");
		return ids[ordinal];
	}
	public boolean hasHydrophobicity(){return hydrophobicities!=null;}
	public double median(int ordinal){
		require(ordinal>=0 && ordinal<lengths.length, "Member median diagnostics require a live family ordinal");
		return lengths[ordinal];
	}
	public double hydrophobicity(int ordinal){
		require(hydrophobicities!=null && ordinal>=0 && ordinal<hydrophobicities.length,
			"Member KD diagnostics require a bound hydrophobicity table and live family ordinal");
		return hydrophobicities[ordinal];
	}

	/** Recompute digests only in the default policy; structural mode checks metadata stability. */
	public void verify(){
		if(medians==null){return;}
		require(validIdentity(medianPin) &&
			Files.isRegularFile(medians) && !Files.isSymbolicLink(medians) &&
			(verificationMode==FamilyCalibration.VerificationMode.STRUCTURAL ||
			medianPin.equals(DigestSuffix.file(medians.toString()))), "Member median binding changed: "+medians);
		if(familyHydro!=null){
			require(validIdentity(hydroPin) &&
				Files.isRegularFile(familyHydro) && !Files.isSymbolicLink(familyHydro) &&
				(verificationMode==FamilyCalibration.VerificationMode.STRUCTURAL ||
				hydroPin.equals(DigestSuffix.file(familyHydro.toString()))), "Member KD binding changed: "+familyHydro);
		}
		if(structural!=null){
			try{structural.track(medians);if(familyHydro!=null){structural.track(familyHydro);}}
			catch(java.io.IOException | IllegalStateException failure){
				throw new IllegalArgumentException("Member metadata changed while bound",failure);
			}
		}
	}

	/** Unmeasured is permitted only under the explicit structural policy. */
	private boolean validIdentity(String pin){
		return pin!=null && (pin.matches("[0-9a-f]{20}") ||
			verificationMode==FamilyCalibration.VerificationMode.STRUCTURAL && pin.equals("unmeasured"));
	}

	private void validateRoster(){
		final HashSet<String> seenIds=new HashSet<String>();
		final IntHashMap seenRanks=new IntHashMap(ranks.length*2+1);
		for(int i=0; i<ids.length; i++){
			require(ids[i]!=null && !ids[i].isEmpty() && ids[i].indexOf('\t')<0 && ids[i].indexOf('\n')<0 &&
				ids[i].indexOf('\r')<0 && seenIds.add(ids[i]), "Live family IDs must be unique for metadata binding");
			require(ranks[i]>=0 && seenRanks.get(ranks[i])<0,
				"Live physical ranks must be unique and nonnegative for metadata binding");
			seenRanks.put(ranks[i], i);
			require(members[i]>0, "Metadata binding requires live nonempty family member counts");
		}
	}

	private double[] readMedians(){
		final HashMap<String,Integer> ordinals=new HashMap<String,Integer>();
		for(int i=0; i<ids.length; i++){
			require(ordinals.put(ids[i], i)==null, "Live family IDs must be unique for median binding");
		}
		final double[] values=new double[ids.length];
		final boolean[] seen=new boolean[values.length];
		final ByteFile input=ByteFile.makeByteFile1(medians.toString(), false);
		final LineParser1 lp=new LineParser1('\t');
		int count=0;
		try{
			require(Arrays.equals(input.nextLine(), MEDIAN_HEADER), "Unexpected member median table schema");
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				lp.set(line); require(lp.terms()==6, "Member median rows need six fields");
				final Integer ordinal=ordinals.get(lp.parseString(0));
				require(ordinal!=null && !seen[ordinal] && lp.parseLong(1)==members[ordinal],
					"Median identity and count must match the exact live family roster once");
				final double value=lp.parseDouble(3);
				require(Double.isFinite(value) && value>=40 && value<=200 && value*2==Math.rint(value*2),
					"Member median must be a finite integer or half integer in40..200; no imputation");
				values[ordinal]=value; seen[ordinal]=true; count++;
			}
		}finally{require(!input.close(), "Member median table read failed");}
		require(count==ids.length, "Member median table is incomplete"); return values;
	}

	private double[] readHydrophobicities(){
		final IntHashMap ordinals=new IntHashMap(ids.length*2+1);
		for(int i=0; i<ranks.length; i++){ordinals.put(ranks[i], i);}
		final double[] values=new double[ids.length];
		final boolean[] seen=new boolean[values.length];
		final HashSet<String> rowIds=new HashSet<String>();
		final ByteFile input=ByteFile.makeByteFile1(familyHydro.toString(), false);
		final LineParser1 lp=new LineParser1('\t');
		int rows=0, count=0;
		try{
			require(Arrays.equals(input.nextLine(), HYDRO_HEADER), "Unexpected member KD table schema");
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				lp.set(line); require(lp.terms()==4, "Member KD rows need rank, identity, count and mean");
				final int rank=lp.parseInt(0); final String id=lp.parseString(1); final long n=lp.parseLong(2);
				require(rank==rows++ && n>=0 && !id.isEmpty() && rowIds.add(id),
					"Member KD rows must use dense nonnegative physical ranks and unique family identities");
				final int ordinal=ordinals.get(rank);
				if(n==0){
					require(ordinal<0 && lp.termEquals("NA", 3), "Empty KD rows must be non-live and explicitly undefined");
				}else{
					require(ordinal>=0 && !seen[ordinal] && id.equals(ids[ordinal]) && n==members[ordinal],
						"Member KD rank, identity and count must match the exact live family roster");
					values[ordinal]=lp.parseDouble(3);
					FamilyHeuristicRanking.hydrophobicitySimilarity(values[ordinal], values[ordinal]);
					seen[ordinal]=true; count++;
				}
			}
		}finally{require(!input.close(), "Member KD table read failed");}
		require(count==ids.length, "Member KD table is incomplete; no consensus fallback"); return values;
	}

	private static String[] ids(FamilyTailCascade.SearchConfig config){
		require(config!=null && config.families!=null, "Metadata binding needs a classifier live roster");
		final String[] out=new String[config.families.length];
		for(int i=0; i<out.length; i++){out[i]=family(config, i).id;}
		return out;
	}

	private static int[] ranks(FamilyTailCascade.SearchConfig config){
		require(config!=null && config.families!=null, "Metadata binding needs a classifier live roster");
		final int[] out=new int[config.families.length];
		for(int i=0; i<out.length; i++){out[i]=family(config, i).rank;}
		return out;
	}

	private static long[] members(FamilyTailCascade.SearchConfig config){
		require(config!=null && config.families!=null, "Metadata binding needs a classifier live roster");
		final long[] out=new long[config.families.length];
		for(int i=0; i<out.length; i++){out[i]=family(config, i).members;}
		return out;
	}

	private static FamilyReference family(FamilyTailCascade.SearchConfig config, int ordinal){
		final FamilyReference family=config.families[ordinal];
		require(family!=null, "Metadata binding needs a complete classifier live roster");
		return family;
	}

	private static String[] ids(FamilyMemberPreparation.Prepared[] prepared){
		require(prepared!=null, "Metadata binding needs prepared live families");
		final String[] out=new String[prepared.length];
		for(int i=0; i<out.length; i++){out[i]=prepared(prepared, i).id;}
		return out;
	}

	private static int[] ranks(FamilyMemberPreparation.Prepared[] prepared){
		require(prepared!=null, "Metadata binding needs prepared live families");
		final int[] out=new int[prepared.length];
		for(int i=0; i<out.length; i++){out[i]=prepared(prepared, i).rank;}
		return out;
	}

	private static long[] members(FamilyMemberPreparation.Prepared[] prepared){
		require(prepared!=null, "Metadata binding needs prepared live families");
		final long[] out=new long[prepared.length];
		for(int i=0; i<out.length; i++){out[i]=prepared(prepared, i).members;}
		return out;
	}

	private static FamilyMemberPreparation.Prepared prepared(FamilyMemberPreparation.Prepared[] prepared, int ordinal){
		final FamilyMemberPreparation.Prepared family=prepared[ordinal];
		require(family!=null, "Metadata binding needs complete prepared live families");
		return family;
	}

	private static void require(boolean valid,String why){if(!valid){throw new IllegalArgumentException(why);}}

	private final Path medians, familyHydro;
	private final FamilyCalibration.VerificationMode verificationMode;
	private final StructuralInputs structural;
	private final String medianPin, hydroPin;
	private final String[] ids;
	private final int[] ranks;
	private final long[] members;
	private final double[] lengths, hydrophobicities;
	private static final byte[] MEDIAN_HEADER=("family\tmembers\tlength_p10\tlength_median"+
		"\tlength_p90\tcount_under_200").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] HYDRO_HEADER=("family_rank\tfamily_id\tN\thydrophobicity_mean")
		.getBytes(StandardCharsets.US_ASCII);
}
