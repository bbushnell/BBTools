package prot;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Map;

import ddl.DDLPentamerBlacklist;
import ddl.ProteinKmerGeometry;
import ddl.ProteinBlacklistPolicy;
import ddl.ProteinBlacklistPolicy.Descriptor;

/** Explicit blacklist policy shared by resident preparation and classifier loading.
 * An unmeasured label describes provenance; it is never proof that two lists match.
 * Resident callers additionally retain and compare the immutable filter object.
 * @author Keqing
 */
public final class FamilyBlacklistBinding{

	/** New producers materialize their exact geometry's selected policy in the recipe. */
	public static void applyDefaults(Map<String,String> options){
		final Descriptor selected=officialPolicy(options);
		if(selected==null){return;}
		int fields=0;for(String key:POLICY_KEYS){if(options.containsKey(key)){fields++;}}
		require(fields==0 || fields==POLICY_KEYS.length,"Partial blacklist policy cannot be replaced by defaults");
		if(!options.containsKey("blacklist")){
			require(!options.containsKey("blacklistsha80") && fields==0,"Blacklist metadata has no input path");
			options.put("blacklist",selected.resource().toString());
			options.put("blacklistsha80",selected.sha80);
		}
		if(!hasPolicy(options)){
			options.put("blacklistk",Integer.toString(selected.filterK));options.put("blacklistmode",selected.mode.label());
			options.put("blacklistwords",Integer.toString(selected.words));
			options.put("blacklistprovenance",selected.provenance);
		}
		requireOfficialOptions(options);
	}

	/** Existing file-backed models must declare their actual filter; never silently retrofit them. */
	public static void requireOfficialOptions(Map<String,String> options){
		final Descriptor selected=officialPolicy(options);
		if(selected==null){return;}
		require(options.containsKey("blacklist") && selected.sha80.equals(options.get("blacklistsha80")) &&
			Integer.toString(selected.filterK).equals(options.get("blacklistk")) && selected.mode.label().equals(options.get("blacklistmode")) &&
			Integer.toString(selected.words).equals(options.get("blacklistwords")) && selected.provenance.equals(options.get("blacklistprovenance")),
			selected.label+" policy is required; rebuild historical profiles and clustering instead of mixing filters");
	}

	static boolean officialGeometry(Map<String,String> options){return officialPolicy(options)!=null;}
	/** Metadata-only admission shared by file-backed classifier and caller loading. */
	public static boolean structuralFileGeometry(Map<String,String> options){return officialPolicy(options)!=null || explicitReduced(options);}
	/** Descriptor selection is metadata-only and must not load an unrelated resource. */
	private static Descriptor officialPolicy(Map<String,String> options){
		require(options!=null,"Blacklist policy selection requires explicit geometry options");
		final String k=options.get("k");
		if(!"5".equals(k) && !"6".equals(k)){return null;}
		return ProteinBlacklistPolicy.policy(ProteinKmerGeometry.make(Integer.parseInt(k),
			options.containsKey("aabits") ? options.get("aabits") : options.get("aapartition")));
	}

	/** Check only configuration shape and domains; no files or digests are read. */
	public static void validateOptions(Map<String,String> options,FamilyCalibration.VerificationMode verification){
		require(options!=null && verification!=null,"Blacklist options require an explicit verification policy");
		requireOfficialOptions(options);
		final boolean enabled=options.containsKey("blacklist"),structural=verification==FamilyCalibration.VerificationMode.STRUCTURAL;
		int policyFields=0;for(String key:POLICY_KEYS){if(options.containsKey(key)){policyFields++;}}
		require(enabled==options.containsKey("blacklistsha80"),"Blacklist path and declared identity must be supplied together");
		require(enabled || policyFields==0,"Blacklist policy fields cannot silently enable a missing list");
		if(!enabled){return;}
		require(value(options,"blacklist").length()>0,"An enabled blacklist needs a nonempty path");
		final String identity=value(options,"blacklistsha80");
		require(identity.matches("[0-9a-f]{20}") || structural && identity.equals("unmeasured"),
			"Only structural verification permits an explicitly unmeasured blacklist identity");
		require(policyFields==0 || policyFields==POLICY_KEYS.length,"Blacklist length, mode, count and provenance must be specified together");
		require(!structural || policyFields==POLICY_KEYS.length,"Structural filtered preparation requires the complete explicit blacklist policy");
		if(policyFields==0){return;}
		final int k=positive(options,"blacklistk"),words=positive(options,"blacklistwords");
		final DDLPentamerBlacklist.MatchMode mode=DDLPentamerBlacklist.MatchMode.parse(value(options,"blacklistmode"));
		if(mode==DDLPentamerBlacklist.MatchMode.REDUCED_LITERAL){
			final Descriptor selected=officialPolicy(options);
			require(k==positive(options,"k"),"Reduced blacklist and sketch word lengths must agree");
			require((k==5 || k==6) && (selected!=null && selected.mode==DDLPentamerBlacklist.MatchMode.REDUCED_LITERAL ||
				structural && explicitReduced(options)),
				"Reduced filtering requires an official policy or explicit structural full-word monomer/nonoverlapping K6 dimer geometry");
		}else{
			require((k==5 || k==6) && (mode==DDLPentamerBlacklist.MatchMode.LITERAL ||
				mode==DDLPentamerBlacklist.MatchMode.PREFIX && k==5),
				"Supported exclusion policies are raw literal K5/K6 or raw K5 prefix");
		}
		requireOfficialOptions(options);
		require(words>0,"A required blacklist must contain at least one unique word");
		final String provenance=value(options,"blacklistprovenance");
		require(!provenance.trim().isEmpty() && provenance.indexOf('\t')<0 && provenance.indexOf('\n')<0 && provenance.indexOf('\r')<0,
			"Blacklist provenance must be explicit and fit one metadata field");
	}

	/** Load once using the requested verification policy; legacy measured calls keep K5-prefix defaults. */
	public static DDLPentamerBlacklist load(Map<String,String> options,FamilyCalibration.VerificationMode verification)throws Exception{
		validateOptions(options,verification);
		requireOfficialOptions(options);
		if(!options.containsKey("blacklist")){return null;}
		final Path path=Paths.get(value(options,"blacklist"));
		final String identity=value(options,"blacklistsha80");
		final DDLPentamerBlacklist filter;
		final Descriptor selected=officialPolicy(options);
		if(selected!=null && !(verification==FamilyCalibration.VerificationMode.STRUCTURAL && explicitDimer(options))){
			final DDLPentamerBlacklist loaded=selected.load(path,verification==FamilyCalibration.VerificationMode.STRUCTURAL);
			// Structural identity alone is not equivalence: resolve compares all immutable keys.
			filter=ProteinBlacklistPolicy.resolve(selected.geometry,loaded);
		}else if(verification==FamilyCalibration.VerificationMode.STRUCTURAL){
			require(path.isAbsolute(),"Structural blacklist paths must be explicit absolute paths");
			final DDLPentamerBlacklist.MatchMode mode=DDLPentamerBlacklist.MatchMode.parse(value(options,"blacklistmode"));
			filter=mode==DDLPentamerBlacklist.MatchMode.REDUCED_LITERAL ?
				DDLPentamerBlacklist.loadReducedStructural(path,identity,positive(options,"blacklistwords"),
					explicitGeometry(options),value(options,"blacklistprovenance")) :
				DDLPentamerBlacklist.loadStructural(path,identity,positive(options,"blacklistwords"),positive(options,"blacklistk"),
					mode,value(options,"blacklistprovenance"));
		}else if(hasPolicy(options)){
			filter=DDLPentamerBlacklist.load(path,identity,positive(options,"blacklistk"),
				DDLPentamerBlacklist.MatchMode.parse(value(options,"blacklistmode")),value(options,"blacklistprovenance"));
		}else{filter=DDLPentamerBlacklist.load(path,identity);}
		verifyLoaded(options,filter,verification);return filter;
	}

	/** Verify a directly transferred object without reloading its path or computing a digest. */
	public static void verifyLoaded(Map<String,String> options,DDLPentamerBlacklist filter,FamilyCalibration.VerificationMode verification){
		validateOptions(options,verification);
		requireOfficialOptions(options);
		require(options.containsKey("blacklist")== (filter!=null),"Resident blacklist presence differs from the declared policy");
		if(filter==null){return;}
		require(filter.sha80().equals(value(options,"blacklistsha80")),"Resident blacklist declared identity differs from the recipe");
		if(hasPolicy(options)){
			require(filter.size()==positive(options,"blacklistwords") && filter.wordLength()==positive(options,"blacklistk") &&
				filter.matchMode()==DDLPentamerBlacklist.MatchMode.parse(value(options,"blacklistmode")) &&
				filter.provenance().equals(value(options,"blacklistprovenance")),
				"Resident blacklist count, word length, mode or provenance differs from the recipe");
		}else{
			require(filter.wordLength()==5 && filter.matchMode()==DDLPentamerBlacklist.MatchMode.PREFIX,
				"A historical blacklist recipe denotes raw K5-prefix filtering only");
		}
		final Descriptor selected=officialPolicy(options);
		if(selected!=null && !(verification==FamilyCalibration.VerificationMode.STRUCTURAL && explicitDimer(options))){
			ProteinBlacklistPolicy.resolve(selected.geometry,filter);
		}
		else if(filter.matchMode()==DDLPentamerBlacklist.MatchMode.REDUCED_LITERAL){
			filter.validate(explicitGeometry(options));
		}
	}

	private static boolean explicitDimer(Map<String,String> options){
		final String alphabet=options.get("aabits");
		if(!"6".equals(options.get("k")) || alphabet==null || !alphabet.startsWith("dimer:1,3,5:")){return false;}
		return ProteinKmerGeometry.make(6,alphabet).name().startsWith("dimer6+1,3,5+");
	}

	/** Explicit full-word geometry controls reduced projection; no implicit alphabet. */
	private static ProteinKmerGeometry explicitGeometry(Map<String,String> options){
		final String alphabet=options.containsKey("aabits") ? options.get("aabits") : options.get("aapartition");
		require(alphabet!=null,"Reduced filtering requires an explicit alphabet");
		return ProteinKmerGeometry.make(positive(options,"k"),alphabet);
	}
	private static boolean explicitReduced(Map<String,String> options){
		if(!options.containsKey("aabits") && !options.containsKey("aapartition")){return false;}
		final ProteinKmerGeometry geometry=explicitGeometry(options);
		return geometry.isReduced() || geometry.k()==6 && geometry.name().startsWith("dimer6+1,3,5+");
	}

	/** Geometry compatibility is independent of whether a list is measured or structurally verified. */
	public static void validateGeometry(DDLPentamerBlacklist filter,ProteinKmerGeometry geometry){
		require(geometry!=null,"Blacklist geometry validation needs an explicit sketch geometry");
		if(filter!=null){filter.validate(geometry);}
	}
	public static boolean hasPolicy(Map<String,String> options){
		require(options!=null,"Blacklist policy lookup requires options");return options.containsKey("blacklistk");
	}
	private static int positive(Map<String,String> options,String key){
		final String text=value(options,key);final int number;
		try{number=Integer.parseInt(text);}catch(NumberFormatException e){throw new IllegalArgumentException("Blacklist "+key+" must be a positive integer: "+text,e);}
		require(number>0,"Blacklist "+key+" must be positive: "+number);return number;
	}
	private static String value(Map<String,String> options,String key){
		final String value=options.get(key);require(value!=null && !value.isEmpty(),"Missing blacklist option: "+key);return value;
	}
	private static void require(boolean valid,String why){if(!valid){throw new IllegalArgumentException(why);}}
	private static final String[] POLICY_KEYS={"blacklistk","blacklistmode","blacklistwords","blacklistprovenance"};
}
