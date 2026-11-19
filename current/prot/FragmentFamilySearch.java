package prot;

import java.util.Arrays;
import java.util.HashMap;

import fileIO.ByteFile;
import parse.LineParser1;
import prok.ContigEdgeFragment;

/** Standalone, explicitly bound fragment search. No caller bonus or production default.
 * Raw-floor5 admission, DDL top50 plus ties, then raw hits + stored monomer cosine,
 * ascending family-rank ties and at most50 full-query glocal alignments. All shortlisted
 * alignments are measured, even after acceptance, for uncensored calibration observation.
 * @author Qiqi
 */
public final class FragmentFamilySearch{
	FragmentFamilySearch(final FamilyTailCascade.SearchConfig config_, final HashMap<String, String> policy){
		assert(config_!=null && config_.top==50 && config_.ddlReject==5 && !config_.scalarGates &&
			!config_.index.usesPruning() && !config_.index.usesDirectCompletion()) :
			"Fragment observations need the complete raw index and separately bound acceptance gates";
		config=config_; alignment=new FragmentAlignment(Double.parseDouble(policy.get("identityfloor")), Double.parseDouble(policy.get("blosumfloor")));
		policyId=policy.get("policyid"); referenceSha80=policy.get("familyreferencessha80"); profilesSha80=policy.get("ddlprofilessha80");
	}

	/** The package header and all input pins are checked before and after loading.
	 * A caller-tail manifest is intentionally not accepted as a fragment manifest. */
	public static FragmentFamilySearch load(final String manifest, final String sha80, final int threads)throws Exception{
		final String pin=DigestSuffix.requireSuffix(sha80, "fragment manifest sha80");
		require(DigestSuffix.file(manifest).equals(pin), "Fragment manifest pin mismatch");
		final HashMap<String, String> policy=readPolicy(manifest);
		final FragmentFamilySearch result=FamilyTailCascade.loadFragments(policy, threads);
		require(DigestSuffix.file(manifest).equals(pin), "Fragment manifest changed during loading");
		return result;
	}
	/** Metadata-only entry point used to inspect/freeze a recipe before loading references. */
	public static HashMap<String, String> readPolicy(final String manifest){
		final ByteFile input=ByteFile.makeByteFile1(manifest, false);
		final LineParser1 lp=new LineParser1('\t');
		final HashMap<String, String> policy=new HashMap<String, String>();
		Throwable failure=null;
		try{
			final byte[] header=input.nextLine();
			require(header!=null && Arrays.equals(header, HEADER), "Expected a separate #bbgap_edge_fragment version1 package");
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				lp.set(line); require(lp.terms()==2, "Fragment policy rows require exactly key and value");
				final String key=lp.parseString(0), value=lp.parseString(1);
				require(!value.isEmpty() && policy.put(key, value)==null, "Empty or duplicate fragment policy row: "+key);
			}
			validatePolicy(policy);
		}catch(RuntimeException | Error e){failure=e; throw e;}
		finally{if(input.close()){
			final IllegalStateException e=new IllegalStateException("Failed closing fragment manifest");
			if(failure==null){throw e;}failure.addSuppressed(e);
		}}
		return policy;
	}
	static void validatePolicy(final HashMap<String, String> policy){
		require(policy!=null && policy.size()==FIXED.length/2+VARIABLE.length, "Fragment package has missing or unrecognized policy fields");
		for(int i=0; i<FIXED.length; i+=2){require(FIXED[i+1].equals(policy.get(FIXED[i])), "Fragment policy mismatch for "+FIXED[i]);}
		for(String key:VARIABLE){require(policy.get(key)!=null && !policy.get(key).isEmpty(), "Missing fragment policy field "+key);}
		final String id=policy.get("policyid");
		for(int i=0; i<id.length(); i++){require(id.charAt(i)>=32 && id.charAt(i)<127, "Fragment policyid must be printable ASCII");}
		DigestSuffix.requireSuffix(policy.get("familyreferencessha80"), "fragment reference sha80");
		DigestSuffix.requireSuffix(policy.get("ddlprofilessha80"), "fragment profiles sha80");
		new FragmentAlignment(Double.parseDouble(policy.get("identityfloor")), Double.parseDouble(policy.get("blosumfloor")));
	}
	public Worker newWorker(){return new Worker();}
	public int familyCount(){return config.families.length;}
	public String familyIdAt(final int ordinal){return config.families[ordinal].id;}
	public int familyRankAt(final int ordinal){return config.families[ordinal].rank;}
	public long familyMembersAt(final int ordinal){return config.families[ordinal].members;}

	/** Per-thread worker; observations/callback arrays are borrowed until the next query. */
	public final class Worker{
		public void observe(final String id, final byte[] raw, final Observation out, final Observer observer){
			require(!active && raw!=null && raw.length>=40 && raw.length<=100 && out!=null,
				"Fragment observation requires40..100 residues, output storage and a non-reentrant worker");
			active=true;
			try{
				final ProteinSequence query=new ProteinSequence(id, raw);
				require(query.length()==raw.length, "An edge fragment must not carry a terminal stop marker");
				for(byte aa:query.enc){require(Blosum62.isStandard(aa), "Fragment query contains an ambiguous residue");}
				out.clear();
				final boolean admitted=searcher.prepare(raw, query, result, counters, null, false);
				out.rawFamilies=searcher.work.size; out.admittedFamilies=result.ddlActive;
				out.windowFamilies=result.cosineWindow;
				if(observer!=null){observer.rawHits(searcher.work.touched, searcher.work.scores, searcher.work.size);}
				if(!admitted){out.reason=Reason.RAW_REJECT; return;}
				assert(result.dpShortlist==Math.min(50, searcher.shortlist.size())) :
					"Fragment policy disables all scalar filters and direct-accept shortcuts after ranking";
				for(int i=0; i<searcher.shortlist.size(); i++){
					final FamilyTailCascade.Comparison c=searcher.shortlist.get(i);
					if(observer!=null){observer.ranked(c.ordinal, c.ddlHits, c.cosine, i, i<result.dpShortlist);}
					if(i>=result.dpShortlist){continue;}
					final FragmentAlignment.Reason reason=scorer.measure(query, c.family.consensus, measured);
					out.alignedFamilies++;
					if(reason==FragmentAlignment.Reason.ACCEPT){
						out.acceptedFamilies++;
						if(out.familyOrdinal<0){out.familyOrdinal=c.ordinal;}
					}
					if(observer!=null){observer.aligned(c.ordinal, reason, measured);}
				}
				out.reason=out.familyOrdinal>=0 ? Reason.ACCEPT : Reason.ALIGNMENT_REJECT;
			}finally{active=false;}
		}
		private final FamilyTailCascade.Searcher searcher=new FamilyTailCascade.Searcher(config);
		private final FamilyTailCascade.Result result=new FamilyTailCascade.Result();
		private final FamilyTailCascade.Counters counters=new FamilyTailCascade.Counters();
		private final FragmentAlignment.Worker scorer=alignment.newWorker();
		private final FamilyAlignmentMetric.Observation measured=new FamilyAlignmentMetric.Observation();
		private boolean active;
	}
	/** Callbacks are synchronous and cannot mutate, retain or re-enter borrowed worker state.
	 * Raw ordinals + fixed floor and the ranked window expose every retrieval loss. */
	public interface Observer{
		default void rawHits(int[] touched, int[] scores, int size){}
		default void ranked(int ordinal, int hits, float cosine, int position, boolean aligned){}
		default void aligned(int ordinal, FragmentAlignment.Reason reason, FamilyAlignmentMetric.Observation measured){}
	}
	public static final class Observation{
		private void clear(){rawFamilies=admittedFamilies=windowFamilies=alignedFamilies=acceptedFamilies=0; familyOrdinal=-1; reason=Reason.RAW_REJECT;}
		public int rawFamilies, admittedFamilies, windowFamilies, alignedFamilies, acceptedFamilies, familyOrdinal=-1;
		public Reason reason=Reason.RAW_REJECT;
	}
	public enum Reason{RAW_REJECT, ALIGNMENT_REJECT, ACCEPT}
	public final String policyId, referenceSha80, profilesSha80;
	private final FamilyTailCascade.SearchConfig config;
	private final FragmentAlignment alignment;
	private static final String[] VARIABLE={"policyid", "identityfloor", "blosumfloor", "familyreferences", "familyreferencessha80", "ddlprofiles", "ddlprofilessha80"};
	private static final String[] FIXED={
		"querypolicy", ContigEdgeFragment.POLICY, "geneticcode", "explicit_per_input", "minqueryaa", "40", "maxqueryaa", "100",
		"alignmentmode", "glocal_full_query", "pairedcolumns", "28", "gaterule", "identity_and_blosum",
		"ranking", "raw_ddl_plus_monomer_cosine", "k", "6", "width", "512", "seed", "0", "floor", "5", "top", "50",
		"positionmode", "fixed", "partitions", "4", "overlap", "2", "aabits", "all:A/C/DE/FWY/G/IV/LM/NH/P/Q/RK/ST"};
	private static final byte[] HEADER="#bbgap_edge_fragment\t1".getBytes(java.nio.charset.StandardCharsets.US_ASCII);
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
}
