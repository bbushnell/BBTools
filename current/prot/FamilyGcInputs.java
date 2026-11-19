package prot;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.util.Arrays;

import ddl.DDLFamilySearch.StructuralInputs;
import fileIO.ByteFile;
import parse.LineParser1;

/** Loads family-GC reports against their exact collection maps. The current caller
 * schema covers SHORT and LONG sense tails; the older SHORT report remains readable.
 * Both streams retain every original family in physical-rank order. No row is
 * dropped, inferred or joined by rank alone; stable ID, label and count must agree.
 * Structural stability is observable metadata, not a fresh content fingerprint.
 * @author Keqing, Qiqi
 */
public final class FamilyGcInputs{

	/** Compact family resources fold map-validated sense-tail GC rows into the live roster. */
	public static FamilyGcPenalty fromCompact(FamilyReference[] families,int[] stable,long[] measured,
			double[] means,double[] sampleSd,double coefficient,boolean filter,double z){
		require(families!=null && stable!=null && measured!=null && means!=null && sampleSd!=null,
			"Compact GC binding requires families, stable IDs, measured counts, means and sample SD");
		final int n=families.length;
		require(n>0 && stable.length==n && measured.length==n && means.length==n && sampleSd.length==n,
			"Compact GC arrays must cover the same nonempty live roster");
		final int[] ranks=new int[n];
		final String[] ids=new String[n];
		final long[] members=new long[n];
		final boolean[] available=new boolean[n];
		for(int i=0; i<n; i++){
			final FamilyReference family=families[i];
			require(family!=null && family.members>0,"Compact GC can only bind live families");
			require(measured[i]==family.members,"Compact GC measured count differs from family members");
			ranks[i]=family.rank; ids[i]=family.id; members[i]=family.members; available[i]=true;
		}
		return new FamilyGcPenalty(stable,ranks,ids,members,available,means,sampleSd,coefficient,filter,z,30);
	}

	/** Current caller statistics: every assigned member, actual sense DNA, stop excluded.
	 * A coverage deficit is an input failure, never a silently reduced population. */
	public static FamilyGcPenalty loadTailStats(Path report,Path familyMap,int expectedFamilies,long expectedMembers,
		int minAa,int maxAa,double coefficient,boolean filter,double z,StructuralInputs tracked)throws IOException{
		require(report!=null && familyMap!=null && tracked!=null,"Tail GC requires a report, map and stability tracker");
		require(expectedFamilies>0 && expectedMembers>=expectedFamilies,"Tail GC requires declared full-roster counts");
		require(minAa==40 && maxAa==99 || minAa==100 && maxAa==100,"Tail GC is defined for SHORT40..99 or LONG100");
		final Path rp=tracked.track(report),mp=tracked.track(familyMap);
		final int[] stable=new int[expectedFamilies],ranks=new int[expectedFamilies];
		final String[] ids=new String[expectedFamilies];
		final long[] members=new long[expectedFamilies];
		final boolean[] available=new boolean[expectedFamilies];
		final double[] means=new double[expectedFamilies],sd=new double[expectedFamilies];
		final ByteFile rf=ByteFile.makeByteFile1(rp.toString(),false);
		ByteFile mf=null;
		final LineParser1 r=new LineParser1('\t'),m=new LineParser1('\t');
		long total=0;int rows=0;
		try{
			mf=ByteFile.makeByteFile1(mp.toString(),false);
			require(Arrays.equals(rf.nextLine(),TAIL_HEADER),"Unexpected sense-tail GC report header");
			require(Arrays.equals(mf.nextLine(),MAP_HEADER),"Unexpected family map header");
			for(byte[] row=rf.nextLine();row!=null;row=rf.nextLine()){
				require(rows<expectedFamilies,"Excess GC families");
				final byte[] mapping=mf.nextLine();require(mapping!=null,"Map ended before GC report");
				r.set(row);m.set(mapping);require(r.terms()==10 && m.terms()==5,"Tail GC/map require ten/five columns");
				stable[rows]=r.parseInt(0);ranks[rows]=r.parseInt(1);ids[rows]=r.parseString(2);members[rows]=r.parseLong(3);
				final long residues=r.parseLong(4),measured=r.parseLong(5);
				require(stable[rows]>=0 && ranks[rows]==rows && members[rows]>0,"Invalid live GC roster row");
				require(m.parseInt(0)==stable[rows] && m.parseInt(1)==ranks[rows] && m.termEquals(ids[rows],2) &&
					m.parseLong(3)==members[rows] && m.parseLong(4)==residues,
					"GC identity or counts differ from current collection map");
				require(measured==members[rows],"Missing assigned-member GC coverage at family "+ids[rows]);
				require(residues>=Math.multiplyExact((long)minAa,measured) && residues<=Math.multiplyExact((long)maxAa,measured),
					"GC member residues do not match the requested view");
				require(r.parseInt(8)==minAa && r.parseInt(9)==maxAa,"GC interval differs from requested sense-tail view");
				available[rows]=true;means[rows]=exactFraction(r,6);sd[rows]=exactFraction(r,7);
				require(measured>1 || sd[rows]==0,"Singleton SD must use zero storage");
				total=Math.addExact(total,measured);rows++;
			}
			require(mf.nextLine()==null,"Map contains families missing from GC report");
		}finally{
			final boolean failed=rf.close(),mapFailed=mf!=null && mf.close();
			require(!failed && !mapFailed,"Tail GC input read failed");
		}
		require(rows==expectedFamilies && total==expectedMembers,"Tail GC does not conserve declared population");
		tracked.verify();
		return new FamilyGcPenalty(stable,ranks,ids,members,available,means,sd,coefficient,filter,z,30);
	}

	/** The caller retains the tracker and verifies it before/after downstream work.
	 * This released report has complete mean coverage, including singleton means.
	 * Missing query GC remains a separate per-query flag in FamilyGcPenalty.score.
	 */
	public static FamilyGcPenalty load(Path report,Path familyMap,int expectedFamilies,long expectedMembers,
		double coefficient,StructuralInputs tracked)throws IOException{
		require(report!=null && familyMap!=null && tracked!=null,
			"GC loading requires the final report, collection map and shared input tracker");
		require(expectedFamilies>0 && expectedMembers>=expectedFamilies,
			"GC loading requires a nonempty independently counted full family cohort");
		require(Double.isFinite(coefficient) && coefficient>=0,
			"GC coefficient must be finite and nonnegative before loading metadata");
		final Path reportPath=tracked.track(report),mapPath=tracked.track(familyMap);
		final int[] stable=new int[expectedFamilies],ranks=new int[expectedFamilies];
		final String[] ids=new String[expectedFamilies];
		final long[] members=new long[expectedFamilies];
		final boolean[] available=new boolean[expectedFamilies];
		final double[] means=new double[expectedFamilies];
		final ByteFile rf=ByteFile.makeByteFile1(reportPath.toString(),false);
		ByteFile mf=null;
		final LineParser1 r=new LineParser1('\t'),m=new LineParser1('\t');
		long total=0;int rows=0;
		try{
			mf=ByteFile.makeByteFile1(mapPath.toString(),false);
			require(Arrays.equals(rf.nextLine(),REPORT_HEADER),"Unexpected final family-GC report header");
			require(Arrays.equals(mf.nextLine(),MAP_HEADER),"Unexpected collection family-map header");
			for(byte[] row=rf.nextLine();row!=null;row=rf.nextLine()){
				require(rows<expectedFamilies,"GC report contains more than the declared family population");
				final byte[] mapping=mf.nextLine();require(mapping!=null,"Collection map ended before the GC report");
				r.set(row);m.set(mapping);
				require(r.terms()==19 && m.terms()==5,"GC report/map rows require exactly19/5 fields");
				stable[rows]=r.parseInt(0);ranks[rows]=r.parseInt(1);ids[rows]=r.parseString(2);members[rows]=r.parseLong(3);
				require(stable[rows]>=0 && ranks[rows]==rows && members[rows]>0,
					"GC report requires nonnegative IDs, dense physical ranks and nonempty families");
				require(m.parseInt(0)==stable[rows] && m.parseInt(1)==ranks[rows] &&
					m.termEquals(ids[rows],2) && m.parseLong(3)==members[rows],
					"GC stable ID, physical rank, label or member count differs from the accepted collection map");
				final long residues=m.parseLong(4);
				require(residues>=Math.multiplyExact(40L,members[rows]) && residues<=Math.multiplyExact(99L,members[rows]),
					"GC collection must retain the released forty-through-ninety-nine residue SHORT population");
				require(r.termEquals("true",4),
					"The fixed complete GC report must supply every family mean; no silent missing-family substitution");
				available[rows]=true;means[rows]=exactFraction(r,5);
				if(members[rows]==1){
					require(r.termEquals("NA",6) && r.termEquals("false",7) && r.termEquals("singleton_neutral",8),
						"Singleton means are measured but sample SD and the family gate are unavailable");
				}else{
					final double sd=exactFraction(r,6);
					require(r.termEquals("true",7) && r.termEquals(sd==0 ? "zero_sd" : "available",8),
						"Measured family SD and gate availability/status must agree");
				}
				long prior=members[rows];
				for(int field=9;field<19;field++){
					final long excluded=r.parseLong(field);
					require(excluded>=0 && excluded<=prior,"GC exclusion counts must be bounded and nested across the fixed grid");
					require(!r.termEquals("singleton_neutral",8) && !r.termEquals("zero_sd",8) || excluded==0,
						"Singleton-neutral and constant-family retained members cannot be excluded");
					prior=excluded;
				}
				total=Math.addExact(total,members[rows]);rows++;
			}
			require(mf.nextLine()==null,"Collection map contains families absent from the GC report");
		}finally{
			final boolean reportFailed=rf.close(),mapFailed=mf!=null && mf.close();
			require(!reportFailed && !mapFailed,"GC report or collection map input failed");
		}
		require(rows==expectedFamilies && total==expectedMembers,
			"GC report does not conserve the declared family/member cohort");
		tracked.verify();
		return new FamilyGcPenalty(stable,ranks,ids,members,available,means,coefficient);
	}

	/** This small family metadata is loaded once. Java's decimal parser is deliberate:
	 * the native fast decimal parser does not promise exact Double.toString roundtrip.
	 * No process-global parser setting is changed and no query-row strings are retained.
	 */
	private static double exactFraction(LineParser1 parser,int field){
		final String text=parser.parseString(field);
		final double value=Double.parseDouble(text);
		require(Double.isFinite(value) && value>=0 && value<=1 && Double.toString(value).equals(text),
			"Family GC metadata must use the exact finite Double.toString fraction");
		return value;
	}
	private static void require(boolean valid,String message){if(!valid){throw new IllegalArgumentException(message);}}
	private FamilyGcInputs(){}
	private static final byte[] MAP_HEADER=("stable_family_id\troster_rank\tfamily_id\t"+
		"assigned\tresidues").getBytes(StandardCharsets.US_ASCII);
	public static final String TAIL_HEADER_TEXT="stable_family_id\troster_rank\tfamily_id\t"+
		"members\tresidues\tmeasured\tmean_gc\tsample_sd\ttail_min_aa\ttail_max_aa";
	private static final byte[] TAIL_HEADER=TAIL_HEADER_TEXT.getBytes(StandardCharsets.US_ASCII);
	private static final byte[] REPORT_HEADER=("stable_family_id\troster_rank\tfamily_id\t"+
		"members\tmean_available\tmean_gc\tsample_sd\tgate_available\tgate_status"+
		"\texcluded_N0\texcluded_N0.5\texcluded_N1\texcluded_N1.5\texcluded_N2"+
		"\texcluded_N2.5\texcluded_N3\texcluded_N4\texcluded_N5\texcluded_N6")
		.getBytes(StandardCharsets.US_ASCII);
}
