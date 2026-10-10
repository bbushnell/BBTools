package prot;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/**
 * Independently checks completed experimental shards and totals their assignments.
 * Reads the hash-verified compressed calls, checks every row against family/core
 * metadata and per-shard counters, and writes all query IDs for a subsequent
 * global uniqueness check. Its PASS does not certify that separate uniqueness step.
 * @author Keqing
 */
public final class HbmAssignmentCollect {
	public static void main(String[] args){
		try{
			Shared.setThreads(1); PreParser pp=new PreParser(args, HbmAssignmentCollect.class, false);
			try{run(HmmComparisonData.options(pp.args));}finally{Shared.closeStream(pp.outstream);}
		}catch(Throwable error){error.printStackTrace(); System.exit(1);}
	}
	private static void run(HashMap<String,String> o) throws Exception{
		for(String key : o.keySet()){require(Arrays.asList("queries", "queriessha80", "modelconfig", "modelconfigsha80", "root", "out", "expected").contains(key), "Unknown collection parameter");}
		final String queryFile=HmmComparisonData.required(o, "queries"), configFile=HmmComparisonData.required(o, "modelconfig");
		HbmProfilePilot.requireHash(queryFile, HmmComparisonData.required(o, "queriessha80"));
		HbmProfilePilot.requireHash(configFile, HmmComparisonData.required(o, "modelconfigsha80"));
		final HashMap<String,String> model=readConfig(configFile);
		for(String key : MODEL_KEYS){HbmProfilePilot.requireHash(HmmComparisonData.required(model, key), HmmComparisonData.required(model, key+"sha80"));}
		final List<String[]> queries=HmmComparisonData.rows(queryFile), identities=HmmComparisonData.rows(model.get("identity")), cores=HmmComparisonData.rows(model.get("cores"));
		require(queries.size()>1 && String.join("\t", queries.get(0)).equals("shard\tparent\tchild\trecords\traw_residues\tfasta_file\tfasta_sha80\traw_sha80"), "Invalid query shard manifest");
		require(identities.size()>1 && String.join("\t", identities.get(0)).equals("active_index\tfamily_id\tsource_active_index\trep_id\tmembers\traw_residues"), "Invalid identity table");
		require(cores.size()==identities.size(), "Core/identity table sizes differ");
		final int n=identities.size()-1; final String[] reps=new String[n]; final int[] stableIds=new int[n], lengths=new int[n], first=new int[n], last=new int[n];
		for(int i=0; i<n; i++){
			final String[] r=identities.get(i+1), c=cores.get(i+1);
			require(r.length==6 && c.length==9 && r[0].equals(Integer.toString(i)) && c[0].equals(r[0]) && c[1].equals(r[1]) && c[2].equals(r[3]), "Unbound family/core identity");
			reps[i]=r[3]; stableIds[i]=Integer.parseInt(r[1]); lengths[i]=Integer.parseInt(c[4]); require(lengths[i]>0, "Empty family consensus");
			first[i]=Integer.parseInt(c[6]); last[i]=Integer.parseInt(c[7]);
			require(first[i]>=0 && last[i]>=first[i] && last[i]<lengths[i] && Integer.parseInt(c[8])==last[i]-first[i]+1, "Invalid family core interval");
		}
		final long expected=Long.parseLong(HmmComparisonData.required(o, "expected")); require(expected>0, "Collection expected count must be positive");
		final Path root=Paths.get(HmmComparisonData.required(o, "root")), out=Paths.get(HmmComparisonData.required(o, "out"));
		require(!Files.exists(out), "Collection output must be a fresh directory"); Files.createDirectory(out);
		final ByteStreamWriter ids=HmmComparisonData.writer(out.resolve("query_ids.txt").toString()); final ByteBuilder idRow=new ByteBuilder();
		final long[] total=new long[10], families=new long[n], local=new long[n]; final LineParser1 lp=new LineParser1('\t');
		try{
			for(int shard=0; shard<queries.size()-1; shard++){
				final String[] q=queries.get(shard+1); require(q.length==8 && q[0].equals(Integer.toString(shard)), "Noncontiguous query shards");
				final long count=Long.parseLong(q[3]), raw=Long.parseLong(q[4]); require(count>0 && raw>0, "Empty query shard");
				final Path dir=root.resolve("results/"+shard);
				exact(root.resolve("jobs/"+shard+".exit"), "0"); exact(dir.resolve("SHARD_PASS"), "COMPETITIVE_SHARD_PASS"); exact(dir.resolve("PASS"), "HBM_COMPETITIVE_ASSIGN_PASS");
				final List<String[]> artifacts=HmmComparisonData.rows(dir.resolve("artifacts.tsv").toString());
				require(artifacts.size()==4 && String.join("\t", artifacts.get(0)).equals("artifact\traw_sha80\tcompressed_sha80"), "Invalid shard artifact manifest");
				for(int a=0; a<ARTIFACTS.length; a++){
					final String[] r=artifacts.get(a+1); require(r.length==3 && r[0].equals(ARTIFACTS[a]), "Unexpected shard artifact order");
					DigestSuffix.requireSuffix(r[1], "raw output pin"); HbmProfilePilot.requireHash(dir.resolve(r[0]+".gz").toString(), r[2]);
				}
				final HashMap<String,String> recipe=readPairs(dir.resolve("recipe.tsv").toString());
				for(String key : MODEL_KEYS){require(model.get(key+"sha80").equals(recipe.get(key+"_sha80")), "Shard uses a different model dependency");}
				require(q[6].equals(recipe.get("in_sha80")), "Shard query pin differs from manifest");
				for(int i=0; i<SETTINGS.length; i+=2){require(SETTINGS[i+1].equals(recipe.get(SETTINGS[i])), "Shard uses a different scoring recipe");}
				final List<String[]> summary=HmmComparisonData.rows(dir.resolve("summary.tsv").toString());
				require(summary.size()==2 && String.join("\t", summary.get(0)).equals(SUMMARY_HEADER) && summary.get(1).length==10, "Invalid shard summary");
				final long[] s=new long[10]; for(int i=0; i<s.length; i++){s[i]=Long.parseLong(summary.get(1)[i]); require(s[i]>=0, "Negative shard summary count");}
				require(s[0]==count && s[1]+s[2]==count && s[3]==raw && s[4]<=raw && s[5]<=s[6] && s[6]<=raw-s[4] && s[7]+s[9]==50*count && s[8]<=s[7], "Shard summary conservation failed");
				Arrays.fill(local, 0); long rows=0, accepted=0, encoded=0;
				final ByteFile calls=ByteFile.makeByteFile(dir.resolve("assignments.tsv.gz").toString(), false);
				try{
					require(Arrays.equals(calls.nextLine(), CALL_HEADER), "Unexpected assignment header");
					for(byte[] line=calls.nextLine(); line!=null; line=calls.nextLine()){
						lp.set(line); require(lp.terms()==13, "Assignment row needs13 fields");
						int end=0; while(end<line.length && line[end]!='\t'){require(line[end]>32 && line[end]<127, "Non-token query ID in assignment"); end++;}
						require(end>0 && end<line.length, "Missing query ID"); ids.print(idRow.clear().append(line, 0, end).nl());
						final int length=lp.parseInt(10); require(length>0, "Empty assigned/rejected query"); encoded+=length; rows++;
						if(lp.termEquals("ASSIGNED", 1)){
							final int rank=lp.parseInt(2); require(rank>=0 && rank<n, "Unknown assigned family index");
							require(lp.parseInt(3)==stableIds[rank] && lp.termEquals(reps[rank], 4), "Assigned family identity differs");
							lp.parseInt(5); final float identity=lp.parseFloat(6), cq=lp.parseFloat(7), ct=lp.parseFloat(8), mutual=lp.parseFloat(9);
							require(Float.isFinite(identity) && identity>=HbmCompetitiveGate.MIN_IDENTITY && identity<=100, "Assigned row violates identity gate");
							final int start=lp.parseInt(11), stop=lp.parseInt(12); require(start>=0 && stop>=start && stop<lengths[rank], "Assignment outside consensus frame");
							checkCoverage(length, lengths[rank], first[rank], last[rank], start, stop, cq, ct, mutual);
							accepted++; local[rank]++;
						}else{
							require(lp.termEquals("NO_ELIGIBLE_FAMILY", 1) && lp.parseInt(2)==-1 && lp.parseInt(3)==-1 && lp.parseInt(11)==-1 && lp.parseInt(12)==-1, "Invalid rejection row");
							for(int i=4; i<=9; i++){require(lp.termEquals("NA", i), "Rejected row retains winner evidence");}
						}
					}
				}finally{require(!calls.close(), "Assignment input failed");}
				require(rows==count && accepted==s[1] && encoded==s[4], "Recounted assignment rows disagree with summary");
				final ByteFile counts=ByteFile.makeByteFile(dir.resolve("family_counts.tsv.gz").toString(), false); int rank=0;
				try{
					require(Arrays.equals(counts.nextLine(), FAMILY_HEADER), "Invalid family count header");
					for(byte[] line=counts.nextLine(); line!=null; line=counts.nextLine()){
						lp.set(line); require(rank<n && lp.terms()==4 && lp.parseInt(0)==rank && lp.parseInt(1)==stableIds[rank] && lp.termEquals(reps[rank], 2) && lp.parseLong(3)==local[rank], "Recounted family assignments differ from producer table");
						families[rank]+=local[rank]; rank++;
					}
				}finally{require(!counts.close(), "Family count input failed");}
				require(rank==n, "Missing family counts");
				final ByteFile repairs=ByteFile.makeByteFile(dir.resolve("edge_repairs.tsv.gz").toString(), false);
				long repaired=0, removed=0;
				final java.util.HashSet<String> repairedIds=new java.util.HashSet<String>();
				try{
					require(Arrays.equals(repairs.nextLine(), REPAIR_HEADER), "Invalid edge-repair header");
					for(byte[] line=repairs.nextLine(); line!=null; line=repairs.nextLine()){
						lp.set(line); require(lp.terms()==4, "Edge-repair row needs4 fields");
						final String id=lp.parseString(0); final long before=lp.parseLong(1), after=lp.parseLong(2), delta=lp.parseLong(3);
						require(!id.isEmpty() && repairedIds.add(id) && before>after && after>0 && delta==before-after, "Invalid or duplicate edge repair");
						repaired++; removed+=delta;
					}
				}finally{require(!repairs.close(), "Repair input failed");}
				require(repaired==s[5] && removed==s[6], "Recounted edge repairs differ from summary");
				for(int i=0; i<10; i++){total[i]=Math.addExact(total[i], s[i]);}
				System.err.println("HBM_COLLECTION_PROGRESS shards="+(shard+1)+" queries="+total[0]);
			}
		}finally{HmmComparisonData.close(ids);}
		require(total[0]==expected, "Collected corpus count differs from expected");
		HbmProfilePilot.requireHash(queryFile, o.get("queriessha80")); HbmProfilePilot.requireHash(configFile, o.get("modelconfigsha80"));
		final ByteStreamWriter counts=HmmComparisonData.writer(out.resolve("family_counts.tsv").toString()); counts.println("active_index\tfamily_id\trep_id\tassigned");
		final ByteBuilder row=new ByteBuilder(); long sum=0;
		for(int i=0; i<n; i++){sum+=families[i]; counts.println(row.clear().append(i).tab().append(stableIds[i]).tab().append(reps[i]).tab().append(families[i]));}
		HmmComparisonData.close(counts); require(sum==total[1], "Collected families do not conserve assigned queries");
		final ByteStreamWriter summary=HmmComparisonData.writer(out.resolve("summary.tsv").toString()); summary.println(SUMMARY_HEADER);
		row.clear(); for(int i=0; i<10; i++){if(i>0){row.tab();}row.append(total[i]);} summary.println(row); HmmComparisonData.close(summary);
		final ByteStreamWriter pass=HmmComparisonData.writer(out.resolve("PASS").toString()); pass.println("HBM_ASSIGNMENT_COLLECTION_PASS global_uniqueness_pending=true"); HmmComparisonData.close(pass);
	}
	/**
	 * Checks existence of a paired-column count consistent with all serialized ratios.
	 * The positional DP's supported dimensions bound accepted lengths tightly enough
	 * that rounding coverageQ*queryLength recovers its integer paired count exactly.
	 * Paired columns must also fit the intersection of the core and reported span.
	 */
	static void checkCoverage(int queryLength, int referenceLength, int coreStart, int coreEnd,
			int start, int end, float cq, float ct, float mutual){
		require(queryLength>0 && referenceLength>0 && coreStart>=0 && coreEnd>=coreStart && coreEnd<referenceLength && start>=0 && end>=start && end<referenceLength,
			"Invalid coverage geometry coordinates");
		require(((long)queryLength+1)*(referenceLength+1L)<=Integer.MAX_VALUE-8 && ((long)queryLength+referenceLength)*6400<=Integer.MAX_VALUE/2,
			"Assigned dimensions exceed HbmPositionModel.align bounds");
		require(cq>=HbmCompetitiveGate.MIN_COVERAGE && cq<=1 && ct>=HbmCompetitiveGate.MIN_COVERAGE && ct<=1 && mutual==Math.min(cq, ct),
			"Assigned row violates paired-core coverage gate");
		final int coreLength=coreEnd-coreStart+1, intersection=Math.max(0, Math.min(end, coreEnd)-Math.max(start, coreStart)+1);
		final long paired=Math.round((double)cq*queryLength);
		require(paired>0 && paired<=queryLength && paired<=coreLength && paired<=intersection, "Paired coverage cannot fit inside the reported core/span intersection");
		require((float)(paired/(double)queryLength)==cq && (float)(paired/(double)coreLength)==ct && (float)(paired/(double)Math.max(queryLength, coreLength))==mutual,
			"Coverage floats disagree with an integer paired-column count");
	}
	private static HashMap<String,String> readConfig(String file){
		final java.util.ArrayList<String> lines=new java.util.ArrayList<String>(); final ByteFile bf=ByteFile.makeByteFile(file, false);
		try{for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){if(line.length>0 && line[0]!='#'){lines.add(new String(line, java.nio.charset.StandardCharsets.UTF_8));}}}
		finally{require(!bf.close(), "Model config read failed");}
		return HmmComparisonData.options(lines.toArray(new String[0]));
	}
	private static HashMap<String,String> readPairs(String file){
		final HashMap<String,String> out=new HashMap<String,String>();
		for(String[] r : HmmComparisonData.rows(file)){require(r.length==2 && out.put(r[0], r[1])==null, "Invalid or duplicate recipe key");} return out;
	}
	private static void exact(Path file, String expected) throws Exception{require(Files.readAllLines(file).equals(java.util.Collections.singletonList(expected)), "Missing successful shard completion marker: "+file);}
	private static void require(boolean ok, String message){if(!ok){throw new IllegalArgumentException(message);}}
	private static final String[] MODEL_KEYS={"consensus", "identity", "hbm", "provenance", "cores", "background", "sidecar", "kmers"};
	private static final String[] ARTIFACTS={"assignments.tsv", "edge_repairs.tsv", "family_counts.tsv"};
	private static final String[] SETTINGS={"profile","logodds","beta","0.01","clip_min","-4","clip_max","11","gap","4","shortlist","50","ranking","raw_score64_then_ascii_rep_id","identity_floor_float32","40.612846","paired_core_floor_float32","0.8"};
	private static final String SUMMARY_HEADER="queries\tassigned\tunassigned\traw_residues\tencoded_residues\trepaired_records\tremoved_edge_markers\tscored_candidates\trecorded_candidates\tlength_rejected";
	private static final byte[] CALL_HEADER="query\tstatus\tactive_index\tfamily_id\trep_id\tscore64\tidentity_pct\tcoverage_q\tcoverage_core\tmutual_coverage\tquery_length\tref_start\tref_end".getBytes(java.nio.charset.StandardCharsets.US_ASCII);
	private static final byte[] FAMILY_HEADER="active_index\tfamily_id\trep_id\tassigned".getBytes(java.nio.charset.StandardCharsets.US_ASCII);
	private static final byte[] REPAIR_HEADER="query\traw_length\tderived_length\tremoved_edge_markers".getBytes(java.nio.charset.StandardCharsets.US_ASCII);
}
