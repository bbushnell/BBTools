package prot;

import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Extracts one parent group's completed query shards into bounded family buckets.
 * Source and calls are consumed in lockstep, including rejected proteins. Output
 * keeps the caller's derived raw residues, including a normal single terminal stop;
 * derived-raw and encoded residue totals therefore remain distinct. No assignment
 * is recomputed. The separate collector still owns score/core and global-ID checks.
 * Intermediate PASS does not authorize final profile publication.
 * @author Keqing
 */
public final class HbmAssignedBuckets {

	public static void main(String[] args){
		try{
			Shared.setThreads(1); Shared.AMINO_IN=true; Read.VALIDATE_IN_CONSTRUCTOR=false;
			final PreParser pp=new PreParser(args, HbmAssignedBuckets.class, false);
			try{run(HmmComparisonData.options(pp.args));}finally{Shared.closeStream(pp.outstream);}
		}catch(Throwable failure){failure.printStackTrace(); System.exit(1);}
	}

	private static void run(HashMap<String,String> o) throws Exception{
		for(String key : o.keySet()){
			require(Arrays.asList("queries", "queriessha80", "identity", "identitysha80", "root", "out", "parent", "buckets", "expectedshards").contains(key), "Unknown bucket parameter: "+key);
		}
		final String queryFile=HmmComparisonData.required(o, "queries"), identityFile=HmmComparisonData.required(o, "identity");
		HbmProfilePilot.requireHash(queryFile, HmmComparisonData.required(o, "queriessha80"));
		HbmProfilePilot.requireHash(identityFile, HmmComparisonData.required(o, "identitysha80"));
		final List<String[]> queries=HmmComparisonData.rows(queryFile), identities=HmmComparisonData.rows(identityFile);
		require(queries.size()>1 && String.join("\t", queries.get(0)).equals(QUERY_HEADER), "Unexpected query manifest schema");
		require(identities.size()>1 && String.join("\t", identities.get(0)).equals(IDENTITY_HEADER), "Unexpected family identity schema");
		final int parent=Integer.parseInt(HmmComparisonData.required(o, "parent"));
		final int buckets=Integer.parseInt(HmmComparisonData.required(o, "buckets"));
		final int expectedShards=Integer.parseInt(HmmComparisonData.required(o, "expectedshards"));
		require(parent>=0 && buckets>=1 && buckets<=64 && expectedShards>0, "Invalid parent/bucket/shard dimensions");
		final int n=identities.size()-1;
		final int[] stable=new int[n]; final String[] reps=new String[n];
		final HashSet<Integer> seenIds=new HashSet<Integer>(); final HashSet<String> seenReps=new HashSet<String>();
		for(int rank=0; rank<n; rank++){
			final String[] r=identities.get(rank+1);
			require(r.length==6 && r[0].equals(Integer.toString(rank)), "Noncontiguous family identity");
			stable[rank]=Integer.parseInt(r[1]); reps[rank]=r[3];
			require(stable[rank]>=0 && seenIds.add(stable[rank]) && seenReps.add(reps[rank]) && HbmCompetitiveAssign.queryId(reps[rank]).equals(reps[rank]), "Invalid or duplicate permanent family identity");
		}
		final Path root=Paths.get(HmmComparisonData.required(o, "root")), out=Paths.get(HmmComparisonData.required(o, "out"));
		require(!Files.exists(out), "Bucket output must be a fresh directory"); Files.createDirectory(out);
		final OutputStream[] writers=new OutputStream[buckets]; final long[][] counts=new long[n][3];
		final long[] total=new long[7]; int shards=0, previousChild=-1;
		try{
			for(int b=0; b<buckets; b++){writers[b]=ReadWrite.getOutputStream(out.resolve("bucket_"+b+".tsv").toString(), false, true, false);}
			for(int i=1; i<queries.size(); i++){
				final String[] q=queries.get(i);
				require(q.length==8 && q[0].equals(Integer.toString(i-1)), "Noncontiguous query shard manifest");
				if(Integer.parseInt(q[1])!=parent){continue;}
				final int child=Integer.parseInt(q[2]); require(child==previousChild+1, "Parent child shards must be contiguous and ordered"); previousChild=child;
				final long[] local=extract(q, root, o.get("identitysha80"), stable, reps, writers, counts);
				for(int j=0; j<total.length; j++){total[j]=Math.addExact(total[j], local[j]);}
				shards++; System.err.println("HBM_BUCKET_PROGRESS parent="+parent+" shards="+shards+" queries="+total[0]);
			}
		}finally{
			boolean failed=false;
			for(OutputStream writer : writers){if(writer!=null){failed|=ReadWrite.close(writer);}}
			require(!failed, "Bucket output I/O failure");
		}
		require(shards==expectedShards, "Parent shard count differs: "+shards+" != "+expectedShards);
		HbmProfilePilot.requireHash(queryFile, o.get("queriessha80")); HbmProfilePilot.requireHash(identityFile, o.get("identitysha80"));
		final ByteBuilder row=new ByteBuilder();
		final ByteStreamWriter families=HmmComparisonData.writer(out.resolve("family_counts.tsv").toString());
		families.println("active_index\tfamily_id\trep_id\tassigned\tderived_raw_residues\tencoded_residues");
		long sum=0;
		for(int rank=0; rank<n; rank++){
			sum+=counts[rank][0];
			families.println(row.clear().append(rank).tab().append(stable[rank]).tab().append(reps[rank]).tab().append(counts[rank][0]).tab().append(counts[rank][1]).tab().append(counts[rank][2]));
		}
		HmmComparisonData.close(families); require(sum==total[1], "Bucket family totals do not conserve assigned queries");
		final ByteStreamWriter summary=HmmComparisonData.writer(out.resolve("summary.tsv").toString());
		summary.println("queries\tassigned\tsource_raw_residues\tall_encoded_residues\trepaired_records\tremoved_edge_markers\tassigned_derived_raw_residues");
		row.clear(); for(int i=0; i<total.length; i++){if(i>0){row.tab();}row.append(total[i]);} summary.println(row); HmmComparisonData.close(summary);
		final ByteStreamWriter manifest=HmmComparisonData.writer(out.resolve("buckets.tsv").toString());
		manifest.println("bucket\tfile\tsha80");
		for(int b=0; b<buckets; b++){
			final String name="bucket_"+b+".tsv";
			manifest.println(row.clear().append(b).tab().append(name).tab().append(DigestSuffix.file(out.resolve(name).toString())));
		}
		HmmComparisonData.close(manifest);
		final ByteStreamWriter recipe=HmmComparisonData.writer(out.resolve("recipe.tsv").toString());
		recipe.println("queries_sha80\t"+o.get("queriessha80")+"\nidentity_sha80\t"+o.get("identitysha80")+"\nparent\t"+parent+"\nbuckets\t"+buckets+"\nshards\t"+shards);
		HmmComparisonData.close(recipe);
		final ByteStreamWriter pass=HmmComparisonData.writer(out.resolve("PASS").toString());
		pass.println("HBM_ASSIGNED_BUCKETS_PASS global_collection_pending=true"); HmmComparisonData.close(pass);
	}

	/** Checks an accepted shard and routes each assigned source record exactly once. */
	private static long[] extract(String[] q, Path root, String identityPin, int[] stable, String[] reps,
			OutputStream[] writers, long[][] totalFamilies) throws Exception{
		final int shard=Integer.parseInt(q[0]); final Path dir=root.resolve("results/"+shard);
		exact(root.resolve("jobs/"+shard+".exit"), "0"); exact(dir.resolve("SHARD_PASS"), "COMPETITIVE_SHARD_PASS");
		final String input=q[5], callsFile=dir.resolve("assignments.tsv.gz").toString();
		HbmProfilePilot.requireHash(input, q[6]);
		final String manifestFile=dir.resolve("artifacts.tsv").toString(), manifestPin=DigestSuffix.file(manifestFile);
		final List<String[]> artifacts=HmmComparisonData.rows(manifestFile);
		require(artifacts.size()==4 && artifacts.get(1).length==3 && artifacts.get(1)[0].equals("assignments.tsv"), "Invalid shard artifact inventory");
		HbmProfilePilot.requireHash(callsFile, artifacts.get(1)[2]);
		final String recipeFile=dir.resolve("recipe.tsv").toString(), recipePin=DigestSuffix.file(recipeFile);
		final HashMap<String,String> recipe=new HashMap<String,String>();
		for(String[] r : HmmComparisonData.rows(recipeFile)){require(r.length==2 && recipe.put(r[0], r[1])==null, "Invalid shard recipe");}
		require(q[6].equals(recipe.get("in_sha80")) && identityPin.equals(recipe.get("identity_sha80")), "Shard source or identity binding differs");
		final String summaryFile=dir.resolve("summary.tsv").toString(), summaryPin=DigestSuffix.file(summaryFile);
		final List<String[]> summary=HmmComparisonData.rows(summaryFile);
		require(summary.size()==2 && summary.get(1).length==10, "Invalid assignment summary");
		final long[] s=new long[10]; for(int i=0; i<s.length; i++){s[i]=Long.parseLong(summary.get(1)[i]); require(s[i]>=0, "Negative summary counter");}
		require(s[0]==Long.parseLong(q[3]) && s[3]==Long.parseLong(q[4]), "Summary differs from original query manifest");
		final ByteFile calls=ByteFile.makeByteFile(callsFile, false); final LineParser1 lp=new LineParser1('\t');
		final Streamer source=StreamerFactory.makeStreamer(FileFormat.testInput(input, FileFormat.FASTA, null, true, false), null, true, -1);
		final ByteBuilder row=new ByteBuilder(); final long[] c=new long[7], families=new long[stable.length];
		source.start();
		try{
			require(Arrays.equals(calls.nextLine(), CALL_HEADER), "Unexpected assignment header");
			for(ListNum<Read> list=source.nextList(); list!=null && list.size()>0; list=source.nextList()){
				for(Read read : list){
					require(read.mate==null && read.bases!=null, "Expected unpaired source protein");
					final String id=HbmCompetitiveAssign.queryId(read.id);
					final byte[] normalized=HbmCompetitiveAssign.trimEdges(read.bases, id);
					final int encoded=new ProteinSequence(id, normalized).length();
					final byte[] call=calls.nextLine(); require(call!=null, "Assignments ended before source proteins"); lp.set(call);
					if(lp.terms()!=13 || !lp.termEquals(id, 0)){throw new IllegalArgumentException("Source/assignment query ID mismatch: "+id);}
					if(lp.parseInt(10)!=encoded){throw new IllegalArgumentException("Source/assignment encoded length mismatch: "+id);}
					c[0]++; c[2]+=read.bases.length; c[3]+=encoded;
					if(normalized!=read.bases){c[4]++; c[5]+=read.bases.length-normalized.length;}
					if(lp.termEquals("ASSIGNED", 1)){
						final int rank=lp.parseInt(2); require(rank>=0 && rank<stable.length, "Assigned family index outside roster");
						require(lp.parseInt(3)==stable[rank] && lp.termEquals(reps[rank], 4), "Assigned permanent family identity mismatch");
						row.clear().append(rank).tab().append(stable[rank]).tab().append(reps[rank]).tab().append(id).tab().append(normalized).nl();
						writers[rank%writers.length].write(row.array, 0, row.length());
						families[rank]++; totalFamilies[rank][0]++; totalFamilies[rank][1]+=normalized.length; totalFamilies[rank][2]+=encoded;
						c[1]++; c[6]+=normalized.length;
					}else{
						require(lp.termEquals("NO_ELIGIBLE_FAMILY", 1) && lp.parseInt(2)==-1 && lp.parseInt(3)==-1 && lp.parseInt(11)==-1 && lp.parseInt(12)==-1, "Invalid unassigned row");
						for(int j=4; j<=9; j++){require(lp.termEquals("NA", j), "Rejected row carries winner evidence");}
					}
				}
			}
			require(calls.nextLine()==null, "Assignments contain extra rows after source proteins");
		}finally{source.close(); require(!calls.close(), "Assignment reader failed");}
		require(!source.errorState(), "Source protein reader failed");
		require(c[0]==s[0] && c[1]==s[1] && c[0]-c[1]==s[2] && c[2]==s[3] && c[3]==s[4] && c[4]==s[5] && c[5]==s[6], "Extracted source/call counters disagree with summary");
		final String countsFile=dir.resolve("family_counts.tsv.gz").toString();
		require(artifacts.get(3).length==3 && artifacts.get(3)[0].equals("family_counts.tsv"), "Missing family count binding");
		HbmProfilePilot.requireHash(countsFile, artifacts.get(3)[2]);
		final List<String[]> countRows=HmmComparisonData.rows(countsFile);
		require(countRows.size()==stable.length+1, "Family count inventory differs");
		for(int rank=0; rank<stable.length; rank++){
			final String[] r=countRows.get(rank+1);
			require(r.length==4 && r[0].equals(Integer.toString(rank)) && Integer.parseInt(r[1])==stable[rank] && r[2].equals(reps[rank]) && Long.parseLong(r[3])==families[rank], "Extracted family counts differ at rank "+rank);
		}
		HbmProfilePilot.requireHash(input, q[6]); HbmProfilePilot.requireHash(callsFile, artifacts.get(1)[2]);
		HbmProfilePilot.requireHash(countsFile, artifacts.get(3)[2]); HbmProfilePilot.requireHash(manifestFile, manifestPin);
		HbmProfilePilot.requireHash(recipeFile, recipePin); HbmProfilePilot.requireHash(summaryFile, summaryPin);
		return c;
	}

	private static void exact(Path file, String expected) throws Exception{
		require(Files.readAllLines(file).equals(java.util.Collections.singletonList(expected)), "Incomplete assignment shard: "+file);
	}
	private static void require(boolean ok, String reason){if(!ok){throw new IllegalArgumentException(reason);}}
	private static final String QUERY_HEADER="shard\tparent\tchild\trecords\traw_residues\tfasta_file\tfasta_sha80\traw_sha80";
	private static final String IDENTITY_HEADER="active_index\tfamily_id\tsource_active_index\trep_id\tmembers\traw_residues";
	private static final byte[] CALL_HEADER="query\tstatus\tactive_index\tfamily_id\trep_id\tscore64\tidentity_pct\tcoverage_q\tcoverage_core\tmutual_coverage\tquery_length\tref_start\tref_end".getBytes(StandardCharsets.US_ASCII);
}
