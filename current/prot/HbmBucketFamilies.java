package prot;

import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/**
 * Materializes one disjoint rank bucket from all completed parent extractions.
 * Retains source order and derived raw residues, compares every parent's family
 * counts and residue totals, and emits explicit EMPTY rows instead of losing
 * zero-member families. Input buckets are pinned products of HbmAssignedBuckets;
 * global assignment uniqueness/collection remains a separate publication gate.
 * @author Keqing
 */
public final class HbmBucketFamilies {

	public static void main(String[] args){
		try{
			Shared.setThreads(1); final PreParser pp=new PreParser(args, HbmBucketFamilies.class, false);
			try{run(HmmComparisonData.options(pp.args));}finally{Shared.closeStream(pp.outstream);}
		}catch(Throwable failure){failure.printStackTrace(); System.exit(1);}
	}

	private static void run(HashMap<String,String> o) throws Exception{
		for(String key : o.keySet()){
			require(Arrays.asList("identity", "identitysha80", "queriessha80", "root", "out", "parents", "buckets", "bucket").contains(key), "Unknown materialization parameter");
		}
		final String identity=HmmComparisonData.required(o, "identity"), identityPin=HmmComparisonData.required(o, "identitysha80");
		final String queryPin=HmmComparisonData.required(o, "queriessha80"); DigestSuffix.requireSuffix(queryPin, "query manifest pin");
		HbmProfilePilot.requireHash(identity, identityPin);
		final int parents=Integer.parseInt(HmmComparisonData.required(o, "parents")), buckets=Integer.parseInt(HmmComparisonData.required(o, "buckets"));
		final int bucket=Integer.parseInt(HmmComparisonData.required(o, "bucket"));
		require(parents>0 && buckets>0 && buckets<=64 && bucket>=0 && bucket<buckets, "Invalid parent/bucket dimensions");
		final List<String[]> rows=HmmComparisonData.rows(identity);
		require(rows.size()>1 && String.join("\t", rows.get(0)).equals("active_index\tfamily_id\tsource_active_index\trep_id\tmembers\traw_residues"), "Invalid identity schema");
		final int n=rows.size()-1; final int[] stable=new int[n]; final String[] reps=new String[n];
		final HashSet<Integer> ids=new HashSet<Integer>(); final HashSet<String> repIds=new HashSet<String>();
		for(int rank=0; rank<n; rank++){
			final String[] r=rows.get(rank+1); require(r.length==6 && r[0].equals(Integer.toString(rank)), "Noncontiguous family identity");
			stable[rank]=Integer.parseInt(r[1]); reps[rank]=r[3];
			require(stable[rank]>=0 && ids.add(stable[rank]) && repIds.add(reps[rank]) && HbmCompetitiveAssign.queryId(reps[rank]).equals(reps[rank]), "Invalid permanent identity or representative token");
		}
		final Path root=Paths.get(HmmComparisonData.required(o, "root")), out=Paths.get(HmmComparisonData.required(o, "out"));
		require(!Files.exists(out), "Materialization output must be fresh"); Files.createDirectory(out);
		final OutputStream[] writers=new OutputStream[n]; final long[][] counts=new long[n][3];
		try{
			for(int rank=bucket; rank<n; rank+=buckets){writers[rank]=ReadWrite.getOutputStream(out.resolve("rank_"+rank+".faa").toString(), false, true, false);}
			for(int parent=0; parent<parents; parent++){
				readParent(root.resolve("parent_"+parent), parent, bucket, buckets, queryPin, identityPin, stable, reps, writers, counts);
				System.err.println("HBM_FAMILY_MATERIALIZATION_PROGRESS bucket="+bucket+" parents="+(parent+1));
			}
		}finally{
			boolean failed=false; for(OutputStream writer : writers){if(writer!=null){failed|=ReadWrite.close(writer);}}
			require(!failed, "Family FASTA output I/O failure");
		}
		HbmProfilePilot.requireHash(identity, identityPin);
		final ByteStreamWriter manifest=HmmComparisonData.writer(out.resolve("family_manifest.tsv").toString());
		manifest.println("active_index\tfamily_id\trep_id\tassigned_count\tderived_raw_residues\tencoded_residues\tfasta_file\tfasta_sha80\tmembership_status");
		final ByteBuilder row=new ByteBuilder(); long total=0, raw=0, encoded=0; int empty=0, families=0;
		for(int rank=bucket; rank<n; rank+=buckets){
			final String name="rank_"+rank+".faa"; final long[] c=counts[rank];
			require(c[0]>0 || c[1]==0 && c[2]==0 && Files.size(out.resolve(name))==0, "Empty family has residue/output content");
			manifest.println(row.clear().append(rank).tab().append(stable[rank]).tab().append(reps[rank]).tab().append(c[0]).tab().append(c[1]).tab().append(c[2])
				.tab().append(name).tab().append(DigestSuffix.file(out.resolve(name).toString())).tab().append(c[0]==0 ? "EMPTY" : "NONEMPTY"));
			families++; if(c[0]==0){empty++;} total+=c[0]; raw+=c[1]; encoded+=c[2];
		}
		HmmComparisonData.close(manifest);
		final ByteStreamWriter summary=HmmComparisonData.writer(out.resolve("summary.tsv").toString());
		summary.println("bucket\tfamilies\tempty_families\tmembers\tderived_raw_residues\tencoded_residues");
		summary.println(row.clear().append(bucket).tab().append(families).tab().append(empty).tab().append(total).tab().append(raw).tab().append(encoded));
		HmmComparisonData.close(summary);
		final ByteStreamWriter recipe=HmmComparisonData.writer(out.resolve("recipe.tsv").toString());
		recipe.println("identity_sha80\t"+identityPin+"\nqueries_sha80\t"+queryPin+"\nparents\t"+parents+"\nbuckets\t"+buckets+"\nbucket\t"+bucket);
		HmmComparisonData.close(recipe);
		final ByteStreamWriter pass=HmmComparisonData.writer(out.resolve("PASS").toString());
		pass.println("HBM_BUCKET_FAMILIES_PASS global_collection_pending=true"); HmmComparisonData.close(pass);
	}

	/** Verifies and consumes one parent contribution, with one owner for each output family. */
	private static void readParent(Path parentDir, int parent, int bucket, int buckets, String queryPin, String identityPin,
			int[] stable, String[] reps, OutputStream[] writers, long[][] totals) throws Exception{
		require(Files.readAllLines(parentDir.resolve("PASS")).equals(java.util.Collections.singletonList("HBM_ASSIGNED_BUCKETS_PASS global_collection_pending=true")), "Parent extraction is incomplete");
		final String recipeFile=parentDir.resolve("recipe.tsv").toString(), recipePin=DigestSuffix.file(recipeFile);
		final HashMap<String,String> recipe=new HashMap<String,String>();
		for(String[] r : HmmComparisonData.rows(recipeFile)){require(r.length==2 && recipe.put(r[0], r[1])==null, "Malformed parent recipe");}
		require(identityPin.equals(recipe.get("identity_sha80")) && queryPin.equals(recipe.get("queries_sha80"))
			&& Integer.toString(parent).equals(recipe.get("parent")) && Integer.toString(buckets).equals(recipe.get("buckets")), "Parent identity/query/order/bucket binding differs");
		final String manifestFile=parentDir.resolve("buckets.tsv").toString(), manifestPin=DigestSuffix.file(manifestFile);
		final List<String[]> entries=HmmComparisonData.rows(manifestFile);
		require(entries.size()==buckets+1 && String.join("\t", entries.get(0)).equals("bucket\tfile\tsha80"), "Invalid parent bucket inventory");
		for(int b=0; b<buckets; b++){
			final String[] r=entries.get(b+1); require(r.length==3 && r[0].equals(Integer.toString(b)) && r[1].equals("bucket_"+b+".tsv"), "Parent bucket inventory order/path differs");
			DigestSuffix.requireSuffix(r[2], "parent bucket pin");
		}
		final String input=parentDir.resolve(entries.get(bucket+1)[1]).toString(), inputPin=entries.get(bucket+1)[2];
		HbmProfilePilot.requireHash(input, inputPin);
		final String countsFile=parentDir.resolve("family_counts.tsv").toString(), countsPin=DigestSuffix.file(countsFile);
		final List<String[]> expected=HmmComparisonData.rows(countsFile);
		require(expected.size()==stable.length+1 && String.join("\t", expected.get(0)).equals("active_index\tfamily_id\trep_id\tassigned\tderived_raw_residues\tencoded_residues"), "Invalid parent family counts");
		final long[][] actual=new long[stable.length][3]; final LineParser1 lp=new LineParser1('\t'); final ByteBuilder fasta=new ByteBuilder();
		final ByteFile bf=ByteFile.makeByteFile(input, false);
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lp.set(line); require(lp.terms()==5, "Bucket row requires five fields");
				final int rank=lp.parseInt(0); require(rank>=0 && rank<stable.length && rank%buckets==bucket, "Family occurs in the wrong bucket");
				require(lp.parseInt(1)==stable[rank] && lp.termEquals(reps[rank], 2), "Bucket permanent family identity differs");
				int idLength=lp.length(3); require(idLength>0, "Empty bucket query ID");
				for(int i=0; i<idLength; i++){final byte b=lp.parseByteFromCurrentField(i); require(b>32 && b<127, "Non-token bucket query ID");}
				final int raw=lp.length(4); require(raw>0, "Empty bucket sequence");
				// These are hash-bound, constructor-validated extractor bytes. Count the
				// conventional terminal marker without allocating a second protein object.
				final int encoded=raw-(lp.parseByteFromCurrentField(raw-1)=='*' ? 1 : 0);
				require(encoded>0, "Empty encoded bucket sequence");
				for(int i=0; i<encoded; i++){
					final byte b=lp.parseByteFromCurrentField(i);
					require((b>='A' && b<='Z') || (b>='a' && b<='z'), "Malformed derived residue or nonterminal stop in bucket");
				}
				fasta.clear().append('>'); lp.appendTerm(fasta, 3).nl(); lp.appendTerm(fasta, 4).nl();
				writers[rank].write(fasta.array, 0, fasta.length());
				actual[rank][0]++; actual[rank][1]+=raw; actual[rank][2]+=encoded;
			}
		}finally{require(!bf.close(), "Bucket reader failed");}
		for(int rank=0; rank<stable.length; rank++){
			final String[] r=expected.get(rank+1);
			require(r.length==6 && r[0].equals(Integer.toString(rank)) && Integer.parseInt(r[1])==stable[rank] && r[2].equals(reps[rank]), "Parent family count identity differs");
			for(int j=0; j<3; j++){
				final long count=Long.parseLong(r[j+3]); require(count>=0, "Negative parent family count");
				if(rank%buckets==bucket){require(actual[rank][j]==count, "Materialized parent family count/residue total differs"); totals[rank][j]=Math.addExact(totals[rank][j], count);}
			}
		}
		HbmProfilePilot.requireHash(input, inputPin); HbmProfilePilot.requireHash(manifestFile, manifestPin);
		HbmProfilePilot.requireHash(recipeFile, recipePin); HbmProfilePilot.requireHash(countsFile, countsPin);
	}

	private static void require(boolean ok, String reason){if(!ok){throw new IllegalArgumentException(reason);}}
}
