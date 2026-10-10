package prot;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Single-worker experimental full-shortlist positional assignment of one query shard.
 * Every input receives one output row, including explicit no-eligible-family rows.
 * Models, consensus, core coordinates and shortlist inputs must match supplied pins.
 * Parallel corpus execution uses independent shards rather than shared mutable state.
 * @author Keqing
 */
public final class HbmCompetitiveAssign {
	public static void main(String[] args){
		try{
			Shared.setThreads(1); Shared.AMINO_IN=true; Read.VALIDATE_IN_CONSTRUCTOR=false;
			final PreParser pp=new PreParser(args, HbmCompetitiveAssign.class, false);
			try{run(HmmComparisonData.options(pp.args));}finally{Shared.closeStream(pp.outstream);}
		}catch(Throwable failure){failure.printStackTrace(); System.exit(1);}
	}

	private static void run(HashMap<String,String> o) throws Exception{
		final HashSet<String> allowed=new HashSet<String>(Arrays.asList("out", "expected", "familylimit"));
		for(String key : PINNED){allowed.add(key); allowed.add(key+"sha80");}
		for(String key : o.keySet()){require(allowed.contains(key), "Unknown competitive parameter: "+key);}
		for(String key : PINNED){HbmProfilePilot.requireHash(HmmComparisonData.required(o, key), HmmComparisonData.required(o, key+"sha80"));}
		final long expected=Long.parseLong(HmmComparisonData.required(o, "expected"));
		require(expected>0, "Expected query count must be positive");
		final Path out=Paths.get(HmmComparisonData.required(o, "out"));
		require(!Files.exists(out), "Competitive output must be a fresh directory");
		final Library library=new Library(o);
		final int requestedLimit=parse.Parse.parseIntKMG(o.getOrDefault("familylimit", "0"));
		final int familyLimit=requestedLimit==0 ? library.refs.length : requestedLimit;
		require(familyLimit>=50 && familyLimit<=library.refs.length,
			"familylimit must be zero (all families) or between50 and the loaded roster size for a full top50");
		Files.createDirectory(out);
		final ByteStreamWriter calls=HmmComparisonData.writer(out.resolve("assignments.tsv").toString());
		final ByteStreamWriter repairs=HmmComparisonData.writer(out.resolve("edge_repairs.tsv").toString());
		calls.println("query\tstatus\tactive_index\tfamily_id\trep_id\tscore64\tidentity_pct\tcoverage_q\tcoverage_core\tmutual_coverage\tquery_length\tref_start\tref_end");
		repairs.println("query\traw_length\tderived_length\tremoved_edge_markers");
		final Streamer input=StreamerFactory.makeStreamer(FileFormat.testInput(o.get("in"), FileFormat.FASTA, null, true, false), null, true, -1);
		final ProteinSearcher.ShortlistScratch shortlist=new ProteinSearcher.ShortlistScratch(library.sidecar.dims, library.refs.length, 50);
		final HbmCompetitiveGate.Scratch selection=new HbmCompetitiveGate.Scratch(library.refs.length, 50);
		final ByteBuilder row=new ByteBuilder(); final long[] familyCounts=new long[library.refs.length];
		long count=0, assigned=0, raw=0, encoded=0, repaired=0, removed=0, scored=0, recorded=0, lengthRejected=0;
		input.start();
		try{
			for(ListNum<Read> list=input.nextList(); list!=null && list.size()>0; list=input.nextList()){
				for(Read read : list){
					require(read.mate==null && read.bases!=null && read.bases.length>0, "Expected a nonempty unpaired query protein");
					final String id=queryId(read.id); final byte[] normalized=trimEdges(read.bases, id);
					if(normalized!=read.bases){
						repaired++; removed+=read.bases.length-normalized.length;
						repairs.print(row.clear().append(id).tab().append(read.bases.length).tab().append(normalized.length).tab().append(read.bases.length-normalized.length).nl());
					}
					final ProteinSequence query=new ProteinSequence(id, normalized);
					ProteinSearcher.scoreFamiliesF4(query, library.sidecar, shortlist);
					top(shortlist.f4, shortlist.topIdx, familyLimit);
					final int rank=library.gate.select(query.enc, shortlist.topIdx, 50, selection);
					row.clear().append(id).tab();
					if(rank>=0){
						assigned++; familyCounts[rank]++;
						row.append("ASSIGNED").tab().append(rank).tab().append(library.familyIds[rank]).tab().append(library.ids[rank]).tab().append(selection.score64)
							.tab().append(selection.metrics.identity, 8).tab().append(selection.metrics.coverageQ, 8)
							.tab().append(selection.metrics.coverageT, 8).tab().append(selection.metrics.coverage, 8)
							.tab().append(query.length()).tab().append(selection.start).tab().append(selection.end);
					}else{row.append("NO_ELIGIBLE_FAMILY\t-1\t-1\tNA\tNA\tNA\tNA\tNA\tNA").tab().append(query.length()).append("\t-1\t-1");}
					calls.print(row.nl());
					count++; raw+=read.bases.length; encoded+=query.length();
					scored+=selection.scored; recorded+=selection.recorded; lengthRejected+=selection.lengthRejected;
					if(count%10000==0){System.err.println("HBM_COMPETITIVE_PROGRESS queries="+count+" assigned="+assigned);}
				}
			}
		}finally{input.close(); HmmComparisonData.close(calls); HmmComparisonData.close(repairs);}
		require(!input.errorState(), "Competitive query reader failed");
		require(count==expected, "Competitive input count differs: expected="+expected+" actual="+count);
		require(scored+lengthRejected==50*count && recorded<=scored, "Every query must score or length-reject all50 shortlisted families");
		for(String key : PINNED){HbmProfilePilot.requireHash(o.get(key), o.get(key+"sha80"));}
		final ByteStreamWriter counts=HmmComparisonData.writer(out.resolve("family_counts.tsv").toString());
		counts.println("active_index\tfamily_id\trep_id\tassigned"); long sum=0;
		for(int i=0; i<familyCounts.length; i++){sum+=familyCounts[i]; counts.println(row.clear().append(i).tab().append(library.familyIds[i]).tab().append(library.ids[i]).tab().append(familyCounts[i]));}
		HmmComparisonData.close(counts); require(sum==assigned && assigned<=count, "Competitive family counts do not conserve assigned queries");
		final ByteStreamWriter summary=HmmComparisonData.writer(out.resolve("summary.tsv").toString());
		summary.println("queries\tassigned\tunassigned\traw_residues\tencoded_residues\trepaired_records\tremoved_edge_markers\tscored_candidates\trecorded_candidates\tlength_rejected");
		summary.println(row.clear().append(count).tab().append(assigned).tab().append(count-assigned).tab().append(raw).tab().append(encoded).tab().append(repaired).tab().append(removed).tab().append(scored).tab().append(recorded).tab().append(lengthRejected));
		HmmComparisonData.close(summary);
		final ByteStreamWriter recipe=HmmComparisonData.writer(out.resolve("recipe.tsv").toString());
		recipe.println("profile\tlogodds\nbeta\t0.01\nclip_min\t-4\nclip_max\t11\ngap\t4\nshortlist\t50\nranking\traw_score64_then_ascii_rep_id\nidentity_floor_float32\t40.612846\npaired_core_floor_float32\t0.8");
		recipe.println("candidate_families\t"+familyLimit);
		for(String key : PINNED){recipe.println(key+"_sha80\t"+o.get(key+"sha80"));}
		HmmComparisonData.close(recipe);
		final ByteStreamWriter pass=HmmComparisonData.writer(out.resolve("PASS").toString()); pass.println("HBM_COMPETITIVE_ASSIGN_PASS"); HmmComparisonData.close(pass);
		System.err.println("HBM_COMPETITIVE_ASSIGN_PASS queries="+count+" assigned="+assigned);
	}

	/** Same first-whitespace token contract as ProteinSearch.readFasta; retain original source file. */
	static String queryId(String header){
		require(header!=null && !header.isEmpty(), "Empty query header"); int end=0;
		while(end<header.length() && !Character.isWhitespace(header.charAt(end))){
			final char c=header.charAt(end); require(c>32 && c<127, "Query identifiers must be printable ASCII tokens"); end++;
		}
		require(end>0, "Empty query ID before description"); return end==header.length() ? header : header.substring(0, end);
	}

	/** Repairs only malformed edge stops; strict ProteinSequence encoding rejects internal/unsupported residues. */
	static byte[] trimEdges(byte[] raw, String id){
		require(raw!=null && raw.length>0, "Empty query residues: "+id); int first=0, end=raw.length;
		while(first<end && raw[first]=='*'){first++;}
		while(end>first && raw[end-1]=='*'){end--;}
		require(first<end, "Empty query after edge stops: "+id);
		for(int i=first; i<end; i++){if(raw[i]=='*'){throw new IllegalArgumentException("Internal query stop in "+id+" at "+i);}}
		return first==0 && raw.length-end<=1 ? raw : Arrays.copyOfRange(raw, first, end);
	}

	/** Scalar top50: score descending, dense index ascending; nonfinite F4 scores fail rather than silently pad. */
	static void top(float[] scores, int[] out){
		top(scores, out, scores==null ? 0 : scores.length);
	}

	/** Selects within a roster prefix before top50; the full sidecar/F4 calculation stays unchanged. */
	static void top(float[] scores, int[] out, int familyLimit){
		require(scores!=null && out!=null && out.length>0 && out.length<=familyLimit && familyLimit<=scores.length,
			"Invalid shortlist dimensions: the candidate prefix must contain the whole requested shortlist"); int size=0;
		for(int i=0; i<familyLimit; i++){
			if(!Float.isFinite(scores[i])){throw new IllegalArgumentException("Nonfinite F4 score at family "+i);}
			if(size==out.length && scores[i]<=scores[out[size-1]]){continue;}
			int at=Math.min(size, out.length-1);
			while(at>0 && scores[i]>scores[out[at-1]]){at--;}
			if(size<out.length){size++;}
			for(int j=size-1; j>at; j--){out[j]=out[j-1];} out[at]=i;
		}
	}

	/** Loads all dependency bindings before exposing a selectable library. */
	private static final class Library{
		Library(HashMap<String,String> o) throws Exception{
			final byte[][] provenance=HbmBundleLoader.loadSemanticProvenance(o.get("provenance"));
			// HbmFamilyUnion.writeProvenance binds its identity table as cluster_membership.
			require(o.get("identitysha80").equals(DigestSuffix.fromDigest(provenance[13])), "Union identity table differs from HBM cluster-membership provenance");
			final List<ProteinSequence> sequences=ProteinSearch.readFasta(o.get("consensus")); final int n=sequences.size();
			require(n>=50, "Competitive experiment requires at least50 families"); refs=new byte[n][]; ids=new String[n]; familyIds=new int[n];
			final ArrayList<String> roster=new ArrayList<String>(); final HashMap<String,byte[]> byId=new HashMap<String,byte[]>();
			for(int i=0; i<n; i++){final ProteinSequence p=sequences.get(i); refs[i]=p.enc; ids[i]=p.id; roster.add(p.id); require(byId.put(p.id, p.enc)==null, "Duplicate consensus ID");}
			sidecar=FamilyShortlistSidecar.load(o.get("sidecar"), o.get("consensus"), o.get("kmers"));
			require(sidecar.schemaVersion==1 && Arrays.equals(ids, sidecar.repIds), "Competitive shortlist must have the same roster and no production thresholds");
			final List<String[]> identities=HmmComparisonData.rows(o.get("identity"));
			require(identities.size()==n+1 && String.join("\t", identities.get(0)).equals("active_index\tfamily_id\tsource_active_index\trep_id\tmembers\traw_residues"), "Unsupported union identity table");
			final List<String[]> cores=HmmComparisonData.rows(o.get("cores")); final int[] first=new int[n], last=new int[n];
			require(cores.size()==n+1 && String.join("\t", cores.get(0)).equals("active_index\tfamily_id\trep_id\tmembers\tconsensus_length\tcoverage_threshold\tcore_start\tcore_end\tcore_length"), "Unsupported core table");
			final HashSet<Integer> stable=new HashSet<Integer>();
			for(int i=0; i<n; i++){
				final String[] a=identities.get(i+1), c=cores.get(i+1);
				require(a.length==6 && c.length==9 && a[0].equals(Integer.toString(i)) && c[0].equals(a[0]) && a[1].equals(c[1]) && a[3].equals(ids[i]) && c[2].equals(ids[i]) && a[4].equals(c[3]), "Core/union identity mismatch at "+i);
				familyIds[i]=Integer.parseInt(a[1]); final long members=Long.parseLong(a[4]);
				require(familyIds[i]>=0 && stable.add(familyIds[i]) && members>0 && Long.parseLong(c[5])==(members+1)/2 && Integer.parseInt(c[4])==refs[i].length, "Invalid family ID/count/core frame at "+i);
				first[i]=Integer.parseInt(c[6]); last[i]=Integer.parseInt(c[7]); require(Integer.parseInt(c[8])==last[i]-first[i]+1, "Core width differs at "+i);
			}
			final HashSet<String> requiredHeaders=new HashSet<String>(Arrays.asList(
				"#profile\tlogodds\tbeta\t0.01\tclipmin\t-4\tclipmax\t11\tgap\t4",
				"#coverage\tpaired_m_only\tthreshold\tceil(member_count/2)\tcoordinates\tzero_based_inclusive",
				"#consensus_sha80\t"+o.get("consensussha80"), "#background_sha80\t"+o.get("backgroundsha80"),
				// HbmFamilyUnion records the singleton model registry as source_corpus.
				"#models_sha80\t"+DigestSuffix.fromDigest(provenance[12])));
			final ByteFile bf=ByteFile.makeByteFile(o.get("cores"), false);
			try{for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){if(line.length==0 || line[0]!='#'){break;} requiredHeaders.remove(new String(line, java.nio.charset.StandardCharsets.UTF_8));}}
			finally{require(!bf.close(), "Core metadata read failed");}
			require(requiredHeaders.isEmpty(), "Core metadata does not bind the requested model registry, profile, consensus and background");
			final List<String[]> bgRows=HmmComparisonData.rows(o.get("background"));
			require(bgRows.size()==21 && String.join("\t", bgRows.get(0)).equals("residue_code\tprobability"), "Invalid frozen background table"); final double[] bg=new double[20];
			for(int i=0; i<20; i++){final String[] row=bgRows.get(i+1); require(row.length==2 && Integer.parseInt(row[0])==i, "Invalid background residue ordering"); bg[i]=Double.parseDouble(row[1]);}
			final HbmBundleLoader.Loaded loaded=HbmBundleLoader.load(Paths.get(o.get("hbm")), roster, id->byId.get(id), provenance);
			gate=new HbmCompetitiveGate(ids, refs, loaded.positionModels("logodds", .01, true, -4, bg), first, last);
		}
		final byte[][] refs;
		final String[] ids;
		final int[] familyIds;
		final FamilyShortlistSidecar sidecar;
		final HbmCompetitiveGate gate;
	}
	private static void require(boolean ok, String message){if(!ok){throw new IllegalArgumentException(message);}}
	private static final String[] PINNED={"in", "consensus", "identity", "hbm", "provenance", "cores", "background", "sidecar", "kmers"};
}
