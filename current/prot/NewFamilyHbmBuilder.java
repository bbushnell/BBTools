package prot;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

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
 * Builds one experimental HBM from every member of an original cluster. BLOSUM
 * placement on the original representative supplies the seed graph, followed by
 * two frozen positional-profile passes in HbmProfileRefiner. Input identity,
 * counts, raw FASTA residue bytes and a common background are pinned explicitly.
 * Terminal '*' markers are retained in the raw input count but are not encoded
 * amino acids. Internal stops fail; no member is silently skipped or sampled.
 * @author Brian Bushnell, Keqing
 */
public final class NewFamilyHbmBuilder {

	public static void main(String[] args){
		try{
			Shared.setThreads(1); Shared.AMINO_IN=true;
			PreParser pp=new PreParser(args, NewFamilyHbmBuilder.class, false);
			try{run(HmmComparisonData.options(pp.args));}finally{Shared.closeStream(pp.outstream);}
		}catch(Throwable failure){failure.printStackTrace(); System.exit(1);}
	}

	private static void run(HashMap<String,String> o) throws Exception{
		for(String key : o.keySet()){
			if(!Arrays.asList("manifest", "manifestsha80", "source", "rank", "background", "backgroundsha80", "out", "runtime", "padding").contains(key)){
				throw new IllegalArgumentException("Unknown seed-builder parameter: "+key);
			}
		}
		final String manifest=HmmComparisonData.required(o, "manifest"), manifestPin=HmmComparisonData.required(o, "manifestsha80");
		final String backgroundFile=HmmComparisonData.required(o, "background"), backgroundPin=HmmComparisonData.required(o, "backgroundsha80");
		final String runtime=HmmComparisonData.required(o, "runtime");
		final int rank=Integer.parseInt(HmmComparisonData.required(o, "rank"));
		final Path out=Paths.get(HmmComparisonData.required(o, "out"));
		require(!Files.exists(out), "Seed output must be a fresh directory: "+out);
		HbmProfilePilot.requireHash(manifest, manifestPin); HbmProfilePilot.requireHash(backgroundFile, backgroundPin);
		final String[] entry=findEntry(manifest, rank);
		final Path input=Paths.get(HmmComparisonData.required(o, "source")).resolve(entry[5]);
		final int expected=Integer.parseInt(entry[3]);
		final long expectedRaw=Long.parseLong(entry[4]);
		HbmProfilePilot.requireHash(input.toString(), entry[6]);
		final double[] background=readBackground(backgroundFile);
		final ArrayList<ProteinSequence> membersList=new ArrayList<ProteinSequence>();
		final HashSet<String> ids=new HashSet<String>();
		long raw=0, residues=0;
		byte[] representative=null;
		int maxLength=0;
		final FileFormat ff=FileFormat.testInput(input.toString(), FileFormat.FASTA, null, true, false);
		final Streamer st=StreamerFactory.makeStreamer(ff, null, true, -1);
		st.start();
		try{
			for(ListNum<Read> list=st.nextList(); list!=null && list.size()>0; list=st.nextList()){
				for(Read read : list){
					require(read.mate==null && ids.add(read.id), "Paired or duplicate seed protein: "+read.id);
					raw=Math.addExact(raw, read.bases.length);
					final ProteinSequence p=new ProteinSequence(read.id, read.bases);
					membersList.add(p); residues=Math.addExact(residues, p.length()); maxLength=Math.max(maxLength, p.length());
					if(p.id.equals(entry[2])){representative=p.enc;}
				}
			}
		}finally{st.close(); require(!st.errorState(), "Seed FASTA reader reported an I/O error");}
		require(membersList.size()==expected && raw==expectedRaw && representative!=null,
			"Seed membership differs from manifest: count="+membersList.size()+" expected="+expected+" raw="+raw+" expectedRaw="+expectedRaw+" representative="+(representative!=null));
		final byte[][] members=new byte[membersList.size()][];
		for(int i=0; i<members.length; i++){members[i]=membersList.get(i).enc;}
		final AAGraph seed=new AAGraph(representative, 0);
		long seedExcluded=0;
		for(byte[] member : members){
			final AAAlignment alignment=GlocalAminoLinear.align(member, representative, true);
			require(alignment.qStart==0 && alignment.qStop==member.length-1, "Seed alignment must consume the full query for addTrace");
			seedExcluded=Math.addExact(seedExcluded, seed.addTrace(member, alignment.tStart, alignment.match));
		}
		HbmBundleBuilder.validateGraph(new HbmBundleBuilder.FamilyInput(entry[2], representative, seed));
		checkResidues(seed, residues-seedExcluded);
		int padding=Integer.parseInt(o.getOrDefault("padding", "20")), retries=0;
		require(padding>=0, "Negative growth padding");
		HbmProfileRefiner.Result result;
		while(true){
			try{result=HbmProfileRefiner.refine(seed, members, background, padding); break;}
			catch(HbmProfileRefiner.PaddingException small){
				if(padding>=maxLength){throw small;}
				padding=(int)Math.min(maxLength, Math.max(2L*padding, (long)padding+small.excluded)); retries++;
				System.err.println("SEED_PADDING_RETRY rank="+rank+" padding="+padding);
			}
		}
		HbmProfilePilot.requireHash(input.toString(), entry[6]); HbmProfilePilot.requireHash(manifest, manifestPin);
		HbmProfilePilot.requireHash(backgroundFile, backgroundPin);
		Files.createDirectory(out);
		ByteStreamWriter bw=HmmComparisonData.writer(out.resolve("consensus.faa").toString());
		HmmComparisonData.writeFasta(bw, entry[2], result.graph.pivot); HmmComparisonData.close(bw);
		bw=HmmComparisonData.writer(out.resolve("roster.tsv").toString());
		bw.println("rank\trep_id\tmembers"); bw.println(new ByteBuilder().append(0).tab().append(entry[2]).tab().append(expected)); HmmComparisonData.close(bw);
		Files.copy(Paths.get(backgroundFile), out.resolve("background.tsv"));
		bw=HmmComparisonData.writer(out.resolve("recipe.txt").toString());
		bw.println("experimental_new_family_hbm_v1\nseed=BLOSUM62_glocal_gap4_original_representative\nmembership=all_no_cap\nbeta=0.01\nclipmin=-4\nclipmax=11\ngap=4\ntrimdepth=0.1\nprofile_passes=2\nfinal_padding=0\nterminal_policy=count_exclusions");
		bw.println("source_active_index="+rank+"\nfamily_id="+entry[1]+"\nrep_id="+entry[2]+"\ninput_sha80="+entry[6]+"\nmanifest_sha80="+manifestPin);
		bw.println("background_sha80="+backgroundPin+"\ngrowth_padding="+padding+"\npadding_retries="+retries);
		for(String dependency : new String[]{"HbmPositionModel", "HbmProfileRefiner", "HbmProfilePilot"}){
			bw.println(dependency+"_source_sha80="+DigestSuffix.file(runtime+"/current/prot/"+dependency+".java"));
			bw.println(dependency+"_class_sha80="+DigestSuffix.file(runtime+"/current/prot/"+dependency+".class"));
		}
		HmmComparisonData.close(bw);
		final byte[][] provenance=HbmProfilePilot.writeProvenance(out, runtime, input.toString(), input.toString(), manifest, "NewFamilyHbmBuilder");
		final Path bundle=out.resolve("refined.mqhb");
		HbmBundleBuilder.build(bundle, Collections.singletonList(new HbmBundleBuilder.FamilyInput(entry[2], result.graph.pivot, result.graph)), provenance);
		final byte[] consensus=result.graph.pivot;
		final HbmBundleLoader.Loaded loaded=HbmBundleLoader.load(bundle, Collections.singletonList(entry[2]), id->consensus,
			HbmBundleLoader.loadSemanticProvenance(out.resolve("provenance.tsv").toString()));
		loaded.assertStructuralMatch(0, result.graph);
		bw=HmmComparisonData.writer(out.resolve("summary.tsv").toString());
		bw.println("active_index\tfamily_id\trep_id\tmembers\traw_residue_bytes\tencoded_residues\tseed_length\tfinal_length\tseed_terminal_excluded\tfinal_terminal_excluded\tbundle_sha80\tconsensus_sha80");
		bw.println(new ByteBuilder().append(rank).tab().append(entry[1]).tab().append(entry[2]).tab().append(expected).tab().append(raw).tab().append(residues)
			.tab().append(representative.length).tab().append(consensus.length).tab().append(seedExcluded).tab().append(result.excludedTerminalResidues)
			.tab().append(DigestSuffix.file(bundle.toString())).tab().append(DigestSuffix.file(out.resolve("consensus.faa").toString())));
		HmmComparisonData.close(bw);
		bw=HmmComparisonData.writer(out.resolve("PASS").toString()); bw.println("NEW_FAMILY_HBM_PASS"); HmmComparisonData.close(bw);
		System.err.println("NEW_FAMILY_HBM_PASS rank="+rank+" members="+expected+" length="+consensus.length);
	}

	private static String[] findEntry(String manifest, int rank){
		final List<String[]> rows=HmmComparisonData.rows(manifest);
		require(!rows.isEmpty() && String.join("\t", rows.get(0)).equals("active_index\tfamily_id\trep_id\tassigned_count\tsequence_bytes\tfasta_file\tfasta_sha80"), "Unsupported seed manifest header");
		String[] result=null;
		final HashSet<String> ranks=new HashSet<String>(), ids=new HashSet<String>(), reps=new HashSet<String>();
		for(int i=1; i<rows.size(); i++){
			String[] row=rows.get(i);
			require(row.length==7, "Seed manifest must have seven fields");
			int active=Integer.parseInt(row[0]), id=Integer.parseInt(row[1]);
			require(active>=0 && id>=0 && row[0].equals(Integer.toString(active)) && row[1].equals(Integer.toString(id)), "Noncanonical seed manifest IDs");
			require(ranks.add(row[0]) && ids.add(row[1]) && reps.add(row[2]) && Integer.parseInt(row[3])>0 && Long.parseLong(row[4])>0,
				"Duplicate or invalid seed manifest identity/count at active index "+active);
			require(row[5].equals("rank_"+active+".faa"), "Unexpected seed FASTA filename: "+row[5]);
			DigestSuffix.requireSuffix(row[6], "member FASTA sha80");
			if(active==rank){result=row;}
		}
		require(result!=null, "Requested active index absent from seed manifest: "+rank);
		return result;
	}

	private static double[] readBackground(String file){
		final List<String[]> rows=HmmComparisonData.rows(file);
		require(rows.size()==21 && String.join("\t", rows.get(0)).equals("residue_code\tprobability"), "Background requires twenty ordered residue probabilities");
		final double[] background=new double[20]; double sum=0;
		for(int i=0; i<20; i++){
			String[] row=rows.get(i+1);
			require(row.length==2 && Integer.parseInt(row[0])==i, "Background residue order differs at "+i);
			background[i]=Double.parseDouble(row[1]);
			require(Double.isFinite(background[i]) && background[i]>0 && background[i]<=1, "Invalid background probability at "+i);
			sum+=background[i];
		}
		require(Math.abs(sum-1)<1e-9, "Background probabilities must sum to one: "+sum);
		return background;
	}

	/** Counts each placed query residue once, excluding the graph's one scaffold observation per column. */
	private static void checkResidues(AAGraph graph, long expected){
		long actual=0;
		for(int i=0; i<graph.ref.length; i++){
			actual+=graph.ref[i].countSum-1L;
			for(AAGraphNode node=graph.ref[i].insEdge; node!=null; node=node.insEdge){actual+=node.countSum;}
			for(AAGraphNode node=graph.del[i].insEdge; node!=null; node=node.insEdge){actual+=node.countSum;}
		}
		require(actual==expected, "Seed graph residue conservation failed: actual="+actual+" expected="+expected);
	}
	private static void require(boolean ok, String message){if(!ok){throw new IllegalArgumentException(message);}}
}
