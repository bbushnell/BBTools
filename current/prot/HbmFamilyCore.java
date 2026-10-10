package prot;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
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
 * Recomputes half-member family cores in the final consensus coordinate frame.
 * Each member follows the final positional profile with the pinned background;
 * only paired m columns contribute depth. The core spans the first through last
 * positions with depth at least ceil(memberCount/2), matching the C0 core rule.
 * @author Brian Bushnell, Keqing
 */
public final class HbmFamilyCore {

	public static void main(String[] args){
		try{
			Shared.setThreads(1); Shared.AMINO_IN=true;
			PreParser pp=new PreParser(args, HbmFamilyCore.class, false);
			try{run(HmmComparisonData.options(pp.args));}finally{Shared.closeStream(pp.outstream);}
		}catch(Throwable failure){failure.printStackTrace(); System.exit(1);}
	}

	private static void run(HashMap<String,String> o) throws Exception{
		for(String key : o.keySet()){
			require(Arrays.asList("manifest", "manifestsha80", "models", "modelssha80", "consensus", "consensussha80",
				"background", "backgroundsha80", "ranks", "out").contains(key), "Unknown core parameter: "+key);
		}
		for(String key : new String[]{"manifest", "models", "consensus", "background"}){
			HbmProfilePilot.requireHash(HmmComparisonData.required(o, key), HmmComparisonData.required(o, key+"sha80"));
		}
		final Path out=Paths.get(HmmComparisonData.required(o, "out"));
		require(!Files.exists(out), "Core output must be a fresh directory: "+out);
		final List<String[]> members=HmmComparisonData.rows(o.get("manifest")), models=HmmComparisonData.rows(o.get("models"));
		require(members.size()>1 && String.join("\t", members.get(0)).equals(MANIFEST_HEADER), "Unsupported member manifest");
		require(models.size()==members.size() && String.join("\t", models.get(0)).equals(MODEL_HEADER), "Model/member manifests differ");
		final List<ProteinSequence> references=ProteinSearch.readFasta(o.get("consensus"));
		require(references.size()+1==members.size(), "Consensus/member family counts differ");
		final HashSet<String> ids=new HashSet<String>(), reps=new HashSet<String>();
		for(int i=0; i<references.size(); i++){
			final String[] m=members.get(i+1), h=models.get(i+1);
			require(m.length==7 && h.length==6 && m[0].equals(Integer.toString(i)) && h[0].equals(m[0]), "Noncontiguous family indexes at "+i);
			final int id=Integer.parseInt(m[1]);
			require(id>=0 && m[1].equals(Integer.toString(id)) && ids.add(m[1]) && reps.add(m[2]) && references.get(i).id.equals(m[2]) && h[1].equals(m[2]), "Family identity differs at "+i);
			require(Long.parseLong(m[3])>0 && Long.parseLong(m[4])>0 && h[5].equals(m[6]) && Paths.get(m[5]).isAbsolute(), "Invalid count/path/input binding at "+i);
		}
		final double[] background=background(o.get("background"));
		final int[] ranks=HbmProfileLibrary.readRanks(HmmComparisonData.required(o, "ranks"), references.size()); Arrays.sort(ranks);
		Files.createDirectory(out);
		final ByteStreamWriter cores=HmmComparisonData.writer(out.resolve("cores.tsv").toString());
		cores.println("#profile\tlogodds\tbeta\t0.01\tclipmin\t-4\tclipmax\t11\tgap\t4");
		cores.println("#coverage\tpaired_m_only\tthreshold\tceil(member_count/2)\tcoordinates\tzero_based_inclusive");
		for(String key : new String[]{"manifest", "models", "consensus", "background"}){cores.println("#"+key+"_sha80\t"+o.get(key+"sha80"));}
		cores.println("active_index\tfamily_id\trep_id\tmembers\tconsensus_length\tcoverage_threshold\tcore_start\tcore_end\tcore_length");
		long totalMembers=0, totalRaw=0;
		for(int rank : ranks){
			final String[] m=members.get(rank+1), h=models.get(rank+1);
			final Path dir=Paths.get(h[2]); final ProteinSequence reference=references.get(rank);
			HbmProfilePilot.requireHash(m[5], m[6]);
			HbmProfilePilot.requireHash(dir.resolve("refined.mqhb").toString(), h[3]);
			HbmProfilePilot.requireHash(dir.resolve("consensus.faa").toString(), h[4]);
			HbmProfilePilot.requireHash(dir.resolve("background.tsv").toString(), o.get("backgroundsha80"));
			final HbmBundleLoader.Loaded loaded=HbmBundleLoader.load(dir.resolve("refined.mqhb"), Collections.singletonList(m[2]), id->reference.enc,
				HbmBundleLoader.loadSemanticProvenance(dir.resolve("provenance.tsv").toString()));
			final HbmPositionModel profile=loaded.positionModels("logodds", 0.01, true, -4, background)[0];
			profile.requireConsensus(reference.enc);
			final long[] depth=new long[reference.length()];
			long count=0, raw=0;
			final FileFormat ff=FileFormat.testInput(m[5], FileFormat.FASTA, null, true, false);
			final Streamer input=StreamerFactory.makeStreamer(ff, null, true, -1); input.start();
			try{
				for(ListNum<Read> list=input.nextList(); list!=null && list.size()>0; list=input.nextList()){
					for(Read read : list){
						require(read.mate==null && read.bases!=null && read.bases.length>0, "Empty or paired core input in "+m[5]);
						final ProteinSequence query=new ProteinSequence(read.id, read.bases);
						accumulate(profile.align(query.enc, true), query.length(), depth);
						count++; raw=Math.addExact(raw, read.bases.length);
					}
				}
			}finally{input.close(); require(!input.errorState(), "Core member reader failed: "+m[5]);}
			require(count==Long.parseLong(m[3]) && raw==Long.parseLong(m[4]), "Core input count/residue mismatch at "+rank);
			HbmProfilePilot.requireHash(m[5], m[6]);
			final int[] core=findCore(depth, count); final long threshold=(count+1)/2;
			cores.println(new ByteBuilder().append(rank).tab().append(m[1]).tab().append(m[2]).tab().append(count).tab().append(reference.length())
				.tab().append(threshold).tab().append(core[0]).tab().append(core[1]).tab().append(core[1]-core[0]+1));
			totalMembers=Math.addExact(totalMembers, count); totalRaw=Math.addExact(totalRaw, raw);
			System.err.println("HBM_CORE_FAMILY_PASS rank="+rank+" members="+count+" core="+core[0]+".."+core[1]);
		}
		HmmComparisonData.close(cores);
		for(String key : new String[]{"manifest", "models", "consensus", "background"}){HbmProfilePilot.requireHash(o.get(key), o.get(key+"sha80"));}
		final ByteStreamWriter summary=HmmComparisonData.writer(out.resolve("summary.tsv").toString());
		summary.println("families\tmembers\traw_residues\tcores_sha80");
		summary.println(new ByteBuilder().append(ranks.length).tab().append(totalMembers).tab().append(totalRaw).tab().append(DigestSuffix.file(out.resolve("cores.tsv").toString())));
		HmmComparisonData.close(summary);
		final ByteStreamWriter pass=HmmComparisonData.writer(out.resolve("PASS").toString()); pass.println("HBM_FAMILY_CORE_PASS"); HmmComparisonData.close(pass);
		System.err.println("HBM_FAMILY_CORE_PASS families="+ranks.length+" members="+totalMembers);
	}

	/** Paired columns alone count; insertions/deletions cannot manufacture core coverage. */
	static void accumulate(HbmPositionModel.Result alignment, int queryLength, long[] depth){
		require(alignment!=null && alignment.path!=null && queryLength>0 && depth.length>0, "Core coverage requires a full recorded positional alignment");
		int q=0, r=alignment.start;
		require(r>=0 && r<depth.length, "Core alignment starts outside its reference");
		for(byte op : alignment.path){
			if(op=='m'){
				require(q<queryLength && r<depth.length, "Paired core column exceeds its sequence bounds"); depth[r]++; q++; r++;
			}else if(op=='I'){q++;}
			else if(op=='D'){r++;}
			else{throw new IllegalArgumentException("Unknown core path operation: "+(char)op);}
			require(q<=queryLength && r<=depth.length, "Core path exceeds sequence bounds");
		}
		require(q==queryLength && r==alignment.end+1, "Core path does not consume the declared query/reference spans");
	}

	/** C0 half-member endpoint rule, including low-depth internal positions between endpoints. */
	static int[] findCore(long[] depth, long members){
		require(members>0 && depth.length>0, "Core requires members and a nonempty reference");
		final long threshold=Math.addExact(members, 1)/2; int start=-1, end=-1;
		for(int i=0; i<depth.length; i++){
			require(depth[i]>=0 && depth[i]<=members, "Paired depth cannot exceed member count: a glocal path visits each reference column at most once");
			if(depth[i]>=threshold){if(start<0){start=i;} end=i;}
		}
		require(start>=0, "No reference column reaches the half-member core depth "+threshold);
		return new int[]{start, end};
	}

	private static double[] background(String path){
		final List<String[]> rows=HmmComparisonData.rows(path);
		require(rows.size()==21 && String.join("\t", rows.get(0)).equals("residue_code\tprobability"), "Background requires twenty ordered residue rows");
		final double[] out=new double[20]; double sum=0;
		for(int i=0; i<20; i++){
			final String[] row=rows.get(i+1); require(row.length==2 && row[0].equals(Integer.toString(i)), "Background residue order differs");
			out[i]=Double.parseDouble(row[1]); require(Double.isFinite(out[i]) && out[i]>0 && out[i]<=1, "Invalid background probability"); sum+=out[i];
		}
		require(Math.abs(sum-1)<1e-9, "Background probabilities must sum to one"); return out;
	}
	private static void require(boolean ok, String message){if(!ok){throw new IllegalArgumentException(message);}}
	private static final String MANIFEST_HEADER="active_index\tfamily_id\trep_id\tassigned_count\tsequence_bytes\tfasta_file\tfasta_sha80";
	private static final String MODEL_HEADER="active_index\trep_id\tdirectory\tbundle_sha80\tconsensus_sha80\tinput_sha80";
}
