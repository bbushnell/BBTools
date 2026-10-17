package prot;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import shared.Shared;
import structures.ByteBuilder;

/** Builds contiguous legacy HBM subset packs from the accepted 8,537-family package.
 * Slices preserve parent order exactly, renumber only package-local ranks, and bind the
 * derived model to a truthful subset lineage instead of inventing per-family artifacts.
 * The emitted files are raw; the launcher compresses package members and writes artifacts.tsv.
 */
public final class LegacySubsetPackBuilder {
	public static void main(String[] args){
		try{Shared.setThreads(1); run(HmmComparisonData.options(args));}
		catch(Throwable failure){failure.printStackTrace(); System.exit(1);}
	}

	private static void run(HashMap<String,String> o) throws Exception{
		for(String key : o.keySet()){
			require(Arrays.asList("in","out","ranges","runtime").contains(key),
				"Unknown subset-pack parameter: "+key);
		}
		final Path in=Paths.get(HmmComparisonData.required(o, "in"));
		final Path out=Paths.get(HmmComparisonData.required(o, "out"));
		final Path runtime=Paths.get(HmmComparisonData.required(o, "runtime"));
		final Pack[] packs=parseRanges(HmmComparisonData.required(o, "ranges"));
		require(Files.isDirectory(in), "Parent package directory is missing: "+in);
		require(!Files.exists(out), "Subset output must be fresh: "+out);
		final Parent parent=new Parent(in);
		require(packs[packs.length-1].end==parent.families()-1, "Subset ranges must cover through the last parent family");
		Files.createDirectory(out);
		final ByteStreamWriter manifest=HmmComparisonData.writer(out.resolve("packs.tsv").toString());
		manifest.println("pack\tstart\tend\tfamilies\tpath");
		for(Pack pack : packs){
			writePack(parent, pack, out.resolve(pack.name), runtime);
			manifest.println(new ByteBuilder().append(pack.name).tab().append(pack.start).tab()
				.append(pack.end).tab().append(pack.size()).tab().append(out.resolve(pack.name).toString()));
		}
		HmmComparisonData.close(manifest);
		final ByteStreamWriter pass=HmmComparisonData.writer(out.resolve("JAVA_PASS").toString());
		pass.println("LEGACY_SUBSET_PACKS_RAW_PASS packs="+packs.length);
		HmmComparisonData.close(pass);
		System.err.println("LEGACY_SUBSET_PACKS_RAW_PASS packs="+packs.length);
	}

	private static void writePack(Parent parent, Pack pack, Path out, Path runtime) throws Exception{
		require(!Files.exists(out), "Pack output must be fresh: "+out);
		require(pack.end<parent.families(), "Pack range exceeds parent family count: "+pack.end);
		Files.createDirectory(out);
		final ArrayList<String> roster=new ArrayList<String>(pack.size());
		final HashMap<String,byte[]> consensus=new HashMap<String,byte[]>();
		writeConsensusAndRows(parent, pack, out, roster, consensus);
		Files.copy(parent.in.resolve("background.tsv.gz"), out.resolve("background.tsv.gz"));
		copyKmers(parent, out);
		final String consensusPin=DigestSuffix.file(out.resolve("consensus.faa").toString());
		final String identityPin=DigestSuffix.file(out.resolve("family_identity.tsv").toString());
		final String artifactsPin=writeSourceCorpus(parent, pack, out);
		writeCores(parent, pack, out, consensusPin, artifactsPin);
		writeSidecar(parent, pack, out, consensusPin);
		writeRecipe(parent, pack, out);
		HbmProfilePilot.writeProvenance(out, runtime.toString(),
			out.resolve("family_artifacts.tsv").toString(),
			out.resolve("family_identity.tsv").toString(),
			parent.in.resolve("artifacts.tsv").toString(),
			"LegacySubsetPackBuilder");
		require(identityPin.equals(DigestSuffix.fromDigest(HbmBundleLoader.loadSemanticProvenance(
			out.resolve("provenance.tsv").toString())[13])), "Subset identity is not bound as cluster_membership");
		final HbmBundleLoader.Loaded subset=parent.loaded.subset(pack.start, pack.end+1);
		subset.writeCompactBundle(out.resolve("model.hbmc"), readProvenancePins(out.resolve("provenance.tsv")));
		final HbmBundleLoader.Loaded reloaded=HbmBundleLoader.load(out.resolve("model.hbmc"), roster, id->consensus.get(id),
			HbmBundleLoader.loadSemanticProvenance(out.resolve("provenance.tsv").toString()));
		reloaded.assertStructuralEquivalent(subset);
		require(Arrays.equals(FamilyShortlistSidecarBuilder.loadConsensusFamilyOrder(out.resolve("consensus.faa").toString()),
			roster.toArray(new String[0])), "Subset consensus order differs");
		final HashMap<String,String> opts=new HashMap<String,String>();
		opts.put("hbm", out.resolve("model.hbmc").toString());
		opts.put("consensus", out.resolve("consensus.faa").toString());
		opts.put("provenance", out.resolve("provenance.tsv").toString());
		opts.put("background", out.resolve("background.tsv.gz").toString());
		opts.put("identity", out.resolve("family_identity.tsv").toString());
		opts.put("cores", out.resolve("cores.tsv").toString());
		opts.put("kmers", out.resolve("kmers.tsv.gz").toString());
		opts.put("sidecar", out.resolve("sidecar.tsv").toString());
		opts.put("identitysha80", identityPin);
		opts.put("consensussha80", consensusPin);
		opts.put("backgroundsha80", parent.artifacts.get("background.tsv").rawPin);
		require(new HbmCompetitiveAssign.Library(opts).ids.length==pack.size(), "Produced pack is not loadable");
	}

	private static void writeConsensusAndRows(Parent parent, Pack pack, Path out,
			ArrayList<String> roster, HashMap<String,byte[]> consensus){
		final ByteStreamWriter fasta=HmmComparisonData.writer(out.resolve("consensus.faa").toString());
		final ByteStreamWriter identity=HmmComparisonData.writer(out.resolve("family_identity.tsv").toString());
		final ByteStreamWriter ranks=HmmComparisonData.writer(out.resolve("roster.tsv").toString());
		identity.println("active_index\tfamily_id\tsource_active_index\trep_id\tmembers\traw_residues");
		ranks.println("rank\trep_id\tmembers");
		for(int source=pack.start, local=0; source<=pack.end; source++, local++){
			final String[] id=parent.identity.get(source+1), rank=parent.roster.get(source+1);
			require(id.length==6 && rank.length==3 && id[0].equals(Integer.toString(source))
				&& rank[0].equals(id[0]) && rank[1].equals(id[3]) && rank[2].equals(id[4]),
				"Parent identity/roster mismatch at "+source);
			final ProteinSequence sequence=parent.consensus.get(source);
			require(sequence.id.equals(id[3]), "Parent consensus order differs at "+source);
			HmmComparisonData.writeFasta(fasta, sequence.id+" active_index="+local+" family_id="+id[1], sequence.enc);
			identity.println(new ByteBuilder().append(local).tab().append(id[1]).tab().append(id[2]).tab()
				.append(id[3]).tab().append(id[4]).tab().append(id[5]));
			ranks.println(new ByteBuilder().append(local).tab().append(id[3]).tab().append(id[4]));
			roster.add(id[3]); consensus.put(id[3], sequence.enc);
		}
		HmmComparisonData.close(fasta); HmmComparisonData.close(identity); HmmComparisonData.close(ranks);
	}

	private static String writeSourceCorpus(Parent parent, Pack pack, Path out){
		final ByteStreamWriter w=HmmComparisonData.writer(out.resolve("family_artifacts.tsv").toString());
		w.println("active_index\tsource_active_index\tfamily_id\trep_id\tmembers\traw_residues\tparent_model_sha80\tparent_consensus_sha80\tparent_identity_sha80");
		for(int source=pack.start, local=0; source<=pack.end; source++, local++){
			final String[] id=parent.identity.get(source+1);
			w.println(new ByteBuilder().append(local).tab().append(source).tab().append(id[1]).tab()
				.append(id[3]).tab().append(id[4]).tab().append(id[5]).tab()
				.append(parent.artifacts.get("model.hbmc").rawPin).tab()
				.append(parent.artifacts.get("consensus.faa").rawPin).tab()
				.append(parent.artifacts.get("family_identity.tsv").rawPin));
		}
		HmmComparisonData.close(w);
		return DigestSuffix.file(out.resolve("family_artifacts.tsv").toString());
	}

	private static void writeCores(Parent parent, Pack pack, Path out, String consensusPin, String artifactsPin){
		final ByteStreamWriter w=HmmComparisonData.writer(out.resolve("cores.tsv").toString());
		for(String line : parent.cores){
			if(line.startsWith("#consensus_sha80\t")){w.println("#consensus_sha80\t"+consensusPin);}
			else if(line.startsWith("#models_sha80\t")){w.println("#models_sha80\t"+artifactsPin);}
			else if(line.startsWith("#")){w.println(line);}
			else{
				final String[] row=line.split("\t", -1);
				if(row.length>0 && row[0].equals("active_index")){w.println(line); continue;}
				require(row.length==9, "Malformed parent cores row");
				final int source=Integer.parseInt(row[0]);
				if(source>=pack.start && source<=pack.end){
					final String[] local=row.clone();
					local[0]=Integer.toString(source-pack.start); w.println(String.join("\t", local));
				}
			}
		}
		HmmComparisonData.close(w);
	}

	private static void writeSidecar(Parent parent, Pack pack, Path out, String consensusPin){
		final ByteFile bf=ByteFile.makeByteFile(parent.in.resolve("sidecar.tsv.gz").toString(), false);
		final ByteStreamWriter w=HmmComparisonData.writer(out.resolve("sidecar.tsv").toString());
		int source=0, copied=0; boolean header=true;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				final String s=new String(line, java.nio.charset.StandardCharsets.US_ASCII);
				if(header && s.startsWith("#consensus_sha80\t")){w.println("#consensus_sha80\t"+consensusPin);}
				else if(header && s.startsWith("#n_families\t")){w.println("#n_families\t"+pack.size());}
				else if(header && (s.startsWith("#") || s.startsWith("rank\t"))){w.println(s); if(s.startsWith("rank\t")){header=false;}}
				else{
					if(source>=pack.start && source<=pack.end){
						final int tab=s.indexOf('\t');
						require(tab>0, "Malformed parent sidecar row");
						require(Integer.parseInt(s.substring(0, tab))==source, "Parent sidecar rank order differs");
						require(s.startsWith("\t"+parent.consensus.get(source).id+"\t", tab), "Parent sidecar ID differs at "+source);
						w.println(new ByteBuilder().append(source-pack.start).append(s.substring(tab)));
						copied++;
					}
					source++;
				}
			}
		}finally{require(!bf.close(), "Parent sidecar read failed");}
		HmmComparisonData.close(w);
		require(copied==pack.size(), "Sidecar subset count differs");
		require(source==parent.families(), "Parent sidecar count differs");
	}

	private static void copyKmers(Parent parent, Path out) throws Exception{
		Files.copy(parent.in.resolve("kmers.tsv.gz"), out.resolve("kmers.tsv.gz"));
	}

	private static void writeRecipe(Parent parent, Pack pack, Path out){
		final ByteStreamWriter w=HmmComparisonData.writer(out.resolve("recipe.txt").toString());
		w.println("legacy_subset_pack_v1");
		w.println("parent_package="+parent.in.toString());
		w.println("start="+pack.start);
		w.println("end="+pack.end);
		w.println("families="+pack.size());
		w.println("parent_model_sha80="+parent.artifacts.get("model.hbmc").rawPin);
		w.println("parent_consensus_sha80="+parent.artifacts.get("consensus.faa").rawPin);
		w.println("parent_identity_sha80="+parent.artifacts.get("family_identity.tsv").rawPin);
		HmmComparisonData.close(w);
	}

	private static String[] readProvenancePins(Path path){
		final List<String[]> rows=HmmComparisonData.rows(path.toString());
		require(rows.size()==HbmBundleFormat.PROVENANCE_COUNT, "Expected16 provenance rows in "+path);
		final String[] out=new String[HbmBundleFormat.PROVENANCE_COUNT];
		for(int i=0; i<out.length; i++){out[i]=DigestSuffix.requireSuffix(rows.get(i)[1], rows.get(i)[0]);}
		return out;
	}

	private static Pack[] parseRanges(String text){
		final String[] parts=text.split(",");
		final Pack[] out=new Pack[parts.length];
		require(out.length>0, "At least one pack range is required");
		final HashSet<String> names=new HashSet<String>();
		int previous=-1;
		for(int i=0; i<parts.length; i++){
			final String[] nameRange=parts[i].split(":", -1);
			require(nameRange.length==2 && names.add(nameRange[0]) && nameRange[0].matches("[A-Za-z0-9._-]+"),
				"Bad or duplicate pack range: "+parts[i]);
			final int dash=nameRange[1].indexOf('-');
			require(dash>0, "Range must use inclusive start-end: "+parts[i]);
			final int start=Integer.parseInt(nameRange[1].substring(0, dash));
			final int end=Integer.parseInt(nameRange[1].substring(dash+1));
			require(start==previous+1 && end>=start, "Pack ranges must be contiguous in source order: "+parts[i]);
			out[i]=new Pack(nameRange[0], start, end); previous=end;
		}
		return out;
	}

	private static final class Parent{
		Parent(Path in_) throws Exception{
			in=in_;
			ProteinFamilyAssigner.load(in);
			artifacts=readArtifacts(in.resolve("artifacts.tsv"));
			require(DigestSuffix.file(in.resolve("model.hbmc.gz").toString()).equals(artifacts.get("model.hbmc").compressedPin), "Parent model stored bytes differ");
			require(DigestSuffix.file(in.resolve("provenance.tsv.gz").toString()).equals(artifacts.get("provenance.tsv").compressedPin), "Parent provenance stored bytes differ");
			consensus=ProteinSearch.readFasta(in.resolve("consensus.faa.gz").toString());
			identity=HmmComparisonData.rows(in.resolve("family_identity.tsv.gz").toString());
			roster=HmmComparisonData.rows(in.resolve("roster.tsv.gz").toString());
			cores=lines(in.resolve("cores.tsv.gz"));
			require(consensus.size()>0 && identity.size()==consensus.size()+1 && roster.size()==identity.size(),
				"Parent package row counts differ");
			final ArrayList<String> ids=new ArrayList<String>(consensus.size());
			final HashMap<String,byte[]> byId=new HashMap<String,byte[]>();
			for(ProteinSequence sequence : consensus){
				ids.add(sequence.id); require(byId.put(sequence.id, sequence.enc)==null, "Duplicate parent consensus: "+sequence.id);
			}
			loaded=HbmBundleLoader.load(in.resolve("model.hbmc.gz"), ids, id->byId.get(id),
				HbmBundleLoader.loadSemanticProvenance(in.resolve("provenance.tsv.gz").toString()));
			require(loaded.familyCount()==consensus.size(), "Parent HBM family count differs");
		}
		int families(){return consensus.size();}
		final Path in;
		final HashMap<String,Artifact> artifacts;
		final List<ProteinSequence> consensus;
		final List<String[]> identity, roster;
		final List<String> cores;
		final HbmBundleLoader.Loaded loaded;
	}

	private static List<String> lines(Path path){
		final ArrayList<String> out=new ArrayList<String>();
		final ByteFile bf=ByteFile.makeByteFile(path.toString(), false);
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				out.add(new String(line, StandardCharsets.UTF_8));
			}
		}finally{require(!bf.close(), "Input read failed: "+path);}
		return out;
	}

	private static HashMap<String,Artifact> readArtifacts(Path path){
		final List<String[]> rows=HmmComparisonData.rows(path.toString());
		final HashMap<String,Artifact> out=new HashMap<String,Artifact>();
		require(rows.size()==10 && String.join("\t", rows.get(0)).equals("artifact\traw_bytes\traw_sha80\tcompressed_bytes\tcompressed_sha80\tcodec"),
			"Unsupported artifact manifest");
		for(int i=1; i<rows.size(); i++){
			final String[] r=rows.get(i);
			require(r.length==6, "Malformed artifact row");
			out.put(r[0], new Artifact(DigestSuffix.requireSuffix(r[2], r[0]), DigestSuffix.requireSuffix(r[4], r[0])));
		}
		return out;
	}

	private static final class Artifact{
		Artifact(String rawPin_, String compressedPin_){rawPin=rawPin_; compressedPin=compressedPin_;}
		final String rawPin, compressedPin;
	}
	private static final class Pack{
		Pack(String name_, int start_, int end_){name=name_; start=start_; end=end_;}
		int size(){return end-start+1;}
		final String name;
		final int start, end;
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
}
