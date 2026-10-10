package prot;

import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import fileIO.ByteStreamWriter;
import structures.ByteBuilder;

/** Experimental one-family driver on pinned members and an existing HBM library; never installs outputs. */
public final class HbmProfilePilot {
	public static void main(String[] args){
		try{run(HmmComparisonData.options(args));}
		catch(Throwable failure){failure.printStackTrace();System.exit(1);}
	}
	private static void run(HashMap<String,String> o) throws Exception{
		for(String key : o.keySet()){
			if(!Arrays.asList("in","resources","manifest","manifestsha80","rank","out","runtime","padding").contains(key)){
				throw new IllegalArgumentException("Unknown pilot parameter: "+key);
			}
		}
		final Session session=new Session(o,true);
		runFamily(o,session,true);
	}
	/** A shard keeps the immutable parent library and its original background in memory once. */
	static final class Session{
		final String manifest,manifestPin,runtime,parent,parentPin;
		final List<String[]> rows;
		final ArrayList<String> roster=new ArrayList<String>();
		final HashMap<String,byte[]> sequence=new HashMap<String,byte[]>();
		final HbmBundleLoader.Loaded loaded;
		final double[] background;
		Session(HashMap<String,String> o,boolean loadParent) throws Exception{
			final String resources=HmmComparisonData.required(o,"resources");
			manifest=HmmComparisonData.required(o,"manifest");manifestPin=HmmComparisonData.required(o,"manifestsha80");
			runtime=HmmComparisonData.required(o,"runtime");requireHash(manifest,manifestPin);
			rows=HmmComparisonData.rows(manifest);
			if(rows.isEmpty() || !String.join("\t",rows.get(0)).equals("active_index\tfamily_id\trep_id\tassigned_count\tsequence_bytes\tfasta_file\tfasta_sha80")){
				throw new IllegalArgumentException("Unexpected pilot family-manifest schema");
			}
			for(ProteinSequence reference : ProteinSearch.readFasta(resources+"/consensus_reps_round2.faa.gz")){
				if(sequence.put(reference.id,reference.enc)!=null){throw new IllegalArgumentException("Duplicate HBM reference");}
				roster.add(reference.id);
			}
			if(rows.size()!=roster.size()+1){throw new IllegalArgumentException("Family manifest and parent roster sizes differ");}
			for(int i=0; i<roster.size(); i++){
				final String[] row=rows.get(i+1);
				if(row.length!=7 || Integer.parseInt(row[0])!=i || !row[2].equals(roster.get(i)) || Long.parseLong(row[3])<1 || Long.parseLong(row[4])<1
					|| !row[5].equals(String.format(java.util.Locale.ROOT,"rank_%04d.faa",i))){
					throw new IllegalArgumentException("Family manifest and parent roster differ at rank "+i);
				}
			}
			parent=resources+"/magqc_hbm_v1.rare01.hbmt.gz";parentPin=DigestSuffix.file(parent);
			loaded=loadParent ? HbmBundleLoader.load(Paths.get(parent),roster,id->sequence.get(id),
				HbmBundleLoader.loadSemanticProvenance(resources+"/PROVENANCE_MANIFEST.tsv.gz")) : null;
			background=loadParent ? loaded.positionBackground() : null;
		}
		void verify(){requireHash(manifest,manifestPin);requireHash(parent,parentPin);}
	}
	static void runFamily(HashMap<String,String> o,Session session,boolean verifySession) throws Exception{
		if(session.loaded==null){throw new IllegalArgumentException("Family refinement requires a loaded parent library");}
		final String in=HmmComparisonData.required(o,"in"),manifest=session.manifest,runtime=session.runtime;
		final List<String[]> rows=session.rows;
		final Path out=Paths.get(HmmComparisonData.required(o,"out"));
		if(Files.exists(out)){throw new IllegalArgumentException("Pilot output must be a fresh directory: "+out);}
		final int rank=Integer.parseInt(HmmComparisonData.required(o,"rank"));
		int padding=Integer.parseInt(o.getOrDefault("padding","20"));
		if(rank<0 || padding<0){throw new IllegalArgumentException("Rank and growth padding must be nonnegative");}
		if(rank+1>=rows.size()){throw new IllegalArgumentException("Pilot rank outside manifest");}
		final String[] entry=rows.get(rank+1);
		if(entry.length!=7 || Integer.parseInt(entry[0])!=rank){throw new IllegalArgumentException("Noncontiguous pilot family rank");}
		final String rep=entry[2];final int expected=Integer.parseInt(entry[3]);final long expectedResidues=Long.parseLong(entry[4]);
		requireHash(in,entry[6]);
		final List<ProteinSequence> proteins=ProteinSearch.readFasta(in);
		if(proteins.size()!=expected){throw new IllegalArgumentException("Pilot member count differs from pinned manifest");}
		final byte[][] members=new byte[proteins.size()][];final HashSet<String> seen=new HashSet<String>();
		long residues=0;int maxLength=0;
		for(int i=0; i<members.length; i++){
			final ProteinSequence p=proteins.get(i);
			if(!seen.add(p.id)){throw new IllegalArgumentException("Duplicate family member: "+p.id);}
			members[i]=p.enc;residues=Math.addExact(residues,p.enc.length);maxLength=Math.max(maxLength,p.enc.length);
		}
		if(residues!=expectedResidues){throw new IllegalArgumentException("Pilot residue count differs from pinned manifest: "+residues+" != "+expectedResidues);}
		final HbmBundleLoader.Loaded loaded=session.loaded;
		final double[] background=session.background;
		HbmProfileRefiner.Result result;int retries=0;
		while(true){
			try{result=loaded.refineFamily(rank,members,background,padding);break;}
			catch(HbmProfileRefiner.PaddingException small){
				if(padding>=maxLength){throw small;}
				padding=(int)Math.min(maxLength,Math.max(2L*padding,(long)padding+small.excluded));
				retries++;System.err.println("PROFILE_PADDING_RETRY rank="+rank+" padding="+padding);
			}
		}
		requireHash(in,entry[6]);if(verifySession){session.verify();}
		Files.createDirectory(out);
		final ByteStreamWriter fasta=HmmComparisonData.writer(out.resolve("consensus.faa").toString());
		HmmComparisonData.writeFasta(fasta,rep,result.graph.pivot);HmmComparisonData.close(fasta);
		final ByteStreamWriter ranks=HmmComparisonData.writer(out.resolve("roster.tsv").toString());
		ranks.println("rank\trep_id\tmembers");ranks.println(new ByteBuilder().append(0).tab().append(rep).tab().append(expected));HmmComparisonData.close(ranks);
		final ByteStreamWriter bg=HmmComparisonData.writer(out.resolve("background.tsv").toString());
		bg.println("residue_code\tprobability");
		for(int a=0; a<20; a++){bg.println(new ByteBuilder().append(a).tab().append(Double.toString(background[a])));}
		HmmComparisonData.close(bg);
		final ByteStreamWriter recipe=HmmComparisonData.writer(out.resolve("recipe.txt").toString());
		recipe.println("experimental_profile_refinement_v1\nbeta=0.01\nclipmin=-4\nclipmax=11\ngap=4\ntrimdepth=0.1\npasses=2");
		recipe.println("growth_padding="+padding+"\nfinal_padding=0\nbackground=parent_library_ref_counts_plus_one\nfinal_terminal_policy=count_exclusions");
		recipe.println("parent_hbm_sha80="+session.parentPin+"\ninput_sha80="+entry[6]+"\nsource_rank="+rank+"\nsource_rep_id="+rep);
		recipe.println("runtime_source_manifest_sha80="+DigestSuffix.file(runtime+"/POSITION_SOURCE.tsv"));
		recipe.println("background_sha80="+DigestSuffix.file(out.resolve("background.tsv").toString()));
		HmmComparisonData.close(recipe);
		final byte[][] provenance=writeProvenance(out,runtime,in,in,manifest,"HbmProfilePilot");
		final Path bundle=out.resolve("refined.mqhb");
		HbmBundleBuilder.build(bundle,Collections.singletonList(new HbmBundleBuilder.FamilyInput(rep,result.graph.pivot,result.graph)),provenance);
		final byte[] consensus=result.graph.pivot;
		final HbmBundleLoader.Loaded roundtrip=HbmBundleLoader.load(bundle,Collections.singletonList(rep),id->consensus,
			HbmBundleLoader.loadSemanticProvenance(out.resolve("provenance.tsv").toString()));
		roundtrip.assertStructuralMatch(0,result.graph);
		final ByteStreamWriter summary=HmmComparisonData.writer(out.resolve("summary.tsv").toString());
		summary.println("source_rank\trep_id\tmembers\tresidues\told_length\tnew_length\tchanged_consensus\tpadding\tpadding_retries\tfinal_terminal_excluded\tbundle_sha80\tconsensus_sha80");
		summary.println(new ByteBuilder().append(rank).tab().append(rep).tab().append(expected).tab().append(residues).tab().append(session.sequence.get(rep).length)
			.tab().append(consensus.length).tab().append(Boolean.toString(!Arrays.equals(session.sequence.get(rep),consensus))).tab().append(padding).tab().append(retries)
			.tab().append(result.excludedTerminalResidues).tab().append(DigestSuffix.file(bundle.toString())).tab().append(DigestSuffix.file(out.resolve("consensus.faa").toString())));
		HmmComparisonData.close(summary);
		final ByteStreamWriter pass=HmmComparisonData.writer(out.resolve("PASS").toString());pass.println("PROFILE_PILOT_PASS");HmmComparisonData.close(pass);
		System.err.println("PROFILE_PILOT_PASS rank="+rank+" members="+expected);
	}
	/** Binds each new artifact to its actual producer, inputs and coordinate-frame outputs. */
	static byte[][] writeProvenance(Path out,String runtime,String corpus,String membership,String manifest,String producer) throws Exception{
		final String[] keys={"format_spec","loader_source","loader_class","aagraph","aagraphnode","aagraphscorer","glocalaminolinear","blosum62",
			"builder_source","builder_class","roster","consensus_ref","source_corpus","cluster_membership","member_policy","member_manifest"};
		final String[] files=provenanceFiles(out,runtime,corpus,membership,manifest,producer);
		final byte[][] provenance=new byte[keys.length][];
		final ByteStreamWriter prov=HmmComparisonData.writer(out.resolve("provenance.tsv").toString());
		for(int i=0; i<keys.length; i++){
			provenance[i]=digest(files[i]);
			prov.println(new ByteBuilder().append(keys[i]).append("_sha80\t").append(DigestSuffix.fromDigest(provenance[i])).tab().append(files[i]));
		}
		HmmComparisonData.close(prov);
		return provenance;
	}
	static void verifyProvenance(Path out,String runtime,String corpus,String membership,String manifest,String producer){
		final byte[][] declared=HbmBundleLoader.loadSemanticProvenance(out.resolve("provenance.tsv").toString());
		final String[] files=provenanceFiles(out,runtime,corpus,membership,manifest,producer);
		for(int i=0; i<files.length; i++){
			if(!DigestSuffix.file(files[i]).equals(DigestSuffix.fromDigest(declared[i]))){throw new IllegalArgumentException("Refined family provenance does not match actual input/artifact at field "+i+": "+out);}
		}
	}
	private static String[] provenanceFiles(Path out,String runtime,String corpus,String membership,String manifest,String producer){
		return new String[]{runtime+"/current/prot/HbmBundleFormat.java",runtime+"/current/prot/HbmBundleLoader.java",runtime+"/current/prot/HbmBundleLoader.class",
			runtime+"/current/prot/AAGraph.java",runtime+"/current/prot/AAGraphNode.java",runtime+"/current/prot/AAGraphScorer.java",
			runtime+"/current/prot/GlocalAminoLinear.java",runtime+"/current/prot/Blosum62.java",runtime+"/current/prot/"+producer+".java",runtime+"/current/prot/"+producer+".class",
			out.resolve("roster.tsv").toString(),out.resolve("consensus.faa").toString(),corpus,membership,out.resolve("recipe.txt").toString(),manifest};
	}
	static void requireHash(String file,String expected){
		DigestSuffix.requireSuffix(expected,"input sha80");
		if(!DigestSuffix.file(file).equals(expected)){throw new IllegalArgumentException("Pilot input changed or has wrong sha80: "+file);}
	}
	/** Binary digests are required by the existing bundle format; only sha80 text is emitted. */
	private static byte[] digest(String file) throws Exception{
		final MessageDigest md=MessageDigest.getInstance("SHA-256");final byte[] buffer=new byte[1<<16];
		try(InputStream in=Files.newInputStream(Paths.get(file))){for(int n=in.read(buffer); n>=0; n=in.read(buffer)){md.update(buffer,0,n);}}
		return md.digest();
	}
}
