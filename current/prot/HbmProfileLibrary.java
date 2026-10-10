package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;

import fileIO.ByteStreamWriter;
import structures.ByteBuilder;

/** Offline sharded refinement and strict composition of complete, explicitly named model sets. */
public final class HbmProfileLibrary {
	public static void main(String[] args){
		try{
			shared.Shared.setThreads(1);
			final HashMap<String,String> o=HmmComparisonData.options(args);
			for(String key : o.keySet()){
				if(!Arrays.asList("mode","resources","manifest","manifestsha80","runtime","out","ranks","shards","source","families","a","b").contains(key)){
					throw new IllegalArgumentException("Unknown library parameter: "+key);
				}
			}
			final String mode=HmmComparisonData.required(o,"mode");
			if(mode.equals("compare")){compare(Paths.get(HmmComparisonData.required(o,"a")),Paths.get(HmmComparisonData.required(o,"b")));return;}
			final HbmProfilePilot.Session s=new HbmProfilePilot.Session(o,!mode.equals("plan"));
			if(mode.equals("plan")){plan(o,s);}else if(mode.equals("shard")){shard(o,s);}else if(mode.equals("combine")){combine(o,s);}
			else{throw new IllegalArgumentException("Unknown library mode: "+mode);}
		}catch(Throwable failure){failure.printStackTrace();System.exit(1);}
	}
	/** Deterministic longest-estimated-work-first assignment; the cost is a proxy, not measured runtime. */
	private static void plan(HashMap<String,String> o,HbmProfilePilot.Session s) throws Exception{
		final Path out=Paths.get(HmmComparisonData.required(o,"out"));
		final int n=Integer.parseInt(o.getOrDefault("shards","32"));
		final int[] ranks=readRanks(o.get("ranks"),s.roster.size());
		if(n<1 || n>ranks.length){throw new IllegalArgumentException("Shard count must be within the selected family count");}
		final ArrayList<Work> work=new ArrayList<Work>();
		for(int rank : ranks){
			final long residues=Long.parseLong(s.rows.get(rank+1)[4]);
			work.add(new Work(rank,Math.multiplyExact(residues,s.sequence.get(s.roster.get(rank)).length)));
		}
		work.sort(Comparator.comparingLong((Work w)->w.cells).reversed().thenComparingInt(w->w.rank));
		Files.createDirectory(out);Files.createDirectory(out.resolve("shards"));
		final ByteStreamWriter[] writers=new ByteStreamWriter[n];final long[] cells=new long[n];final int[] counts=new int[n];
		for(int i=0; i<n; i++){writers[i]=HmmComparisonData.writer(out.resolve(String.format(Locale.ROOT,"shards/%02d.txt",i)).toString());}
		for(Work w : work){
			int best=0;for(int i=1; i<n; i++){if(cells[i]<cells[best]){best=i;}}
			writers[best].println(Integer.toString(w.rank));cells[best]=Math.addExact(cells[best],w.cells);counts[best]++;
		}
		final ByteStreamWriter summary=HmmComparisonData.writer(out.resolve("shards.tsv").toString());
		summary.println("shard\tfamilies\testimated_cells");
		for(int i=0; i<n; i++){HmmComparisonData.close(writers[i]);summary.println(new ByteBuilder().append(i).tab().append(counts[i]).tab().append(cells[i]));}
		HmmComparisonData.close(summary);s.verify();
		System.err.println("PROFILE_LIBRARY_PLAN_PASS families="+ranks.length+" shards="+n);
	}
	private static final class Work{
		Work(int rank_,long cells_){rank=rank_;cells=cells_;}
		final int rank;final long cells;
	}
	static int[] readRanks(String file,int count){
		if(count<1){throw new IllegalArgumentException("Empty source roster");}
		if(file==null){final int[] ranks=new int[count];for(int i=0; i<count; i++){ranks[i]=i;}return ranks;}
		final List<String[]> rows=HmmComparisonData.rows(file);final int[] ranks=new int[rows.size()];final boolean[] seen=new boolean[count];
		if(rows.isEmpty()){throw new IllegalArgumentException("Empty rank selection");}
		for(int i=0; i<ranks.length; i++){
			if(rows.get(i).length!=1){throw new IllegalArgumentException("Rank selection needs one integer per line");}
			final int rank=Integer.parseInt(rows.get(i)[0]);
			if(rank<0 || rank>=count || seen[rank]){throw new IllegalArgumentException("Invalid or duplicate selected family rank: "+rank);}
			seen[rank]=true;ranks[i]=rank;
		}
		return ranks;
	}
	private static Path family(Path root,int rank){return root.resolve(String.format(Locale.ROOT,"rank_%04d",rank));}
	private static void shard(HashMap<String,String> o,HbmProfilePilot.Session s) throws Exception{
		final Path root=Paths.get(HmmComparisonData.required(o,"out"));Files.createDirectories(root);
		final String source=HmmComparisonData.required(o,"source");
		final int[] ranks=readRanks(HmmComparisonData.required(o,"ranks"),s.roster.size());
		for(int rank : ranks){
			final Path dir=family(root,rank);
			if(Files.exists(dir)){validateFamily(dir,rank,s,source);continue;}
			final HashMap<String,String> one=new HashMap<String,String>();
			one.put("in",Paths.get(source,s.rows.get(rank+1)[5]).toString());one.put("rank",Integer.toString(rank));one.put("out",dir.toString());one.put("padding","20");
			HbmProfilePilot.runFamily(one,s,false);
		}
		s.verify();System.err.println("PROFILE_LIBRARY_SHARD_PASS families="+ranks.length);
	}
	/** Validates actual source bytes, output pins, producer provenance and all family identity/count bindings. */
	static HbmBundleLoader.Loaded validateFamily(Path dir,int rank,HbmProfilePilot.Session s,String source) throws Exception{
		if(!new String(Files.readAllBytes(dir.resolve("PASS")),StandardCharsets.UTF_8).trim().equals("PROFILE_PILOT_PASS")){throw new IllegalArgumentException("Incomplete refined family: "+dir);}
		final String[] original=s.rows.get(rank+1);final String rep=original[2];
		final List<String[]> order=HmmComparisonData.rows(dir.resolve("roster.tsv").toString());
		if(order.size()!=2 || order.get(1).length!=3 || !order.get(1)[0].equals("0") || !order.get(1)[1].equals(rep) || !order.get(1)[2].equals(original[3])){
			throw new IllegalArgumentException("Refined singleton roster/count differs: "+dir);
		}
		final List<String[]> summary=HmmComparisonData.rows(dir.resolve("summary.tsv").toString());
		if(summary.size()!=2 || summary.get(1).length!=12){throw new IllegalArgumentException("Refined family summary shape differs: "+dir);}
		final String[] row=summary.get(1);
		if(Integer.parseInt(row[0])!=rank || !row[1].equals(rep) || !row[2].equals(original[3]) || !row[3].equals(original[4])
			|| Integer.parseInt(row[4])!=s.sequence.get(rep).length || Long.parseLong(row[9])<0 || Long.parseLong(row[9])>Long.parseLong(row[3])){
			throw new IllegalArgumentException("Refined family member/count/identity binding differs: "+dir);
		}
		HbmProfilePilot.requireHash(dir.resolve("refined.mqhb").toString(),row[10]);
		HbmProfilePilot.requireHash(dir.resolve("consensus.faa").toString(),row[11]);
		final List<ProteinSequence> consensus=ProteinSearch.readFasta(dir.resolve("consensus.faa").toString());
		if(consensus.size()!=1 || !consensus.get(0).id.equals(rep) || consensus.get(0).length()!=Integer.parseInt(row[5])){throw new IllegalArgumentException("Refined consensus binding differs: "+dir);}
		if(!row[6].equals(Boolean.toString(!Arrays.equals(s.sequence.get(rep),consensus.get(0).enc)))){throw new IllegalArgumentException("Consensus change flag differs: "+dir);}
		final HashMap<String,String> recipe=recipe(dir.resolve("recipe.txt"));
		check(recipe,"parent_hbm_sha80",s.parentPin);check(recipe,"input_sha80",original[6]);check(recipe,"source_rank",Integer.toString(rank));check(recipe,"source_rep_id",rep);
		check(recipe,"beta","0.01");check(recipe,"clipmin","-4");check(recipe,"clipmax","11");check(recipe,"gap","4");check(recipe,"passes","2");check(recipe,"trimdepth","0.1");
		check(recipe,"growth_padding",row[7]);check(recipe,"final_padding","0");check(recipe,"background","parent_library_ref_counts_plus_one");check(recipe,"final_terminal_policy","count_exclusions");
		check(recipe,"runtime_source_manifest_sha80",DigestSuffix.file(s.runtime+"/POSITION_SOURCE.tsv"));
		check(recipe,"background_sha80",DigestSuffix.file(dir.resolve("background.tsv").toString()));
		if(!Arrays.equals(Files.readAllBytes(dir.resolve("background.tsv")),backgroundBytes(s.background))){throw new IllegalArgumentException("Refined family changed its frozen background: "+dir);}
		final String input=Paths.get(source,original[5]).toString();HbmProfilePilot.requireHash(input,original[6]);
		HbmProfilePilot.verifyProvenance(dir,s.runtime,input,input,s.manifest,"HbmProfilePilot");
		return HbmBundleLoader.load(dir.resolve("refined.mqhb"),Collections.singletonList(rep),id->consensus.get(0).enc,
			HbmBundleLoader.loadSemanticProvenance(dir.resolve("provenance.tsv").toString()));
	}
	private static HashMap<String,String> recipe(Path path) throws Exception{
		final List<String> rows=Files.readAllLines(path,StandardCharsets.UTF_8);final HashMap<String,String> map=new HashMap<String,String>();
		if(rows.isEmpty() || !rows.get(0).equals("experimental_profile_refinement_v1")){throw new IllegalArgumentException("Wrong family recipe schema: "+path);}
		for(int i=1; i<rows.size(); i++){
			final String line=rows.get(i);final int eq=line.indexOf('=');
			if(eq<1 || map.put(line.substring(0,eq),line.substring(eq+1))!=null){throw new IllegalArgumentException("Malformed/duplicate recipe field: "+path);}
		}
		return map;
	}
	private static void check(HashMap<String,String> map,String key,String expected){if(!expected.equals(map.get(key))){throw new IllegalArgumentException("Refined recipe differs at "+key);}}
	private static byte[] backgroundBytes(double[] background){
		final ByteBuilder b=new ByteBuilder().append("residue_code\tprobability\n");
		for(int a=0; a<20; a++){b.append(a).tab().append(Double.toString(background[a])).nl();}
		return b.toBytes();
	}
	private static void combine(HashMap<String,String> o,HbmProfilePilot.Session s) throws Exception{
		final Path root=Paths.get(HmmComparisonData.required(o,"families")),out=Paths.get(HmmComparisonData.required(o,"out"));
		if(Files.exists(out)){throw new IllegalArgumentException("Combined output must be a fresh directory");}
		final String source=HmmComparisonData.required(o,"source");
		final int[] ranks=readRanks(o.get("ranks"),s.roster.size());Arrays.sort(ranks);
		final ArrayList<HbmBundleLoader.Loaded> models=new ArrayList<HbmBundleLoader.Loaded>();
		final ArrayList<String> roster=new ArrayList<String>();final HashMap<String,byte[]> sequences=new HashMap<String,byte[]>();
		long members=0,residues=0;int changed=0;
		for(int rank : ranks){
			final Path dir=family(root,rank);models.add(validateFamily(dir,rank,s,source));
			final String rep=s.roster.get(rank);roster.add(rep);
			final ProteinSequence consensus=ProteinSearch.readFasta(dir.resolve("consensus.faa").toString()).get(0);sequences.put(rep,consensus.enc);
			members=Math.addExact(members,Long.parseLong(s.rows.get(rank+1)[3]));residues=Math.addExact(residues,Long.parseLong(s.rows.get(rank+1)[4]));
			if(!Arrays.equals(consensus.enc,s.sequence.get(rep))){changed++;}
		}
		final HbmBundleLoader.Loaded combined=HbmBundleLoader.Loaded.combine(models,roster);s.verify();
		Files.createDirectory(out);
		final ByteStreamWriter refs=HmmComparisonData.writer(out.resolve("consensus.faa").toString()),order=HmmComparisonData.writer(out.resolve("roster.tsv").toString());
		final ByteStreamWriter artifacts=HmmComparisonData.writer(out.resolve("family_artifacts.tsv").toString());
		order.println("rank\trep_id\tmembers");artifacts.println("source_rank\trep_id\tmembers\tresidues\tinput_sha80\tbundle_sha80\tconsensus_sha80");
		for(int i=0; i<ranks.length; i++){
			final int rank=ranks[i];final String[] entry=s.rows.get(rank+1);final Path dir=family(root,rank);
			HmmComparisonData.writeFasta(refs,entry[2],sequences.get(entry[2]));
			order.println(new ByteBuilder().append(i).tab().append(entry[2]).tab().append(entry[3]));
			artifacts.println(new ByteBuilder().append(rank).tab().append(entry[2]).tab().append(entry[3]).tab().append(entry[4]).tab().append(entry[6])
				.tab().append(DigestSuffix.file(dir.resolve("refined.mqhb").toString())).tab().append(DigestSuffix.file(dir.resolve("consensus.faa").toString())));
		}
		HmmComparisonData.close(refs);HmmComparisonData.close(order);HmmComparisonData.close(artifacts);
		Files.write(out.resolve("background.tsv"),backgroundBytes(s.background));
		final ByteStreamWriter recipe=HmmComparisonData.writer(out.resolve("recipe.txt").toString());
		recipe.println("experimental_refined_library_v1\nbeta=0.01\nclipmin=-4\nclipmax=11\ngap=4\npasses=2\nbackground=parent_library_ref_counts_plus_one");
		recipe.println("parent_hbm_sha80="+s.parentPin+"\nmember_manifest_sha80="+s.manifestPin+"\nruntime_source_manifest_sha80="+DigestSuffix.file(s.runtime+"/POSITION_SOURCE.tsv"));
		recipe.println("family_artifacts_sha80="+DigestSuffix.file(out.resolve("family_artifacts.tsv").toString())+"\nbackground_sha80="+DigestSuffix.file(out.resolve("background.tsv").toString()));
		HmmComparisonData.close(recipe);
		final byte[][] provenance=HbmProfilePilot.writeProvenance(out,s.runtime,s.manifest,out.resolve("family_artifacts.tsv").toString(),s.manifest,"HbmProfileLibrary");
		combined.writeBundle(out.resolve("refined.mqhb"),provenance);
		final HbmBundleLoader.Loaded reloaded=HbmBundleLoader.load(out.resolve("refined.mqhb"),roster,id->sequences.get(id),
			HbmBundleLoader.loadSemanticProvenance(out.resolve("provenance.tsv").toString()));
		reloaded.assertStructuralEquivalent(combined);
		final ByteStreamWriter summary=HmmComparisonData.writer(out.resolve("summary.tsv").toString());
		summary.println("families\tmembers\tresidues\tchanged_consensuses\tbundle_sha80\tconsensus_sha80\tbackground_sha80");
		summary.println(new ByteBuilder().append(ranks.length).tab().append(members).tab().append(residues).tab().append(changed)
			.tab().append(DigestSuffix.file(out.resolve("refined.mqhb").toString())).tab().append(DigestSuffix.file(out.resolve("consensus.faa").toString()))
			.tab().append(DigestSuffix.file(out.resolve("background.tsv").toString())));HmmComparisonData.close(summary);
		final ByteStreamWriter pass=HmmComparisonData.writer(out.resolve("PASS").toString());pass.println("PROFILE_LIBRARY_COMBINE_PASS");HmmComparisonData.close(pass);
		System.err.println("PROFILE_LIBRARY_COMBINE_PASS families="+ranks.length+" members="+members+" residues="+residues);
	}
	private static void compare(Path a,Path b) throws Exception{
		final List<ProteinSequence> left=ProteinSearch.readFasta(a.resolve("consensus.faa").toString()),right=ProteinSearch.readFasta(b.resolve("consensus.faa").toString());
		if(left.size()!=1 || right.size()!=1 || !left.get(0).id.equals(right.get(0).id) || !Arrays.equals(left.get(0).enc,right.get(0).enc)){throw new IllegalArgumentException("Compared singleton consensuses differ");}
		final String id=left.get(0).id;final byte[] consensus=left.get(0).enc;
		final HbmBundleLoader.Loaded x=HbmBundleLoader.load(a.resolve("refined.mqhb"),Collections.singletonList(id),name->consensus,HbmBundleLoader.loadSemanticProvenance(a.resolve("provenance.tsv").toString()));
		final HbmBundleLoader.Loaded y=HbmBundleLoader.load(b.resolve("refined.mqhb"),Collections.singletonList(id),name->consensus,HbmBundleLoader.loadSemanticProvenance(b.resolve("provenance.tsv").toString()));
		x.assertStructuralEquivalent(y);System.err.println("PROFILE_LIBRARY_FAMILY_PARITY_PASS");
	}
}
