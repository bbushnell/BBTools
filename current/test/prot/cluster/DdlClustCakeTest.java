package test.prot.cluster;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import ddl.FamilySketchBuilder;
import ddl.ProteinKmerGeometry;
import ddl.ProteinPositionMode;
import fileIO.ByteFile;
import fileIO.ReadWrite;
import prot.ProteinSequence;
import prot.cluster.ClusterSketchCodec;
import prot.cluster.DdlClust;
import prot.cluster.DdlClustCake;
import prot.cluster.DdlClustStageB;
import prot.cluster.Family;

/** Complete tiny-arm oracle: unchanged setup, empty-template retention, explicit
 * lambda changes, and final top-budget reassignment. Run only on Dori.
 * @author Brian Bushnell, Keqing
 */
public final class DdlClustCakeTest{

	public static void main(final String[] args)throws Exception{
		check(args.length==1,"Supply one fresh test output directory");
		final Path root=Paths.get(args[0]);check(!Files.exists(root),"Do not overwrite test evidence");
		Files.createDirectories(root);
		PhaseCostLogTest.main(new String[]{root.resolve("phase_cost_test").toString()});
		final ProteinSequence[] sequences={sequence('A'),sequence('W'),sequence('C')};
		final FamilySketchBuilder builder=new FamilySketchBuilder(64,2,7,false,0,0,null);
		final Family[] families=new Family[3];
		for(int i=0; i<families.length; i++){
			families[i]=new Family(7+i*12,"source"+i,i,sequences[i].enc);
			families[i].refreshStats(sequences,builder);
		}
		final Path codec=root.resolve("source.pddl");
		ClusterSketchCodec.write(codec.toString(),false,families,2,64,7,false,0,0);
		final Path input=root.resolve("tails.faa");
		Files.write(input,(">a0\n"+repeat('A',100)+"\n>a1\n"+repeat('A',100)+"\n>a2\n"+repeat('A',100)+
			"\n>w\n"+repeat('W',100)+"\n>short\n"+repeat('C',90)+"\n").getBytes(StandardCharsets.US_ASCII));
		for(int threads:new int[]{1,4}){
			final Path result=root.resolve("arm_t"+threads);
			final ArrayList<String> argv=arguments(result,input,codec,threads);
			argv.add("fractions=.25,1");argv.add("sizelambdas=.25,0");
			DdlClust.main(argv.toArray(new String[argv.size()]));
			for(String recipe:new String[]{"run.config","layer1/run.config","layer2/run.config"}){
				check(setting(result.resolve(recipe),"memberranking").equals("false"),
					"Explicit legacy mode survives default changes and is recorded once for replay");
				check(setting(result.resolve(recipe),"memberhydroweight").equals("0.2"),
					"Default coefficient is recorded once even when legacy ranking leaves it unused");
			}
			PhaseCostLogTest.verifyReport(result.resolve("phase_costs.tsv"),"cake_process","seed_selection",
				"initial_assignment","layer_1","fission_waves","fusion",
				"fusion_candidates","family_index","nearest_neighbor_comparison","selection_publication","final_refinement");
			final List<String> first=Files.readAllLines(result.resolve("layer1/state.tsv"),StandardCharsets.US_ASCII);
			final List<String> second=Files.readAllLines(result.resolve("layer2/state.tsv"),StandardCharsets.US_ASCII);
			check(first.contains("size_lambda\t0.25") && second.contains("size_lambda\t0.0"),"Each layer receives its own strength, including zero");
			check(first.contains("slot_fraction\t0.25") && second.contains("slot_fraction\t1.0"),"Each layer receives its exact slot fraction");
			check(first.contains("setup_slot_cost\t3") && second.contains("setup_slot_cost\t3"),"Both layers retain setup price from actual initial counts");
			check(first.contains("below_mode\tpriorscaled") && second.contains("below_mode\tpriorscaled"),
				"Cake defaults to prior-scaled floor pricing");
			check(Files.readAllLines(result.resolve("run.config")).contains("setupslotcost=auto"),"Default setup cost remains automatic and is recorded");
			checkNoIntermediateModels(result,2);
			check(second.contains("models\t3") && second.contains("active_models\t2") && second.contains("assigned\t4"),
				"No-op resident layers preserve all three templates, including the empty reference");
			final ClusterSketchCodec.Loaded finalModels=ClusterSketchCodec.read(result.resolve("final/families.pddl.gz").toString());
			check(finalModels.families.length==1 && finalModels.families[0].familyId==0 && finalModels.families[0].members==3,
				"Final top-one selection retains only the three A tails");
			final List<String> table=Files.readAllLines(result.resolve("layers.tsv"),StandardCharsets.US_ASCII);
			check(table.size()==5 && table.get(1).equals("initial\t0\tNA\tNA\t3\t2\t4\t5\t3\t3") &&
				table.get(4).equals("final\t2\tNA\tNA\t1\t1\t3\t5\t3\t3"),"Coverage bracket uses all five input records, including the short tail");
			check(Files.exists(result.resolve("final/quality.cohesion.tsv")),"Standard quality output is produced only by final refinement");
			check(Files.readAllLines(result.resolve("final/stats.tsv"),StandardCharsets.US_ASCII)
				.contains("completed_passes\t2"),"Stable membership stops after pass2 despite the configured maximum3");
			checkQuality(result.resolve("final/quality.cohesion.tsv"));
			check(Files.readAllLines(result.resolve("final/family_sizes.tsv"),StandardCharsets.US_ASCII)
				.equals(Arrays.asList("size\tfamilies\tassigned_members","3\t1\t3")),"Family sizes are actual final member counts, not shortlist lengths");
			check(!Files.exists(result.resolve("initial/quality.cohesion.tsv")),"Initial assignment does not do quality work");
		}
		final Path control=root.resolve("control");
		final ArrayList<String> argv=arguments(control,input,codec,1);argv.add("fractions=none");argv.add("lambda=0");
		argv.add("stopchangefraction=0");
		final Path config=root.resolve("control.config");
		Files.write(config,argv,StandardCharsets.US_ASCII);
		DdlClust.main(new String[]{"config="+config});
		check(!hasSetting(control.resolve("run.config"),"bgzip") &&
			!hasSetting(control.resolve("run.config"),"nativebgzip"),"Absent backend options do not change historical recorded defaults");
		check(!Files.exists(control.resolve("layer1")) && Files.readAllLines(control.resolve("layers.tsv")).size()==3,
			"No-cake control bypasses every layer but keeps setup and final refinement");
		check(Files.readAllLines(control.resolve("final/stats.tsv"),StandardCharsets.US_ASCII)
			.contains("completed_passes\t3"),"Disabling early stopping retains all three configured passes");
		checkQuality(control.resolve("final/quality.cohesion.tsv"));
		check(Arrays.equals(Files.readAllBytes(control.resolve("final/quality.cohesion.tsv")),
			Files.readAllBytes(root.resolve("arm_t1/final/quality.cohesion.tsv"))),
			"Early-stop and fixed-length runs publish identical final cohesion for identical members/models");
		check(DdlClustCake.lambdaAt(new double[]{.125},4,5)==.125,"Singleton strength is constant at the last layer");
		check(DdlClustCake.priorSlotCost(0,79999,80000,0,1.01)==0 &&
			DdlClustCake.priorSlotCost(29,79999,80000,0,1.01)==29 &&
			DdlClustCake.priorSlotCost(50,90000,80000,39,1.01)==39 &&
			DdlClustCake.priorSlotCost(39,79999,80000,0,1.01)==39 &&
			DdlClustCake.priorSlotCost(39,90000,80000,50,1.01)==50 &&
			DdlClustCake.priorSlotCost(39,82000,80000,19,1.01)==19 &&
			DdlClustCake.priorSlotCost(19,80800,80000,21,1.01)==19 &&
			DdlClustCake.priorSlotCost(19,79007,80000,0,1.01)==19,
			"Prior pricing must handle initial underfill, reentry, above-edge replacement and edge equality");
		for(int override:new int[]{0,29}){
			final Path result=root.resolve("fixed_setup_"+override);
			final ArrayList<String> fixed=arguments(result,input,codec,override==0 ? 1 : 4);
			fixed.add("familybudget=3");fixed.add("initialfraction=1");
			fixed.add("fractions=.25,1");fixed.add("setupslotcost="+override);fixed.add("belowmode=constant");
			DdlClust.main(fixed.toArray(new String[fixed.size()]));
			for(int layer=1; layer<=2; layer++){
				check(Files.readAllLines(result.resolve("layer"+layer+"/state.tsv"))
					.contains("setup_slot_cost\t"+override),"Explicit setup cost reaches every layer unchanged");
				check(Files.readAllLines(result.resolve("layer"+layer+"/state.tsv"))
					.contains("initial_slot_cost\t"+override),"An underfilled roster uses the explicit fallback as its effective cost");
			}
			check(Files.readAllLines(result.resolve("run.config")).contains("setupslotcost="+override),"Explicit setup cost is recorded for replay");
			final List<String> rows=Files.readAllLines(result.resolve("layers.tsv"));
			for(int i=1; i<rows.size(); i++){
				check(rows.get(i).split("\t",-1)[8].equals(Integer.toString(override)),"Initial, layer and final summaries retain the explicit cost");
			}
			check(rows.get(1).startsWith("initial\t0\tNA\tNA\t3\t2\t4\t5\t"),
				"The cost override does not change initial assignment or empty-template retention");
			checkNoIntermediateModels(result,2);
		}
		final Path explicitPrior=root.resolve("explicit_priorscaled");
		final ArrayList<String> prior=arguments(explicitPrior,input,codec,1);
		prior.add("familybudget=3");prior.add("initialfraction=1");
		prior.add("fractions=.25,1");prior.add("setupslotcost=29");
		DdlClust.main(prior.toArray(new String[prior.size()]));
		for(int layer=1; layer<=2; layer++){
			final List<String> state=Files.readAllLines(explicitPrior.resolve("layer"+layer+"/state.tsv"));
			check(state.contains("setup_slot_cost\t29") && state.contains("initial_slot_cost\t19.333333333333332"),
				"Explicit prior-scaled underfill must carry the unscaled fallback and scale it by active/budget");
		}
		final Path leftPinned=root.resolve("left_pinned");
		final ArrayList<String> left=arguments(leftPinned,input,codec,1);
		left.add("anchorterminus=head");left.add("fractions=.25");
		DdlClust.main(left.toArray(new String[left.size()]));
		check(setting(leftPinned.resolve("run.config"),"anchorterminus").equals("head"),
			"Cake accepts and records explicit left-pinned head clustering geometry");
		for(String recipe:new String[]{"initial/stats.tsv","final/stats.tsv"}){
			check(tabSetting(leftPinned.resolve(recipe),"pinleft").equals("true") &&
				tabSetting(leftPinned.resolve(recipe),"pinright").equals("false"),
				"Initial and final DdlClust calls record normalized head pin flags");
		}
		check(setting(leftPinned.resolve("layer1/run.config"),"pinleft").equals("true") &&
			setting(leftPinned.resolve("layer1/run.config"),"pinright").equals("false"),
			"Cake forwards normalized head pin flags to resident DdlClust and StageB calls");
		final Path replay=root.resolve("left_pinned_replay");
		final ArrayList<String> replayConfig=new ArrayList<String>();
		for(String line:Files.readAllLines(leftPinned.resolve("run.config"),StandardCharsets.US_ASCII)){
			replayConfig.add(line.startsWith("outdir=") ? "outdir="+replay : line);
		}
		final Path replayPath=root.resolve("left_pinned_replay.config");
		Files.write(replayPath,replayConfig,StandardCharsets.US_ASCII);
		DdlClust.main(new String[]{"config="+replayPath});
		for(String recipe:new String[]{"initial/stats.tsv","final/stats.tsv"}){
			check(tabSetting(replay.resolve(recipe),"pinleft").equals("true") &&
				tabSetting(replay.resolve(recipe),"pinright").equals("false"),
				"Recorded head Cake config replays after normalized flags replace contradictory defaults");
		}
		for(String invalid:new String[]{"pinleft=t","pinleft=f","pinright=t","pinright=f"}){
			final Path bad=root.resolve("bad_anchor_"+checks);
			final ArrayList<String> a=arguments(bad,input,codec,1);
			a.add("anchorterminus="+(invalid.equals("pinleft=t") ? "tail" :
				invalid.equals("pinright=f") ? "tail" : "head"));
			a.add(invalid);
			checkRejected(a,bad,"Contradictory anchorterminus and legacy pin flags must fail");
		}
		checkStageBAnchorRejection();
		for(String invalid:new String[]{"sizelambdas=.25,","sizelambdas=.25,.2,.1","sizelambdas=-1","fractions=-1","sizelambdas=Infinity",
			"setupslotcost=-1","setupslotcost=null","setupslotcost=NaN","setupslotcost=2147483648","setupslotcost=Infinity","memberhydroweight=0"}){
			final Path bad=root.resolve("bad"+checks);
			final ArrayList<String> a=arguments(bad,input,codec,1);a.add("fractions=.25,1");a.add(invalid);
			boolean rejected=false;
			try{DdlClust.main(a.toArray(new String[a.size()]));}catch(RuntimeException expected){rejected=true;}
			check(rejected && !Files.exists(bad),"Invalid schedule fails before output creation: "+invalid);
		}
		checkFissionSchedules(root,input,codec);
		checkPositionModes(root);
		checkFractionalGeometry(root);
		checkOrderedTailInputs(root,codec);
		System.out.println("DDLCLUST_CAKE_TEST_PASS checks="+checks);
	}

	/** Exercises the same fused input path through greedy and all resident Cake stages.
	 * Literal pre-cropped proteins are the oracle; no crop utility generates them.
	 */
	private static void checkOrderedTailInputs(final Path root,final Path codec)throws Exception{
		final Path first=root.resolve("crop_first.faa"),second=root.resolve("crop_second.faa");
		final Path reference=root.resolve("crop_reference.faa");
		final String firstRaw=">dup first\n"+repeat('W',70)+repeat('A',100)+"\n>dup second\n"+repeat('A',100)+"\n";
		final String secondRaw=">a2 third\n"+repeat('C',50)+repeat('A',100)+"\n>w fourth\n"+repeat('A',100)+repeat('W',100)+"\n";
		final String referenceRaw=">dup first\n"+repeat('A',100)+"\n>dup second\n"+repeat('A',100)+
			"\n>a2 third\n"+repeat('A',100)+"\n>w fourth\n"+repeat('W',100)+"\n";
		Files.write(first,firstRaw.getBytes(StandardCharsets.US_ASCII));
		Files.write(second,secondRaw.getBytes(StandardCharsets.US_ASCII));
		Files.write(reference,referenceRaw.getBytes(StandardCharsets.US_ASCII));
		final String sources=first+","+second;
		Path greedyOracle=null,cakeOracle=null;
		for(boolean fused:new boolean[]{false,true}){
			final Path greedy=Files.createDirectory(root.resolve("crop_greedy_"+fused));
			final ArrayList<String> greedyArgs=new ArrayList<String>(Arrays.asList(
				"in="+(fused ? sources : reference),"tailcrop="+(fused ? 100 : 0),"memberranking=f",
				"members="+greedy.resolve("members.tsv.gz"),"ddlout="+greedy.resolve("families.pddl.gz"),
				"out="+greedy.resolve("family#.faa"),"stats="+greedy.resolve("stats.tsv"),"refinement=null",
				"passes=0","mincluster=1","topn=0","t=1","seed=1701","k=2","width=64","sketchseed=7",
				"positional=f","minlen=100","minspan=70","maxalignlen=100","minid=1","minblosumratio=0",
				"greedyhits=0","maxhydrodelta=9","greedymaxhydrodelta=9","ow=f"));
			DdlClust.main(greedyArgs.toArray(new String[greedyArgs.size()]));
			final Path cake=root.resolve("crop_cake_"+fused);
			final ArrayList<String> cakeArgs=arguments(cake,reference,codec,1);
			cakeArgs.addAll(Arrays.asList("in="+(fused ? sources : reference),"tailcrop="+(fused ? 100 : 0),
				"fractions=.25,1","sizelambdas=.25,0","quality=null"));
			if(fused){cakeArgs.addAll(Arrays.asList("bgzip=t","nativebgzip=f"));}
			DdlClust.main(cakeArgs.toArray(new String[cakeArgs.size()]));
			if(!fused){greedyOracle=greedy;cakeOracle=cake;}
			else{
				for(String name:new String[]{"families.pddl.gz","members.tsv.gz"}){
					checkFileParity(greedyOracle.resolve(name),greedy.resolve(name));
					checkFileParity(cakeOracle.resolve("final").resolve(name),cake.resolve("final").resolve(name));
				}
				for(ClusterSketchCodec.Entry family:ClusterSketchCodec.read(greedy.resolve("families.pddl.gz").toString()).families){
					checkFileParity(greedyOracle.resolve("family"+family.familyId+".faa"),greedy.resolve("family"+family.familyId+".faa"));
				}
				checkFileParity(cakeOracle.resolve("layers.tsv"),cake.resolve("layers.tsv"));
				check(tabSetting(greedy.resolve("stats.tsv"),"tailcrop").equals("100"),"Greedy records the fused crop");
				check(setting(cake.resolve("run.config"),"tailcrop").equals("100") &&
					setting(cake.resolve("run.config"),"in").equals(sources),"Cake records the ordered original inputs and crop for replay");
				for(String stage:new String[]{"initial","final"}){
					check(tabSetting(cake.resolve(stage+"/stats.tsv"),"tailcrop").equals("100"),"Both assignments retain the population crop contract");
				}
				check(!hasSetting(cake.resolve("layer1/run.config"),"tailcrop"),"Resident StageB borrows cropped inputs without an unsupported crop argument");
				checkNoIntermediateModels(cake,2);
				check(!Files.exists(cake.resolve("population.faa")),"Final Cake does not materialize an input bridge");
				checkCompressionReplay(root,cake);
			}
		}
		final Path bad=root.resolve("bad_crop_setting");
		final ArrayList<String> invalid=arguments(bad,reference,codec,1);
		invalid.addAll(Arrays.asList("fractions=none","tailcrop=-1"));
		checkRejected(invalid,bad,"Negative crop lengths fail before Cake output creation");
		final Path alias=Files.createDirectory(root.resolve("crop_input_alias"));
		Files.createLink(alias.resolve("members.tsv"),second);
		boolean rejected=false;
		try{
			new DdlClust(new String[]{"in="+sources,"members="+alias.resolve("members.tsv"),
				"stats=null","refinement=null","ow=t","memberranking=f"});
		}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected,"Greedy output cannot overwrite a hard link to the second input");
		check(new String(Files.readAllBytes(first),StandardCharsets.US_ASCII).equals(firstRaw) &&
			new String(Files.readAllBytes(second),StandardCharsets.US_ASCII).equals(secondRaw),"Both original inputs remain byte-identical");
	}

	/** Saved recipes must reestablish explicit backend policy in a fresh JVM's default state. */
	private static void checkCompressionReplay(final Path root,final Path original)throws Exception{
		for(String name:new String[]{"run.config","layer1/run.config","layer2/run.config"}){
			check(setting(original.resolve(name),"bgzip").equals("t") &&
				setting(original.resolve(name),"nativebgzip").equals("f"),"Explicit backend flags reach both saved top-level and child recipes");
		}
		final Path replay=root.resolve("crop_compression_replay");
		final ArrayList<String> lines=new ArrayList<String>();
		for(String line:Files.readAllLines(original.resolve("run.config"),StandardCharsets.US_ASCII)){
			lines.add(line.startsWith("outdir=") ? "outdir="+replay : line);
		}
		final Path config=root.resolve("crop_compression_replay.config");
		Files.write(config,lines,StandardCharsets.US_ASCII);
		final boolean allow=ReadWrite.ALLOW_NATIVE_BGZF,prefer=ReadWrite.PREFER_NATIVE_BGZF_OUT,use=ReadWrite.USE_BGZIP;
		try{
			ReadWrite.ALLOW_NATIVE_BGZF=true;ReadWrite.PREFER_NATIVE_BGZF_OUT=true;ReadWrite.USE_BGZIP=false;
			DdlClust.main(new String[]{"config="+config});
			check(!ReadWrite.ALLOW_NATIVE_BGZF && !ReadWrite.PREFER_NATIVE_BGZF_OUT && ReadWrite.USE_BGZIP,
				"Standalone saved-config replay restores explicit external policy without inherited process state");
		}finally{
			ReadWrite.ALLOW_NATIVE_BGZF=allow;ReadWrite.PREFER_NATIVE_BGZF_OUT=prefer;ReadWrite.USE_BGZIP=use;
		}
		for(String name:new String[]{"families.pddl.gz","members.tsv.gz"}){
			checkFileParity(original.resolve("final").resolve(name),replay.resolve("final").resolve(name));
		}
		checkFileParity(original.resolve("layers.tsv"),replay.resolve("layers.tsv"));
	}

	/** Cake must accept and preserve fractional amino-bit geometry at every resident
	 * DdlClust assignment boundary; the selected codec header catches any mismatch.
	 */
	private static void checkFractionalGeometry(final Path root)throws Exception{
		final ProteinKmerGeometry geometry=ProteinKmerGeometry.make(5,"hydro");
		final String high=motif(100,"IVLFCMAGTSG"), low=motif(100,"WYHPNQDEKRE");
		final ProteinSequence[] sequences={new ProteinSequence("h|a0",high),
			new ProteinSequence("h|a1",high),new ProteinSequence("h|b",low)};
		final Path input=root.resolve("fractional.faa");
		Files.write(input,(">h|a0\n"+high+"\n>h|a1\n"+high+"\n>h|b\n"+low+"\n").getBytes(StandardCharsets.US_ASCII));
		final FamilySketchBuilder builder=new FamilySketchBuilder(64,geometry,3,ProteinPositionMode.PLAIN,0,0,null);
		final Family[] families={new Family(7,"hi",0,sequences[0].enc),new Family(19,"lo",2,sequences[2].enc)};
		families[0].addMember(1);
		for(Family family:families){family.refreshStats(sequences,builder);}
		final Path codec=root.resolve("fractional.pddl");
		ClusterSketchCodec.write(codec.toString(),false,families,geometry,64,3,ProteinPositionMode.PLAIN,0,0,null);
		final Path out=root.resolve("fractional_hydro");
		final ArrayList<String> args=arguments(out,input,codec,1);
		args.remove("k=2");args.add("k=5");args.add("aabits=hydro");
		args.remove("sketchseed=7");args.add("sketchseed=3");
		args.remove("familybudget=1");args.add("familybudget=2");
		args.remove("initialfraction=3");args.add("initialfraction=1");
		args.add("fractions=.25,1");args.add("passes=1");args.add("stopchangefraction=0");args.add("quality=null");
		DdlClust.main(args.toArray(new String[args.size()]));
		final ClusterSketchCodec.Loaded actual=ClusterSketchCodec.read(out.resolve("final/families.pddl.gz").toString());
		check(actual.header.geometry.equals(geometry) && actual.families.length==2,
			"Cake accepts fractional amino geometry and retains it in the final codec");
		check(setting(out.resolve("run.config"),"aabits").equals("hydro"),"Cake records the fractional amino-bit setting");
		check(tabSetting(out.resolve("initial/stats.tsv"),"aabits").equals("hydro"),
			"Initial resident assignment receives fractional geometry");
		check(setting(out.resolve("layer1/run.config"),"selectedcount").equals("2"),
			"Fractional geometry reaches the first live-probe layer without parser rejection");
		final Path bad=root.resolve("bad_fractional_"+checks);
		final ArrayList<String> both=arguments(bad,input,codec,1);
		both.add("aabits=hydro");both.add("aapartition=IVLFCMAGTS");
		checkRejected(both,bad,"Cake rejects conflicting amino-bit aliases");
	}

	private static void checkStageBAnchorRejection(){
		for(String invalid:new String[]{"anchorterminus=head:pinleft=f","anchorterminus=head:pinright=t",
				"anchorterminus=tail:pinleft=t","anchorterminus=tail:pinright=f"}){
			final String[] args=invalid.split(":",-1);
			boolean rejected=false;
			try{new DdlClustStageB(args);}catch(RuntimeException expected){rejected=true;}
			check(rejected,"StageB rejects explicit anchor/pin contradictions: "+invalid);
		}
	}

	/** Real cake calls must preserve the chosen seed/query layout through both layers and refinement. */
	private static void checkPositionModes(final Path root)throws Exception{
		final String first=motif(100,"ACDEFGHIKLMNPQRSTVWY"),second=motif(149,"WVPITNQYFECMAHRSDKGL");
		final ProteinSequence[] sequences={new ProteinSequence("g|a0",first),new ProteinSequence("g|a1",first),
			new ProteinSequence("g|b",second)};
		final Path input=root.resolve("position_modes.faa");
		Files.write(input,(">g|a0\n"+first+"\n>g|a1\n"+first+"\n>g|b\n"+second+"\n").getBytes(StandardCharsets.US_ASCII));
		for(ProteinPositionMode mode:ProteinPositionMode.values()){
			final int partitions=mode.positional() ? 4 : 1,overlap=mode.positional() ? 2 : 0;
			final FamilySketchBuilder builder=new FamilySketchBuilder(1024,5,7,mode,partitions,overlap,null);
			final Family[] families={new Family(7,"first",0,sequences[0].enc),new Family(19,"second",2,sequences[2].enc)};
			families[0].addMember(1);
			for(Family family:families){family.refreshStats(sequences,builder);}
			final Path codec=root.resolve("position_"+mode.flag+".pddl");
			ClusterSketchCodec.write(codec.toString(),false,families,5,1024,7,mode,partitions,overlap,null);
			Path serial=null;
			for(int threads:new int[]{1,4}){
				final Path out=root.resolve("position_"+mode.flag+"_t"+threads);
				final ArrayList<String> args=positionArguments(out,input,codec,mode,threads);
				args.add("positionmode="+mode.flag);
				DdlClust.main(args.toArray(new String[args.size()]));
				final ClusterSketchCodec.Loaded actual=ClusterSketchCodec.read(out.resolve("final/families.pddl.gz").toString());
				check(actual.header.positionMode==mode && actual.families.length==2,"Cake keeps explicit layout and both distinct families");
				for(int i=0;i<2;i++){
					check(actual.families[i].familyId==i && actual.families[i].members==(i==0 ? 2 : 1),
						"Resident explicit-mode calls preserve the literal two-family membership oracle");
					check(Arrays.equals(actual.families[i].sketch,families[i].sketch()),"Final family sketch must match the explicitly constructed seed layout");
				}
				check(setting(out.resolve("run.config"),"positionmode").equals(mode.flag),"Saved recipe retains canonical explicit mode");
				check(setting(out.resolve("run.config"),"positional").equals(Boolean.toString(mode.positional())),
					"An omitted legacy boolean follows the explicit mode, including plain");
				checkNoIntermediateModels(out,2);
				if(serial==null){serial=out;}else{
					checkFileParity(serial.resolve("final/families.pddl.gz"),out.resolve("final/families.pddl.gz"));
					checkFileParity(serial.resolve("final/members.tsv.gz"),out.resolve("final/members.tsv.gz"));
				}
			}
			if(mode!=ProteinPositionMode.PLAIN){
				final Path out=root.resolve("position_"+mode.flag+"_compatibility");
				final ArrayList<String> args=positionArguments(out,input,codec,mode,1);
				if(mode==ProteinPositionMode.RELATIVE){args.add("positionalmode=relative");}
				DdlClust.main(args.toArray(new String[args.size()]));
				checkFileParity(serial.resolve("final/families.pddl.gz"),out.resolve("final/families.pddl.gz"));
				checkFileParity(serial.resolve("final/members.tsv.gz"),out.resolve("final/members.tsv.gz"));
			}
			for(String flag:new String[]{"positionmode=unknown","positionmode=multiscale","positionmode=null"}){
				final Path bad=root.resolve("bad_position_"+checks);
				final ArrayList<String> args=positionArguments(bad,input,codec,mode,1);args.add(flag);
				checkRejected(args,bad,"Unsupported coordinate mode must fail before output creation");
			}
			for(boolean booleanFirst:new boolean[]{false,true}){
				final Path bad=root.resolve("contradictory_position_"+checks);
				final ArrayList<String> args=positionArguments(bad,input,codec,mode,1);
				final String bool="positional="+!mode.positional(),explicit="positionmode="+mode.flag;
				args.add(booleanFirst ? bool : explicit);args.add(booleanFirst ? explicit : bool);
				checkRejected(args,bad,"Contradictory explicit mode and boolean must fail in either order");
			}
		}
	}

	private static ArrayList<String> positionArguments(final Path out,final Path input,final Path codec,
			final ProteinPositionMode mode,final int threads){
		final ArrayList<String> args=arguments(out,input,codec,threads);
		args.remove("positional=f");args.remove("partitions=0");args.remove("overlap=0");
		args.addAll(Arrays.asList("k=5","width=1024","maxalignlen=200","familybudget=2","initialfraction=1",
			"partitions="+(mode.positional()?4:1),"overlap="+(mode.positional()?2:0),"fractions=.25,1",
			"passes=1","stopchangefraction=0","quality=null"));
		return args;
	}
	private static String motif(final int length,final String text){
		check(length>0 && !text.isEmpty(),"Position fixture requires a nonempty literal motif");
		final StringBuilder result=new StringBuilder(length);
		for(int i=0;i<length;i++){result.append(text.charAt(i%text.length()));}
		return result.toString();
	}

	private static void checkFissionSchedules(final Path root,final Path input,final Path codec)throws Exception{
		final String[] variants={null,"fissiontopfraction=1","fissiontopfractions=1","fissiontopfractions=1,1,1",
			"fissiontopfraction=.5","fissiontopfractions=.5","fissiontopfraction=0","fissiontopfractions=0"};
		final Path[] results=new Path[variants.length];
		for(int v=0; v<variants.length; v++){
			final Path result=results[v]=root.resolve("fission_constant_"+v);
			final ArrayList<String> args=fissionArguments(result,input,codec);
			if(variants[v]!=null){args.add(variants[v]);}
			DdlClust.main(args.toArray(new String[args.size()]));
			final double expected=v<4 ? 1 : v<6 ? .5 : 0;
			for(int layer=1; layer<=3; layer++){
				check(Double.parseDouble(setting(result.resolve("layer"+layer+"/run.config"),"fissiontopfraction"))==expected,
					"Default, scalar and constant schedule reach every actual child invocation");
			}
			if(v>0 && v<4){checkArmParity(results[0],result);}
			if(v==5 || v==7){checkArmParity(results[v-1],result);}
		}
		final Path result=root.resolve("fission_first_only");
		final ArrayList<String> args=fissionArguments(result,input,codec);
		args.add("fissiontopfractions=1,0,0");
		final Path config=root.resolve("fission_first_only.config");
		Files.write(config,args,StandardCharsets.US_ASCII);
		DdlClust.main(new String[]{"config="+config});
		check(setting(result.resolve("run.config"),"fissiontopfractions").equals("1.0,0.0,0.0"),
			"The parent recipe preserves the complete fission schedule");
		PhaseCostLogTest.verifyReport(result.resolve("phase_costs.tsv"),"residual_index","proposal_setup","proposal_evaluation");
		for(String line:Files.readAllLines(result.resolve("run.config"),StandardCharsets.US_ASCII)){
			check(!line.startsWith("fissiontopfraction="),"Scheduled parent recipes must not contain a conflicting scalar");
		}
		// Constructor validation of the recorded recipe catches scalar/schedule replay conflicts.
		final ArrayList<String> replay=new ArrayList<String>();
		for(String line:Files.readAllLines(result.resolve("run.config"),StandardCharsets.US_ASCII)){
			// This low-level constructor receives parsed arguments, not config-file comments.
			if(!line.isEmpty() && !line.startsWith("#")){replay.add(line);}
		}
		for(int i=0; i<replay.size(); i++){
			if(replay.get(i).startsWith("outdir=")){replay.set(i,"outdir="+root.resolve("fission_recipe_replay"));}
		}
		new DdlClustCake(replay.toArray(new String[replay.size()]),System.err);
		final Path explicitOff=root.resolve("fission_first_only_explicit_off");
		final ArrayList<String> off=fissionArguments(explicitOff,input,codec);
		off.add("fissiontopfractions=1,0,0");off.add("fissionfusionproductsonly=f");
		DdlClust.main(off.toArray(new String[off.size()]));
		checkArmParity(result,explicitOff);
		final Path restricted=root.resolve("fission_products_only");
		final ArrayList<String> productArgs=fissionArguments(restricted,input,codec);
		productArgs.add("fissionfusionproductsonly=t");
		DdlClust.main(productArgs.toArray(new String[productArgs.size()]));
		check(!hasSetting(restricted.resolve("layer1/run.config"),"fissionparents"),
			"The fusion-product-only Cake flag leaves the first layer unrestricted");
		for(int layer=1; layer<=3; layer++){
			final Path child=restricted.resolve("layer"+layer);
			check(!hasSetting(child.resolve("run.config"),"fusionproductsout") &&
				!Files.exists(child.resolve("fusion_products.tsv")),
				"Fusion-product eligibility stays in the resident state instead of an intermediate file");
			check(Files.readAllLines(child.resolve("state.tsv"),StandardCharsets.US_ASCII)
				.contains("fusion_product_families\t0"),
				"Restricted layer state records how many surviving fusion products were persisted");
			if(layer>1){
				check(!hasSetting(child.resolve("run.config"),"fissionparents"),
					"Post-first restricted layers use resident fusion tags without a disk input");
			}
		}
		for(int layer=1; layer<=3; layer++){
			final Path child=result.resolve("layer"+layer);
			final Path recipe=child.resolve("run.config");
			check(Double.parseDouble(setting(recipe,"fissiontopfraction"))==(layer==1 ? 1 : 0),
				"First-layer-only fission resolves independently at each actual child invocation");
			check(setting(recipe,"fusionbottomfraction").equals("1") && setting(recipe,"belowmode").equals("constant") &&
				setting(recipe,"sizelambda").equals("0.125") && setting(recipe,"setupslotcost").equals("3"),
				"Disabling fission does not alter fusion or layer pricing");
			final List<String> state=Files.readAllLines(child.resolve("state.tsv"),StandardCharsets.US_ASCII);
			check(state.contains("fission_waves\t0")== (layer>1),"Only the first layer evaluates fission requests");
			boolean fusion=false;
			for(String row:Files.readAllLines(child.resolve("waves.tsv"),StandardCharsets.US_ASCII)){
				if(row.startsWith("fusion\t")){fusion=true;}
			}
			check(fusion,"Every layer still executes the fusion phase");
		}
		for(String invalid:new String[]{"","null","none","1,0","1,0,0,0","1,0,","-0.01","1.01","NaN","Infinity","1,NaN,0"}){
			final Path bad=root.resolve("bad_fission_"+checks);
			final ArrayList<String> a=fissionArguments(bad,input,codec);a.add("fissiontopfractions="+invalid);
			checkRejected(a,bad,"Invalid fission schedule: "+invalid);
		}
		for(boolean scalarFirst:new boolean[]{false,true}){
			final Path bad=root.resolve("conflicting_fission_"+checks);
			final ArrayList<String> a=fissionArguments(bad,input,codec);
			a.add(scalarFirst ? "fissiontopfraction=1" : "fissiontopfractions=1,0,0");
			a.add(scalarFirst ? "fissiontopfractions=1,0,0" : "fissiontopfraction=1");
			checkRejected(a,bad,"Scalar and schedule conflict in either order");
		}
	}

	private static ArrayList<String> fissionArguments(final Path output,final Path input,final Path codec){
		final ArrayList<String> args=arguments(output,input,codec,1);
		args.remove("fissiontopfraction=0");args.remove("fusionbottomfraction=0");
		args.addAll(Arrays.asList("fractions=.25,.5,1","sizelambda=.125","belowmode=constant","minchildsize=1"));
		return args;
	}
	private static void checkRejected(final ArrayList<String> args,final Path output,final String reason){
		boolean rejected=false;
		try{DdlClust.main(args.toArray(new String[args.size()]));}catch(RuntimeException expected){rejected=true;}
		check(rejected && !Files.exists(output),reason+" must fail before output creation");
	}
	private static String setting(final Path recipe,final String key)throws Exception{
		String value=null;
		for(String line:Files.readAllLines(recipe,StandardCharsets.US_ASCII)){
			if(line.startsWith(key+"=")){check(value==null,"Child settings must occur exactly once: "+key);value=line.substring(key.length()+1);}
		}
		check(value!=null,"Recipe must record "+key);return value;
	}
	private static boolean hasSetting(final Path recipe,final String key)throws Exception{
		for(String line:Files.readAllLines(recipe,StandardCharsets.US_ASCII)){
			if(line.startsWith(key+"=")){return true;}
		}
		return false;
	}
	private static String tabSetting(final Path recipe,final String key)throws Exception{
		String value=null;
		for(String line:Files.readAllLines(recipe,StandardCharsets.US_ASCII)){
			if(line.startsWith(key+"\t")){check(value==null,"Child settings must occur exactly once: "+key);value=line.substring(key.length()+1);}
		}
		check(value!=null,"Recipe must record "+key);return value;
	}
	private static void checkArmParity(final Path first,final Path second){
		checkFileParity(first.resolve("layers.tsv"),second.resolve("layers.tsv"));
		for(String file:new String[]{"families.pddl.gz","members.tsv.gz"}){
			checkFileParity(first.resolve("final").resolve(file),second.resolve("final").resolve(file));
		}
		checkNoIntermediateModels(first,3);checkNoIntermediateModels(second,3);
	}
	private static void checkNoIntermediateModels(final Path root,final int layers){
		check(!Files.exists(root.resolve("selection/selected.pddl.gz")),"Selection must hand seeds directly to assignment");
		check(!Files.exists(root.resolve("final/assignment_seeds.pddl.gz")),"Final refinement must use resident seeds");
		for(int i=0;i<=layers;i++){
			final Path stage=root.resolve(i==0 ? "initial" : "layer"+i);
			for(String name:new String[]{"families.pddl.gz","members.tsv.gz","selected_seeds.tsv.gz","fission.tsv.gz","fusion.tsv.gz"}){
				check(!Files.exists(stage.resolve(name)),"Resident stage must not materialize "+stage.resolve(name));
			}
		}
	}
	private static void checkFileParity(final Path first,final Path second){
		final ByteFile a=ByteFile.makeByteFile(first.toString(),true), b=ByteFile.makeByteFile(second.toString(),true);
		try{
			while(true){
				final byte[] left=a.nextLine(), right=b.nextLine();
				check(Arrays.equals(left,right),"Constant-schedule parity: "+first+" and "+second);
				if(left==null){break;}
			}
		}finally{
			final boolean failedA=a.close(), failedB=b.close();
			check(!failedA && !failedB,"Parity input streams close successfully");
		}
	}

	private static ArrayList<String> arguments(final Path output,final Path input,final Path codec,final int threads){
		return new ArrayList<String>(Arrays.asList("mode=cake","memberranking=f","in="+input,"codecs="+codec,"outdir="+output,
			"familybudget=1","initialfraction=3","passes=3","t="+threads,
			"k=2","width=64","sketchseed=7","positional=f","partitions=0","overlap=0",
			"minlen=100","minspan=70","maxalignlen=100","minid=1","minblosumratio=0",
			"mergeminid=1","mergeminratio=0","trimleft=f","consensustrimleft=f",
			"passhits=0","greedyhits=0","maxhydrodelta=9","greedymaxhydrodelta=9",
			"fissiontopfraction=0","fusionbottomfraction=0"));
	}
	private static void checkQuality(final Path path)throws Exception{
		final List<String> rows=Files.readAllLines(path,StandardCharsets.US_ASCII);
		check(rows.size()==2,"The final quality table contains exactly one retained family");
		final String[] fields=rows.get(1).split("\t",-1);
		check(fields.length==6 && fields[0].equals("0") && fields[1].equals("3") && fields[2].equals("1") &&
			Double.parseDouble(fields[3])==1 && Double.parseDouble(fields[4])==1,
			"Three identical A tails have available identity1 and BLOSUM self-ratio1 at the actual final pass");
	}
	private static ProteinSequence sequence(final char residue){return new ProteinSequence("source",repeat(residue,100));}
	private static String repeat(final char residue,final int length){final char[] text=new char[length];Arrays.fill(text,residue);return new String(text);}
	private static void check(final boolean valid,final String reason){checks++;if(!valid){throw new AssertionError(reason);}}
	private static int checks;
	private DdlClustCakeTest(){}
}
