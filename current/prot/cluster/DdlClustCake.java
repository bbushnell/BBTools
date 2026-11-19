package prot.cluster;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;

import ddl.ProteinPositionMode;
import ddl.DDLKeyTieBreak;
import ddl.DDLFamilySearch.StructuralInputs;
import fileIO.ByteStreamWriter;
import parse.Parse;
import parse.Parser;
import prot.ProteinSubstitutionMatrix;
import shared.Shared;
import structures.ByteBuilder;

/** One command runs a complete saved-greedy structural arm using existing native
 * engines and one resident population. Layers pass frozen state directly; only
 * final bulk artifacts and small progress/configuration summaries are written.
 * @author Brian Bushnell, Keqing
 */
public final class DdlClustCake{

	/** Called after PreParser expansion, so config= and direct flags behave alike. */
	public static boolean requested(final String[] args){
		String mode=null;
		for(String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair[0].equalsIgnoreCase("mode")){mode=pair.length==2 ? pair[1] : null;}
		}
		return "cake".equalsIgnoreCase(mode);
	}

	public DdlClustCake(final String[] args, final PrintStream stream){
		outstream=stream;
		final Parser parser=new Parser();
		defaults();
		boolean scalarSet=false, scheduleSet=false;
		boolean fissionScalarSet=false;
		String aminoPartitionKey=null;
		ProteinPositionMode positionMode=null;
		Boolean positionalOption=null;
		for(String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String key=pair[0].toLowerCase(Locale.ROOT);
			final String value=pair.length==2 && !pair[1].equalsIgnoreCase("null") ? pair[1] : null;
			if(key.equals("mode")){require("cake".equalsIgnoreCase(value), "Expected mode=cake");}
			else if(key.equals("codecs")){codecs=value;}
			else if(key.equals("memberranking")){memberRanking=Parse.parseBoolean(value);}
			else if(key.equals("memberhydroweight")){memberHydroWeight=MemberRankInputs.checkedHydroWeight(Double.parseDouble(value));}
			else if(key.equals("seedmembers")){seedMembers=value;}
			else if(key.equals("seedmemberssha80")){seedMembersPin=value;}
			else if(key.equals("outdir")){outdir=value;}
			else if(key.equals("familybudget")){budget=Parse.parseIntKMG(value);}
			else if(key.equals("initialfraction")){initialFraction=Double.parseDouble(value);}
			else if(key.equals("setupslotcost")){
				require(value!=null,"setupslotcost requires auto or a nonnegative integer");
				if("auto".equalsIgnoreCase(value)){setupSlotCost=-1;}
				else{
					final long cost=Parse.parseKMG(value);
					require(cost>=0 && cost<=Integer.MAX_VALUE,"setupslotcost must lie in [0,Integer.MAX_VALUE]");
					setupSlotCost=(int)cost;
				}
			}
			else if(key.equals("fractions")){
				fractions="none".equalsIgnoreCase(value) ? new double[0] : nonnegativeList(value, "fractions");
			}else if(key.equals("lambda") || key.equals("sizelambda")){
				lambdas=nonnegativeList(value, key);
				require(lambdas.length==1, "lambda/sizelambda accepts one value; use sizelambdas for a schedule");
				scalarSet=true;
			}else if(key.equals("sizelambdas")){lambdas=nonnegativeList(value, key);scheduleSet=true;}
			else if(key.equals("fissiontopfractions")){fissionFractions=nonnegativeList(value,key);}
			else if(key.equals("fissiontopfraction")){layer.put(key,value);fissionScalarSet=true;}
			else if(key.equals("fissionfusionproductsonly")){fissionFusionProductsOnly=Parse.parseBoolean(value);}
			else if(key.equals("passes")){passes=Parse.parseIntKMG(value);}
			else if(key.equals("tailcrop")){tailCrop=Parse.parseIntKMG(value);}
			else if(key.equals("quality")){quality=value;}
			else if(key.equals("stopchangefraction")){stopChangeFraction=Double.parseDouble(value);}
			else if(key.equals("positionmode") || key.equals("positionalmode")){positionMode=ProteinPositionMode.parse(value);}
			else if(key.equals("positional")){positionalOption=Parse.parseBoolean(value);cluster.put(key,value);}
			else if(key.equals("aabits") || key.equals("aapartition")){
				require(aminoPartitionKey==null || aminoPartitionKey.equals(key), "Use either aabits or aapartition, not both");
				aminoPartitionKey=key;cluster.put(key,value);
			}
			else if(key.equals("verification") || key.equals("blacklistwords") || key.equals("blacklist_words") || key.equals("blacklistk") || key.equals("blacklistmode") || key.equals("blacklistprovenance")){common.put(key,value);}
			else if(key.equals("pinleft")){explicitPinLeft=true;common.put(key,value);}
			else if(key.equals("pinright")){explicitPinRight=true;common.put(key,value);}
			else if(key.equals("anchorterminus")){common.put(key,value);}
			else if(key.equals("alignmentmode")){TailAlignmentGate.parseRightMode(value);common.put(key,value);}
			else if(key.equals("minpaired")){common.put(key,value);}
			else if(key.equals("matrix") || key.equals("matrixfile")){common.put(key,value);}
			else if(key.equals("consensusbonus")){
				final long bonus=Parse.parseKMG(value);
				require(bonus>=0 && bonus<=ddl.DDLKeyTieBreak.MAX_VOTES,"Consensus member bonus must fit nonnegative40-bit scores");
				if(bonus>0){common.put(key,value);}else{common.remove(key);}
			}
			else if(key.equals("keytiebreak") || key.equals("rawmanifest") || key.equals("expectedrawtails") ||
				key.equals("rawminlen") || key.equals("rawmaxlen") || key.equals("rawtailcrop")){
				common.put(key,value);
			}
			// Retain explicit backend choices for child calls and standalone recipe replay.
			// Do not add absent options to defaults: historical recipe bytes stay unchanged.
			else if(key.equals("bgzip") || key.equals("nativebgzip")){compression.put(key,value);}
			else if(common.containsKey(key)){common.put(key, value);}
			else if(cluster.containsKey(key)){cluster.put(key, value);}
			else if(layer.containsKey(key)){layer.put(key, value);}
			else if(compression.containsKey(key)){compression.put(key, value);}
			else if(parser.parse(arg, key, value)){}
			else{throw new IllegalArgumentException("Unknown cake parameter: "+arg);}
		}
		if(positionMode!=null){
			positionMode=ProteinPositionMode.resolve(positionMode.flag,null,positionalOption,ProteinPositionMode.FIXED);
			cluster.put("positionmode",positionMode.flag);
			cluster.put("positional",Boolean.toString(positionMode.positional()));
		}
		require(!scalarSet || !scheduleSet, "Specify either lambda/sizelambda or sizelambdas, not both");
		require(memberRanking || memberHydroWeight==MemberRankInputs.DEFAULT_HYDRO_WEIGHT,
			"Nondefault memberhydroweight requires memberranking=t");
		require(!fissionScalarSet || fissionFractions==null,
			"Specify either fissiontopfraction or fissiontopfractions, not both");
		require(parser.in1!=null && parser.in2==null && parser.out1==null && parser.out2==null &&
			!parser.append && parser.maxReads<0, "Cake requires one complete in= population and outdir=; no paired, out=, append or reads limit");
		input=parser.in1;
		require(tailCrop>=0, "tailcrop must be nonnegative; zero preserves complete input proteins");
		require(outdir!=null && !outdir.isEmpty() && codecs!=null && !codecs.isEmpty(), "Cake requires outdir and codecs from saved greedy clustering");
		require(fractions!=null, "Cake requires explicit fractions=list, or fractions=none for the no-cake control");
		validateSchedule(fractions, lambdas);
		if(fissionFractions!=null){
			require(fissionFractions.length==1 || fissionFractions.length==fractions.length,
				"fissiontopfractions needs one value or exactly one value per fraction");
			for(double value:fissionFractions){
				require(Double.isFinite(value) && value>=0 && value<=1,"fissiontopfractions values must be finite and in [0,1]");
			}
		}
		ClusterSeedSelector.requestedCount(budget, initialFraction);
		require(passes>=1 && passes<Integer.MAX_VALUE, "Final refinement requires a positive bounded passes count");
		require(Double.isFinite(stopChangeFraction) && stopChangeFraction>=0 && stopChangeFraction<1,
			"stopchangefraction must be finite and in [0,1)");
		threads=Shared.threads();
		require(threads>0, "Cake needs at least one worker");
		if(common.get("frozenindexthreads")==null){common.put("frozenindexthreads", Integer.toString(threads));}
		final int indexThreads=Parse.parseIntKMG(common.get("frozenindexthreads"));
		require(indexThreads>0 && indexThreads<=threads, "frozenindexthreads must lie in [1,t]");
		validateLayerSettings();
		if(common.get("matrixfile")!=null){
			require(Files.isRegularFile(Paths.get(common.get("matrixfile"))),"Matrix resource must be a regular disk file");
		}
		matrix=ProteinSubstitutionMatrix.load(common.get("matrix"),common.get("matrixfile"));
		if(common.get("matrixfile")==null && matrix.name.equals(ProteinSubstitutionMatrix.DEFAULT_NAME)){
			// Preserve the resolved input for every child's output-alias protection and replay.
			common.put("matrixfile",ProteinSubstitutionMatrix.defaultResourcePath());
		}
		validateEngineSettings();
		normalizeAlignmentSettings();
		final ddl.ProteinKmerGeometry geometry=ddl.ProteinKmerGeometry.make(Parse.parseIntKMG(cluster.get("k")),
			cluster.containsKey("aabits") ? cluster.get("aabits") : cluster.get("aapartition"));
		final ddl.ProteinBlacklistPolicy.Descriptor policy=ddl.ProteinBlacklistPolicy.policy(geometry);
		if(policy!=null){
			if(common.get("blacklist")==null){common.put("blacklist",policy.resource().toString());}
			final String identity=common.get("blacklistsha80");
			require(identity==null || identity.equals(policy.sha80),"Cake requires the "+policy.label+" blacklist");
			for(String key:new String[]{"blacklistwords","blacklist_words"}){
				require(!common.containsKey(key) || Parse.parseIntKMG(common.get(key))==policy.words,
					"Cake requires "+policy.words+" words for "+policy.label);
			}
			common.put("blacklistsha80",policy.sha80);
			common.remove("blacklist_words");common.put("blacklistwords",Integer.toString(policy.words));
		}
		for(Map.Entry<String,String> e:compression.entrySet()){Parser.parseZip(e.getKey()+"="+e.getValue(),e.getKey(),e.getValue());}
		for(Map.Entry<String,String> e:common.entrySet()){
			require(e.getValue()!=null || e.getKey().equals("blacklist") || e.getKey().equals("blacklistsha80") || e.getKey().equals("matrixfile"),
				"Missing value for "+e.getKey());
		}
		sources=codecs.split(",", -1);
		if(memberRanking){
			final boolean hydrating=seedMembers!=null || seedMembersPin!=null;
			require(!hydrating || sources.length==1,
				"Seedmember hydration cannot unambiguously bind multiple ranked codecs");
			require(!hydrating || seedMembers!=null && Files.isRegularFile(Paths.get(seedMembers)) &&
				seedMembersPin!=null && seedMembersPin.matches("[0-9a-f]{20}"),
				"Ranked Cake requires pinned seedmembers from its original complete population");
		}else{require(seedMembers==null && seedMembersPin==null,"Seed hydration requires explicit memberranking=t");}
		for(String source:sources){require(!source.isEmpty() && Files.isRegularFile(Paths.get(source)), "Missing codec: "+source);}
		for(String file:ClusteringInput.inputFiles(input)){
			require(Files.isRegularFile(Paths.get(file)) && Files.isReadable(Paths.get(file)), "Cake input must be an existing readable FASTA: "+file);
		}
		root=Paths.get(outdir).toAbsolutePath().normalize();
		require(!Files.exists(root), "Cake outdir must be new; completed or partial runs are never overwritten: "+root);
	}

	/** Executes dependent stages in this JVM; no shell subprocess runs the biology. */
	public void process(){
		try{
			Files.createDirectories(root.getParent());
			Files.createDirectory(root);
			try(PhaseCostLog costs=new PhaseCostLog(path(root,"phase_costs.tsv"));
					PhaseCostLog.Span total=PhaseCostLog.begin("cake_process")){
				runStages();
			}
		}catch(IOException e){throw new IllegalStateException("Cake output failed: "+root,e);}
		catch(RuntimeException e){throw e;}
		catch(Exception e){throw new IllegalStateException("Cake resident stage failed: "+root,e);}
	}

	/** The outer cost interval excludes JVM startup and CLI parsing, and includes
	 * these nested stages. It must not be added to their inclusive measurements.
	 */
	private void runStages()throws Exception{
		writeConfiguration();
		final ClusterSeedSelector.Selection selection;
		if(memberRanking){selection=null;}// Hydration/selection follows DdlClust normalization and its sole population load.
		else{
			try(PhaseCostLog.Span phase=PhaseCostLog.begin("seed_selection")){
				selection=ClusterSeedSelector.select(codecs.split(",",-1),budget,initialFraction);
			}
		}
		final int cursor;
		final MemberRankInputs memberInputs;
		final Path initial=directory("initial");
		final ClusteringInput population;
		final ddl.FamilySketchBuilder builder;
		LiveProbeState state;
		try(PhaseCostLog.Span phase=PhaseCostLog.begin("initial_assignment")){
			final DdlClust.Assignment initialAssignment=assign(initial,selection==null ? null : selection.models,true,0,null,null,0);
			population=initialAssignment.input;builder=initialAssignment.builder;memberInputs=initialAssignment.memberInputs;
			keyContext=builder.policyContext();keyInputs=initialAssignment.keyInputs;
			require(keyContext!=null && keyInputs!=null,"Initial assignment must retain its key context and input tracker");
			cursor=memberRanking ? initialAssignment.nextFamilyId : selection.nextFamilyId;
			require(memberRanking==(memberInputs!=null) && (memberInputs==null || memberInputs.belongsTo(population.sequences)),
				"Initial assignment must return its ranking mode and original population cache");
			state=LiveProbeState.fromBaseline(LocalProbeBaseline.fromAssignment(initialAssignment.result,
				initialAssignment.header,population.sequences,Parse.parseIntKMG(common.get("minlen")),cursor),population.sequences.length);
			if(memberRanking){state=state.withMemberRanking(population.sequences,memberInputs,currentMemberLayer(state.models()));}
		}
		final Census setup=census(state.models());
		final int setupPrice=setupSlotCost<0 ? setup.budgetSize(budget) : setupSlotCost;
		final LiveProbeLayer.BelowMode belowMode=LiveProbeLayer.BelowMode.valueOf(layer.get("belowmode").toUpperCase(Locale.ROOT));
		final boolean usePrior=belowMode==LiveProbeLayer.BelowMode.PRIOR || belowMode==LiveProbeLayer.BelowMode.PRIORSCALED;
		final double edgeFraction=Double.parseDouble(layer.get("edgefraction"));
		int priorPrice=setupPrice, layerPrice=setupPrice;
		outstream.println("Cake setup_slot_cost="+setupPrice+" source="+(setupSlotCost<0 ? "initial_census" : "explicit"));
		final long records=population.sequences.length;
		int nextId=cursor;
		Census current=setup;
		final ByteStreamWriter counts=open(path(root,"layers.tsv"));
		try{
			counts.println("stage\tlayer\tfraction\tsize_lambda\tmodels\tactive_models\tassigned\tinput_records\tsetup_slot_cost\tnext_family_id");
			writeCount(counts,"initial",0,"NA","NA",current,records,layerPrice,nextId);
			for(int i=0; i<fractions.length; i++){
				final Path next=directory("layer"+(i+1));
				final double lambda=lambdaAt(lambdas,i,fractions.length);
				if(usePrior){
					priorPrice=priorSlotCost(priorPrice,current.positive,budget,current.budgetSize(budget),edgeFraction);
					layerPrice=priorPrice;
				}
				try(PhaseCostLog.Span phase=PhaseCostLog.begin("layer_"+(i+1))){
					keyInputs.verify();
					if(memberRanking){
						state=state.withMemberRanking(population.sequences,memberInputs,
							nextMemberLayerAfter(state.models(),state.memberRankingLayer()));
					}
					state=runLayer(next,state,population,builder,layerPrice,i,lambda).result().state();
					keyInputs.verify();
				}
				require(state.sequenceCount()==records,"Layer changed the input population");
				nextId=state.firstNewId();
				current=census(state.models());
				require(state.assigned()==current.assigned,"Layer state and resident model census differ");
				writeCount(counts,"layer",i+1,Double.toString(fractions[i]),Double.toString(lambda),
					current,records,layerPrice,nextId);
			}
			final Path result=directory("final");
			final ClusterSketchCodec.Loaded finalSeeds=ClusterSketchCodec.resident(state.sourceHeader(),state.models(),true);
			final DdlClust.Assignment finalAssignment;
			try(PhaseCostLog.Span phase=PhaseCostLog.begin("final_refinement")){
				finalAssignment=assign(result,finalSeeds,false,Math.min(budget,finalSeeds.families.length),population,memberInputs,
					memberRanking ? nextMemberLayerAfter(state.models(),state.memberRankingLayer()) : 0);
				require(finalAssignment.input==population && finalAssignment.memberInputs==memberInputs,
					"Final assignment must retain the identical original population and ranking cache");
				require(finalAssignment.builder.policyContext()==keyContext && finalAssignment.keyInputs==keyInputs,
					"Final assignment must reuse the initial key counts and their complete stability tracker");
				require(finalAssignment.builder.packedSelection()==builder.packedSelection(),
					"Final assignment must retain the initial explicit-study or implicit-default selection strategy");
			}
			final Census finished=census(finalAssignment.result.families);
			writeSizes(path(result,"family_sizes.tsv"),finished);
			require(finalAssignment.result.inputs()==records && finalAssignment.result.assigned==finished.assigned,
				"Final refinement must preserve the input population and publish its actual assignments");
			writeCount(counts,"final",fractions.length,"NA","NA",finished,records,layerPrice,nextId);
		}finally{close(counts);}
		outstream.println("Cake arm complete: "+root+"; final outputs are uncalibrated training results.");
	}

	private DdlClust.Assignment assign(final Path dest,final ClusterSketchCodec.Loaded seeds,
			final boolean initial,final int limit,final ClusteringInput input,final MemberRankInputs memberInputs,
			final int memberLayer){
		final ArrayList<String> args=base();
		append(args,cluster);
		if(tailCrop>0){args.add("tailcrop="+tailCrop);}// Only assignment loads; resident layers already borrow cropped proteins.
		args.addAll(Arrays.asList("seedlimit="+limit,"passes="+(initial ? 1 : passes),
			"mincluster=1","topn=0","preserveallseeds="+initial,"skipfinalrebuild="+initial,
			"stopchangefraction="+(initial ? 0 : stopChangeFraction),"refinement=null",
			"ddlout="+(initial ? "null" : path(dest,"families.pddl.gz")),"members="+(initial ? "null" : path(dest,"members.tsv.gz")),
			"selectedseeds=null","stats="+path(dest,"stats.tsv"),
			"passchanges="+(initial ? "null" : path(dest,"passchanges.tsv")),"checkpoints=null",
			"hist="+(initial ? "null" : path(dest,"hist.tsv")),
			"quality="+(initial || quality==null ? "null" : quality.equals("auto") ? path(dest,"quality") : quality),
			"consensus="+(initial ? "null" : path(dest,"consensus.faa.gz"))));
		if(initial){args.add("consensuspasses=2");}// Fixed initial assignment cannot rebuild, regardless of final refinement depth.
		if(memberRanking){
			if(initial){
				if(seedMembers!=null){args.add("seedmembers="+seedMembers);args.add("seedmemberssha80="+seedMembersPin);}
				return DdlClust.assignInitialRanked(args.toArray(new String[args.size()]),sources,budget,initialFraction,matrix);
			}
			return DdlClust.assignResident(args.toArray(new String[args.size()]),seeds,input,memberInputs,memberLayer,matrix,keyContext,keyInputs);
		}
		return DdlClust.assignResident(args.toArray(new String[args.size()]),seeds,input,null,0,matrix,keyContext,keyInputs);
	}

	private DdlClustStageB.LayerResult runLayer(final Path dest,final LiveProbeState state,
			final ClusteringInput input,final ddl.FamilySketchBuilder builder,
			final int setupPrice,final int index,final double lambda)throws Exception{
		require(memberRanking==(state.memberRankingInputs()!=null),
			"A structural layer must preserve the Cake ranking mode and original cache binding");
		require(builder.policyContext()==keyContext,"Every Cake layer must borrow the initial immutable key-policy context");
		final ArrayList<String> args=base();
		// Resolve into a fresh map so a layer cannot change the next layer's defaults.
		final Map<String,String> settings=new LinkedHashMap<String,String>(layer);
		if(fissionFractions!=null){
			settings.put("fissiontopfraction",Double.toString(fissionFractions[fissionFractions.length==1 ? 0 : index]));
		}
		append(args,settings);
		args.addAll(Arrays.asList("mode=liveprobe","phase=layer","familybudget="+budget,
			"slotfraction="+fractions[index],"sizelambda="+lambda,"setupslotcost="+setupPrice,
			"selectedcount="+state.models().length,"report=null","fusionreport=null",
			"passesreport="+path(dest,"waves.tsv"),"metrics="+path(dest,"metrics.tsv"),
			"sizehist="+path(dest,"sizes.tsv"),"stateout="+path(dest,"state.tsv")));
		final ByteStreamWriter recipe=open(path(dest,"run.config"));
		try{for(String arg:args){recipe.println(arg);}}finally{close(recipe);}
		return DdlClustStageB.runResidentLayer(args.toArray(new String[args.size()]),state,input,builder,
			fissionFusionProductsOnly && index>0 ? state.fusionProductIds() : null,matrix);
	}

	/** Continue after imported or newly measured member descriptors without wrapping. */
	private static int nextMemberLayer(final Family[] families,final int floor){
		require(floor>=0,"Member-ranking layer floor must be nonnegative");
		int layer=floor;
		for(Family family:families){
			if(!family.hasMemberRanking()){continue;}
			final int existing=family.memberRankingSnapshot().layer;
			require(existing<Integer.MAX_VALUE,"Member-ranking publication layer space is exhausted");
			layer=Math.max(layer,existing+1);
		}
		return layer;
	}

	private static int nextMemberLayerAfter(final Family[] families,final int current){
		require(current<Integer.MAX_VALUE,"Member-ranking publication layer space is exhausted");
		return nextMemberLayer(families,current+1);
	}

	private static int currentMemberLayer(final Family[] families){
		int layer=0;
		for(Family family:families){
			require(family.hasMemberRanking(),"Ranked Cake publication needs member descriptors");
			layer=Math.max(layer,family.memberRankingSnapshot().layer);
		}
		return layer;
	}

	/** Counts the actual frozen roster without serializing or reparsing model payloads. */
	private static Census census(final Family[] families){
		require(families!=null,"Resident census requires the completed model roster");
		final int[] sizes=new int[families.length];
		long assigned=0;int positive=0;
		for(int i=0; i<sizes.length; i++){
			final int n=families[i].size();sizes[i]=n;assigned=Math.addExact(assigned,n);
			if(n>0){positive++;}
		}
		Arrays.sort(sizes);
		return new Census(sizes,positive,assigned);
	}

	private static void writeCount(final ByteStreamWriter writer,final String stage,final int index,
			final String fraction,final String lambda,final Census census,final long records,final int price,final int cursor){
		require(census.assigned<=records,"A stage cannot assign more tails than its input population");
		writer.println(new ByteBuilder().append(stage).tab().append(index).tab().append(fraction).tab().append(lambda)
			.tab().append(census.models).tab().append(census.positive).tab().append(census.assigned)
			.tab().append(records).tab().append(price).tab().append(cursor));
	}

	/** Reuses the final resident census; no extra sequence read or alignment is needed. */
	private static void writeSizes(final String path,final Census census){
		final ByteStreamWriter writer=open(path);
		try{
			writer.println("size\tfamilies\tassigned_members");
			long total=0;
			for(int i=0; i<census.sizes.length;){
				final int size=census.sizes[i], start=i;
				while(i<census.sizes.length && census.sizes[i]==size){i++;}
				final long members=Math.multiplyExact((long)size,i-start);
				total=Math.addExact(total,members);
				writer.println(new ByteBuilder().append(size).tab().append(i-start).tab().append(members));
			}
			require(total==census.assigned,"Family histogram must conserve the final assigned-tail census");
		}finally{close(writer);}
	}

	private ArrayList<String> base(){
		final ArrayList<String> args=new ArrayList<String>();append(args,common);append(args,compression);
		args.add("in="+input);args.add("t="+threads);args.add("ow=f");
		args.add("memberranking="+memberRanking);
		args.add("memberhydroweight="+memberHydroWeight);
		// DdlClust limits streaming/compression threads internally; explicitly restore
		// the requested worker count at each native entry point.
		return args;
	}

	private static void append(final ArrayList<String> args,final Map<String,String> settings){
		for(Map.Entry<String,String> e:settings.entrySet()){args.add(e.getKey()+"="+(e.getValue()==null ? "null" : e.getValue()));}
	}

	private void writeConfiguration(){
		final ByteStreamWriter writer=open(path(root,"run.config"));
		try{
			writer.println("mode=cake\nin="+input+"\ncodecs="+codecs+"\noutdir="+root+"\nfamilybudget="+budget+
				"\ninitialfraction="+initialFraction+"\nfractions="+csv(fractions)+"\nsizelambdas="+csv(lambdas)+
				"\npasses="+passes+"\nstopchangefraction="+stopChangeFraction+"\nt="+threads);
			writer.println("quality="+quality);
			if(!matrix.same(ProteinSubstitutionMatrix.BLOSUM62)){
				writer.println("# alignment_matrix_sha80="+matrix.sha80+" scale="+matrix.scale+" gap="+matrix.linearGap);
			}
			if(tailCrop>0){writer.println("tailcrop="+tailCrop);}
			writer.println("memberranking="+memberRanking);
			writer.println("memberhydroweight="+memberHydroWeight);
			if(memberRanking){
				if(seedMembers!=null){writer.println("seedmembers="+seedMembers+"\nseedmemberssha80="+seedMembersPin);}
			}
			writer.println("setupslotcost="+(setupSlotCost<0 ? "auto" : Integer.toString(setupSlotCost)));
			writer.println("fissionfusionproductsonly="+fissionFusionProductsOnly);
			if(fissionFractions!=null){writer.println("fissiontopfractions="+csv(fissionFractions));}
			for(Map<String,String> values:Arrays.asList(common,cluster,layer,compression)){
				for(Map.Entry<String,String> e:values.entrySet()){
					if(fissionFractions!=null && e.getKey().equals("fissiontopfraction")){continue;}
					writer.println(e.getKey()+"="+e.getValue());
				}
			}
		}finally{close(writer);}
	}

	private static String csv(final double[] values){
		if(values.length==0){return "none";}
		final StringBuilder text=new StringBuilder();
		for(double value:values){if(text.length()>0){text.append(',');}text.append(value);}
		return text.toString();
	}

	/** One lambda repeats; every other schedule must bind exactly one value per layer. */
	public static void validateSchedule(final double[] fractions,final double[] lambdas){
		require(fractions!=null && lambdas!=null && lambdas.length>0,"Layer fractions and size strengths are required");
		require(lambdas.length==1 || lambdas.length==fractions.length,"sizelambdas needs one value or exactly one value per fraction");
		for(double value:fractions){require(Double.isFinite(value) && value>=0,"Layer fractions must be finite and nonnegative");}
		for(double value:lambdas){require(Double.isFinite(value) && value>=0,"Size strengths must be finite and nonnegative");}
	}

	public static double lambdaAt(final double[] values,final int layerIndex,final int layers){
		require(values!=null && (values.length==1 || values.length==layers) && layerIndex>=0 && layerIndex<layers,
			"A layer must resolve to exactly one configured size strength");
		return values[values.length==1 ? 0 : layerIndex];
	}

	/** Preserve the last above-edge B through equality and underfilled populations. */
	public static int priorSlotCost(final int priorSlotCost,final int positiveFamilies,final int familyBudget,
			final int budgetSlotCost,final double edgeFraction){
		require(priorSlotCost>=0 && budgetSlotCost>=0,"Prior and budget slot costs must be nonnegative");
		return LiveProbeLayer.aboveEdge(positiveFamilies,familyBudget,edgeFraction) ? budgetSlotCost : priorSlotCost;
	}

	private static double[] nonnegativeList(final String text,final String key){
		require(text!=null && !text.isEmpty(),"Missing "+key+" list");
		final String[] terms=text.split(",",-1);
		final double[] values=new double[terms.length];
		for(int i=0; i<terms.length; i++){
			require(!terms[i].isEmpty(),"Empty term in "+key);
			values[i]=Double.parseDouble(terms[i]);
			require(Double.isFinite(values[i]) && values[i]>=0,key+" must contain finite nonnegative values");
		}
		return values;
	}

	private void validateLayerSettings(){
		LiveProbeLayer.PriceMode.valueOf(layer.get("pricemode").toUpperCase(Locale.ROOT));
		LiveProbeLayer.BelowMode.valueOf(layer.get("belowmode").toUpperCase(Locale.ROOT));
		for(String key:new String[]{"fissiontopfraction","fusionbottomfraction","pricedfraction"}){
			final double value=Double.parseDouble(layer.get(key));
			require(Double.isFinite(value) && value>=0 && value<=1,key+" must be finite and in [0,1]");
		}
		final double edge=Double.parseDouble(layer.get("edgefraction"));
		require(Double.isFinite(edge) && edge>0,"edgefraction must be finite and positive");
	}

	/** Reject invalid engine settings before creating a partial arm or recording nonfinite values. */
	private void validateEngineSettings(){
		final int minLength=Parse.parseIntKMG(common.get("minlen")), minSpan=Parse.parseIntKMG(common.get("minspan"));
		require(minLength>=minSpan && minSpan>0,"Cake requires minlen>=minspan>0");
		final TailAlignmentGate.Mode mode=alignmentMode();
		require(mode!=TailAlignmentGate.Mode.FULL_QUERY_UNPINNED,"Cake layers require pinned alignment");
		final String denominator=common.get("denominator");
		require("aligned".equalsIgnoreCase(denominator) || "full".equalsIgnoreCase(denominator),"denominator must be aligned or full");
		final TailAlignmentGate.Denominator denom="aligned".equalsIgnoreCase(denominator) ?
			TailAlignmentGate.Denominator.ALIGNED_QUERY : TailAlignmentGate.Denominator.FULL_QUERY;
		final int maxLength=Parse.parseIntKMG(common.get("maxalignlen"));
		final int minPaired=common.containsKey("minpaired") ? Parse.parseIntKMG(common.get("minpaired")) : 0;
		new TailAlignmentGate(mode,denom,maxLength,minSpan,Double.parseDouble(common.get("minid")),Double.parseDouble(common.get("minblosumratio")),null,matrix,minPaired);
		new TailAlignmentGate(mode,denom,maxLength,minSpan,Double.parseDouble(layer.get("mergeminid")),Double.parseDouble(layer.get("mergeminratio")),null,matrix,minPaired);
		final boolean flex="flex".equalsIgnoreCase(common.get("passhits"));
		final double hydro=Double.parseDouble(common.get("maxhydrodelta"));
		new CandidateFilter(flex ? 4 : Parse.parseIntKMG(common.get("passhits")),0,hydro,hydro,flex)
			.validateWidth(Parse.parseIntKMG(cluster.get("width")));
		new LayerCakeClustering.Options(threads,Parse.parseIntKMG(common.get("frozenindexthreads")),
			Parse.parseIntKMG(layer.get("fissionminsize")),Parse.parseIntKMG(layer.get("minchildsize")),
			Parse.parseIntKMG(layer.get("memberneighbors")),Parse.parseIntKMG(layer.get("memberminhits")),
			Parse.parseIntKMG(layer.get("evidenceneighbors")),Parse.parseIntKMG(common.get("shortlist")),
			Double.parseDouble(layer.get("fusionsizeratio")),Parse.parseBoolean(common.get("consensustrimleft")),
			Float.parseFloat(common.get("trimdepth")),budget,Parse.parseIntKMG(layer.get("localrefinerounds")),
			Double.parseDouble(layer.get("pricedfraction")));
		require(Parse.parseIntKMG(cluster.get("consensuspasses"))>0,"consensuspasses must be positive");
		require(Parse.parseIntKMG(layer.get("fissionminsize"))==1,"Live cake layers require fissionminsize=1");
	}

	private void defaults(){
		put(common,"matrix="+ProteinSubstitutionMatrix.DEFAULT_NAME,"minlen=100","minspan=70","maxalignlen=2147483647","minid=0.35","minblosumratio=0.14",
			"pinright=t","pinleft=f","trimleft=t","denominator=aligned","consensustrimleft=t","trimdepth=0",
			"shortlist=50","passhits=flex","maxhydrodelta=0.55","blacklist=null","blacklistsha80=null","frozenindexthreads=null");
		put(cluster,"k=5","width=1024","sketchseed=0","positional=t","partitions=4","overlap=2",
			"seed=1701","consensuspasses=2","greedyhits=flex","greedymaxhydrodelta=0.53","indexed=t");
		put(layer,"pricemode=floor","belowmode=priorscaled","edgefraction=1.01","pricedfraction=0.99",
			"fissiontopfraction=1","fusionbottomfraction=1","fusionsizeratio=8","fissionminsize=1",
			"minchildsize=5","memberneighbors=8","memberminhits=1","evidenceneighbors=10",
			"mergeminid=0.5","mergeminratio=0.14","localrefinerounds=2");
		put(compression,"zl=9","forcebgzip=t","nativebgzipout=f");
	}

	private static void put(final Map<String,String> map,final String... pairs){
		for(String pair:pairs){final String[] terms=Parse.splitOnFirst(pair,'=');map.put(terms[0],"null".equals(terms[1]) ? null : terms[1]);}
	}
	private TailAlignmentGate.Mode alignmentMode(){
		final String anchor=common.get("anchorterminus");
		final AlignmentAnchor explicitAnchor=anchor==null ? null : parseAlignmentAnchor(anchor);
		final Boolean pinRight=explicitPinRight ? Parse.parseBoolean(common.get("pinright")) : null;
		final Boolean pinLeft=explicitPinLeft ? Parse.parseBoolean(common.get("pinleft")) : null;
		final AlignmentAnchor resolved=explicitAnchor==null ?
			resolveAlignmentAnchor(pinRight,pinLeft) :
			resolveAlignmentAnchor(explicitAnchor,pinRight,pinLeft);
		if(common.containsKey("alignmentmode")){
			require(resolved==AlignmentAnchor.RIGHT && Parse.parseBoolean(common.get("trimleft")) &&
				"aligned".equalsIgnoreCase(common.get("denominator")),
				"Explicit right-pinned alignmentmode requires pinright=t, trimleft=t and denominator=aligned");
			return TailAlignmentGate.parseRightMode(common.get("alignmentmode"));
		}
		return mode(resolved,Parse.parseBoolean(common.get("trimleft")));
	}
	private void normalizeAlignmentSettings(){
		final String anchor=common.get("anchorterminus");
		final AlignmentAnchor explicitAnchor=anchor==null ? null : parseAlignmentAnchor(anchor);
		final Boolean pinRight=explicitPinRight ? Parse.parseBoolean(common.get("pinright")) : null;
		final Boolean pinLeft=explicitPinLeft ? Parse.parseBoolean(common.get("pinleft")) : null;
		final AlignmentAnchor resolved=explicitAnchor==null ?
			resolveAlignmentAnchor(pinRight,pinLeft) :
			resolveAlignmentAnchor(explicitAnchor,pinRight,pinLeft);
		common.put("anchorterminus",resolved==AlignmentAnchor.LEFT ? "head" :
			resolved==AlignmentAnchor.RIGHT ? "tail" : "none");
		common.put("pinleft",Boolean.toString(resolved==AlignmentAnchor.LEFT));
		common.put("pinright",Boolean.toString(resolved==AlignmentAnchor.RIGHT));
	}
	private static TailAlignmentGate.Mode mode(final AlignmentAnchor anchor,final boolean trimLeft){
		if(anchor==AlignmentAnchor.LEFT){return TailAlignmentGate.Mode.PIN_LEFT_TRIM_RIGHT;}
		if(anchor==AlignmentAnchor.RIGHT){return trimLeft ? TailAlignmentGate.Mode.PIN_RIGHT_TRIM_LEFT :
			TailAlignmentGate.Mode.PIN_RIGHT_FULL_QUERY;}
		assert(anchor==AlignmentAnchor.NONE) : "Unhandled anchor terminus: "+anchor;
		return TailAlignmentGate.Mode.FULL_QUERY_UNPINNED;
	}
	private static AlignmentAnchor parseAlignmentAnchor(final String value){
		if(value.equalsIgnoreCase("left") || value.equalsIgnoreCase("start") ||
				value.equalsIgnoreCase("head")){return AlignmentAnchor.LEFT;}
		if(value.equalsIgnoreCase("right") || value.equalsIgnoreCase("end") ||
				value.equalsIgnoreCase("tail")){return AlignmentAnchor.RIGHT;}
		if(value.equalsIgnoreCase("none") || value.equalsIgnoreCase("unpinned")){
			return AlignmentAnchor.NONE;
		}
		throw new IllegalArgumentException("anchorterminus requires left/start/head, right/end/tail, or none/unpinned");
	}
	private static AlignmentAnchor resolveAlignmentAnchor(final Boolean pinRight,final Boolean pinLeft){
		if(Boolean.TRUE.equals(pinLeft) && Boolean.TRUE.equals(pinRight)){
			throw new IllegalArgumentException("pinleft=t and pinright=t are contradictory");
		}
		return Boolean.TRUE.equals(pinLeft) ? AlignmentAnchor.LEFT :
			Boolean.FALSE.equals(pinRight) ? AlignmentAnchor.NONE : AlignmentAnchor.RIGHT;
	}
	private static AlignmentAnchor resolveAlignmentAnchor(final AlignmentAnchor anchor,final Boolean pinRight,final Boolean pinLeft){
		if(Boolean.TRUE.equals(pinLeft) && anchor!=AlignmentAnchor.LEFT){
			throw new IllegalArgumentException("pinleft=t contradicts anchorterminus="+anchor);
		}
		if(Boolean.TRUE.equals(pinRight) && anchor!=AlignmentAnchor.RIGHT){
			throw new IllegalArgumentException("pinright=t contradicts anchorterminus="+anchor);
		}
		if(Boolean.FALSE.equals(pinLeft) && anchor==AlignmentAnchor.LEFT){
			throw new IllegalArgumentException("pinleft=f contradicts anchorterminus=LEFT");
		}
		if(Boolean.FALSE.equals(pinRight) && anchor==AlignmentAnchor.RIGHT){
			throw new IllegalArgumentException("pinright=f contradicts anchorterminus=RIGHT");
		}
		return anchor;
	}
	private Path directory(final String name)throws IOException{final Path dir=root.resolve(name);Files.createDirectory(dir);return dir;}
	private static String path(final Path dir,final String name){return dir.resolve(name).toString();}
	private static ByteStreamWriter open(final String name){final ByteStreamWriter w=new ByteStreamWriter(name,false,false,true);w.start();return w;}
	private static void close(final ByteStreamWriter writer){if(writer.poisonAndWait()){throw new IllegalStateException("Cake output writer failed");}}
	private static void require(final boolean condition,final String reason){if(!condition){throw new IllegalArgumentException(reason);}}

	private static final class Census{
		Census(final int[] sizes_,final int positive_,final long assigned_){sizes=sizes_;models=sizes.length;positive=positive_;assigned=assigned_;}
		int budgetSize(final int budget){return positive<budget ? 0 : sizes[sizes.length-budget];}
		final int[] sizes;
		final int models,positive;
		final long assigned;
	}
	private final Map<String,String> common=new LinkedHashMap<String,String>();
	private final Map<String,String> cluster=new LinkedHashMap<String,String>();
	private final Map<String,String> layer=new LinkedHashMap<String,String>();
	private final Map<String,String> compression=new LinkedHashMap<String,String>();
	private final PrintStream outstream;
	private final String input;
	private final Path root;
	private final int threads;
	private final ProteinSubstitutionMatrix matrix;
	private DDLKeyTieBreak.PolicyContext keyContext;
	private StructuralInputs keyInputs;
	private String codecs,outdir,quality="auto";
	private String[] sources;
	private boolean memberRanking=true;
	private double memberHydroWeight=MemberRankInputs.DEFAULT_HYDRO_WEIGHT;
	private String seedMembers,seedMembersPin;
	private int budget=80000,passes=3,setupSlotCost=-1;
	private int tailCrop=0;
	private double initialFraction=1.5,stopChangeFraction=0.001;
	private double[] fractions,lambdas={0.25};
	private double[] fissionFractions;
	private boolean fissionFusionProductsOnly, explicitPinLeft, explicitPinRight;
	private enum AlignmentAnchor{RIGHT, LEFT, NONE}
}
