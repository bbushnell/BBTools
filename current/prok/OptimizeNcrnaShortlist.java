package prok;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.nio.channels.FileChannel;
import java.nio.channels.FileLock;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import dna.Data;
import map.ObjectIntMap;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaJointPolicy.Member;
import prok.NcrnaPolicySearch.Measurement;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/** Measured local search over coupled ncRNA shortlist controls.
 * Based on template.A_SampleBasic2 for CLI/output ownership; one native sweep
 * child executes each round, not one JVM per genome. Root-relative reference
 * changes stay in the sweep's reports. Small completed rounds can be replayed
 * after interruption without repeating their scientific work.
 *
 * An external run recipe verifies the supplied source/build/resource/data pins
 * and holds the shared build lock. The identity manifest declares those pins;
 * this driver binds their text and its full experiment configuration on resume,
 * but does not authenticate the referenced runtime or bulk data itself.
 * @author Brian Bushnell, Raiden
 */
public final class OptimizeNcrnaShortlist{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final OptimizeNcrnaShortlist tool=new OptimizeNcrnaShortlist(args);
		try{tool.process();}
		finally{Shared.closeStream(tool.outstream);}
		timer.stop();
		System.err.println("Joint shortlist search elapsed: "+timer);
	}

	public OptimizeNcrnaShortlist(final String[] args){
		this(args, null);
	}

	/** Test-only evaluator injection; no CLI option enables it. Its contract and
	 * final results are labeled synthetic_test and never claim caller verification.
	 */
	public OptimizeNcrnaShortlist(final String[] args, final RoundRunner runner_){
		runner=runner_;
		backend=runner==null ? "native" : "synthetic_test";
		final PreParser pp=new PreParser(args, OptimizeNcrnaShortlist.class, false);
		outstream=pp.outstream;
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:pp.args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("caller")){caller=b;}
			else if(a.equals("reference")){reference=b;}
			else if(a.equals("familymap")){familymap=b;}
			else if(a.equals("initial")){initial=b;}
			else if(a.equals("grid")){grid=b;}
			else if(a.equals("identity")){identity=b;}
			else if(a.equals("target")){target=b;}
			else if(a.equals("grading")){grading=b;}
			else if(a.equals("targetrecall")){targetRecall=Double.parseDouble(b);}
			else if(a.equals("cost")){
				require("alignments".equals(b) || "modelbases".equals(b), "cost must be alignments or modelbases");
				modelBaseCost="modelbases".equals(b);
			}
			else if(a.equals("maxlost")){maxLost=Parse.parseKMG(b);}
			else if(a.equals("maxlastloss")){maxLastLoss=Parse.parseKMG(b);}
			else if(a.equals("maxmemberlastloss")){maxMemberLastLoss=Parse.parseKMG(b);}
			else if(a.equals("members")){members=b;}
			else if(a.equals("loci")){loci=b;}
			else if(a.equals("membersha80")){memberSha80=b;}
			else if(a.equals("locisha80")){lociSha80=b;}
			else if(a.equals("maxrounds")){maxRounds=Parse.parseIntKMG(b);}
			else if(a.equals("heapm")){heapMiB=Parse.parseKMG(b);}
			else if(a.equals("cachewindows")){cacheWindows=Parse.parseIntKMG(b);}
			else if(a.equals("cachelength")){cacheLength=Parse.parseIntKMG(b);}
			else if(a.equals("reusecandidates")){reuseCandidates=Parse.parseBoolean(b);}
			else if(a.equals("overlap")){overlap=Float.parseFloat(b);}
			else if(a.equals("resume")){resume=Parse.parseBoolean(b);}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		require(parser.out1!=null && !parser.overwrite && !parser.append, "Supply out=directory with ow=f and append=f");
		require(caller!=null && reference!=null && familymap!=null && initial!=null && grid!=null && identity!=null,
			"Supply caller, reference, familymap, initial, grid and identity files");
		require(target!=null && ("target".equals(grading) || "union".equals(grading)), "Supply target and explicit grading=target|union");
		require(Double.isFinite(targetRecall) && targetRecall>0 && targetRecall<=1 && maxLost>=-1,
			"targetrecall must be explicit in (0,1], relative to shipped TP; maxlost is -1 or a nonnegative root-loss limit");
		require(maxRounds>0 && heapMiB>0 && cacheWindows>=0 && cacheLength>0 && overlap>0 && overlap<=1,
			"Invalid round, heap, cache or overlap bound");
		require(maxLastLoss>=-1 && maxMemberLastLoss>=-1, "Coverage budgets are -1 (disabled) or nonnegative input counts");
		if(maxMemberLastLoss>=0){
			require(members!=null && loci!=null && memberSha80!=null && memberSha80.matches("[0-9a-f]{20}")
				&& lociSha80!=null && lociSha80.matches("[0-9a-f]{20}"),
				"Member protection requires members=, loci=, membersha80= and locisha80= verified by the external recipe");
			require(Tools.testInputFiles(false, true, members, loci), "Member support metadata must be readable");
		}else{
			require(members==null && loci==null && memberSha80==null && lociSha80==null,
				"Member support arguments require enabled maxmemberlastloss; do not silently ignore them");
		}
		require(Tools.testInputFiles(false, true, caller, reference, familymap, initial, grid, identity), "Search input files must be readable");
		root=new File(parser.out1).getAbsoluteFile();
		// HotSpot tokenizes JAVA_TOOL_OPTIONS, which carries the child CPU output.
		require(!root.getPath().matches(".*[\\s\"\\\\].*"), "Output directory cannot contain whitespace, quotes or backslashes");
		final NcrnaJointPolicy[] initialPolicies=NcrnaJointPolicy.load(initial);
		require(initialPolicies.length==2, "initial must contain exactly the shipped root and the adaptive starting center");
		baseline=renamed("shipped", initialPolicies[0]);
		search=new NcrnaPolicySearch(initialPolicies[1], NcrnaPolicySearch.loadGrid(grid, initialPolicies[1]));
		validateIdentity(identity);
		final String installed=Data.findPath("?dev/sweepncrnashortlist.sh", false);
		require(installed!=null, "Cannot locate native shortlist sweep launcher searched as ?dev/sweepncrnashortlist.sh");
		script=new File(installed).getAbsoluteFile();
		callerArguments=Parser.parseConfig(new String[]{"config="+NcrnaCallerProcess.configurationPath(new File(caller))});
		coverageMembers=maxMemberLastLoss<0 ? new String[0] : NcrnaMemberCoverage.names(members,
			new NcrnaShortlistTarget(target, grading, true).referenceFamily);
		coverage=maxLastLoss<0 && maxMemberLastLoss<0 ? null : new NcrnaSearchCoverage(callerArguments, reference, coverageMembers, modelBaseCost);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Search Loop         ----------------*/
	/*--------------------------------------------------------------*/

	/** Executes or resumes under one exclusive directory lock; main owns outstream. */
	public void process(){
		try{
			if(resume){require(root.isDirectory(), "Resume requires an existing search directory");}
			else{Files.createDirectory(root.toPath());}
			try(FileChannel channel=FileChannel.open(new File(root, "search.lock").toPath(), StandardOpenOption.CREATE, StandardOpenOption.WRITE);
				FileLock lock=channel.tryLock()){
				require(lock!=null, "Another process owns this search directory");
				bindContract();
				if(coverage!=null){coverage.bindBaseline(new File(root, "coverage_baseline.tsv").getPath(), resume);}
				validateGridRoundTrip();
				iterate();
			}
		}catch(IOException e){throw new IllegalStateException("Could not preserve search state at "+root, e);}
	}

	private void iterate() throws IOException{
		final ObjectMap<String, Measurement> measured=new ObjectMap<String, Measurement>(String.class, Measurement.class);
		final ObjectIntMap<String> centers=new ObjectIntMap<String>(String.class);
		int[] center=search.initial();
		Measurement rootMeasurement=null, last=null;
		long targetTp=-1;
		for(int round=0; round<maxRounds; round++){
			final String centerKey=search.key(center);
			require(centers.get(centerKey)<0, "Search attempted to revisit a center; stop rather than cycle");
			centers.put(centerKey, round);
			final ArrayList<int[]> neighbors=new ArrayList<int[]>();
			final boolean[] freshNeighbors=new boolean[search.dimensions()*2];
			final ArrayList<int[]> newPoints=new ArrayList<int[]>();
			final ArrayList<NcrnaJointPolicy> policies=new ArrayList<NcrnaJointPolicy>();
			policies.add(baseline);
			policies.add(search.policy("center", center));
			for(int axis=0; axis<search.dimensions(); axis++){
				for(int direction=-1; direction<=1; direction+=2){
					final int[] next=search.neighbor(center, axis, direction);
					if(next==null || centers.get(search.key(next))>=0){continue;}
					neighbors.add(next);
					final boolean fresh=measured.get(search.key(next))==null;
					freshNeighbors[neighbors.size()-1]=fresh;
					if(fresh){
						newPoints.add(next);
						policies.add(search.policy("neighbor"+newPoints.size(), next));
					}
				}
			}
			final File directory=new File(root, "round"+round);
			final Measurement[] result=evaluate(directory, policies);
			if(rootMeasurement==null){
				rootMeasurement=result[0];
				targetTp=BigDecimal.valueOf(rootMeasurement.tp).multiply(BigDecimal.valueOf(targetRecall))
					.setScale(0, RoundingMode.CEILING).longValueExact();
			}else{same(rootMeasurement, result[0], "Shipped objective changed between rounds");}
			last=result[1];
			final Measurement previous=measured.get(centerKey);
			if(previous!=null){same(previous, last, "Center objective changed between rounds");}
			measured.put(centerKey, last);
			for(int i=0; i<newPoints.size(); i++){measured.put(search.key(newPoints.get(i)), result[i+2]);}
			final Measurement[] values=new Measurement[neighbors.size()];
			for(int i=0; i<values.length; i++){values[i]=measured.get(search.key(neighbors.get(i)));}
			final int best=NcrnaPolicySearch.choose(last, values, targetTp, maxLost, maxLastLoss, maxMemberLastLoss);
			final boolean feasible=NcrnaPolicySearch.constraintsMet(last, maxLost, maxLastLoss, maxMemberLastLoss);
			final String stop=best<0 ? (last.tp>=targetTp && feasible ? "LOCAL_POINT_MET_TARGET"
				: coverage!=null && !feasible ? "LOCAL_STOP_CONSTRAINTS_UNMET" : "LOCAL_STOP_TARGET_UNREACHED") : "MOVE";
			final String nextKey=best<0 ? centerKey : search.key(neighbors.get(best));
			bindText(new File(directory, "choice.tsv"), "status\tcenter\tnext\tTP\tFP\t"+costColumn()+"\tlost_root\ttarget_TP"+coverageHeader()+"\n"
				+stop+"\t"+centerKey+"\t"+nextKey+"\t"+last.tp+"\t"+last.fp+"\t"+last.cost+"\t"+last.lost+"\t"+targetTp+coverageText(last, targetTp)+"\n", directory.exists() && resume);
			if(coverage!=null){writeRejections(directory, centerKey, last, neighbors, values, freshNeighbors, best, targetTp);}
			outstream.println("Round "+round+": "+last.tp+" TP, "+last.fp+" apparent FP, "+last.cost
				+(modelBaseCost ? " aligned model bases; " : " logical alignments; ")+stop);
			if(best<0){finish(stop, center, last, targetTp);return;}
			if(round+1==maxRounds){finish("ROUND_LIMIT", neighbors.get(best), values[best], targetTp);return;}
			center=neighbors.get(best);
		}
		throw new AssertionError("A positive finite round bound must terminate with a measured result");
	}

	/** A completed marker is written only after the native child exits and parses.
	 * An existing incomplete directory is never silently rerun or overwritten.
	 */
	private Measurement[] evaluate(final File directory, final ArrayList<NcrnaJointPolicy> policies) throws IOException{
		final boolean existing=directory.exists();
		if(existing){
			require(resume && new File(directory, "complete.txt").isFile(),
				"Incomplete round requires inspection before retry; no completed work was repeated: "+directory);
			require(readText(new File(directory, "complete.txt").getPath()).equals(completedMarker()),
				"Invalid round completion marker: "+directory);
		}else{Files.createDirectory(directory.toPath());}
		final File policyFile=new File(directory, "policies.tsv"), config=new File(directory, "sweep.config");
		bindText(policyFile, policyText(policies), existing);
		final File report=new File(directory, "sweep");
		final String configuration="caller="+new File(root, "caller.config").getPath()+"\nreference="+new File(reference).getAbsolutePath()
			+"\nfamilymap="+new File(familymap).getAbsolutePath()+"\npolicies="+policyFile.getPath()+"\nout="+report.getPath()
			+"\ntarget="+target+"\ngrading="+grading+"\njoint=t\ncachewindows="+cacheWindows+"\ncachelength="+cacheLength
			+"\nreusecandidates="+reuseCandidates+"\noverlap="+overlap+childCoverage()+"\now=f\nappend=f\n";
		bindText(config, configuration, existing);
		if(!existing){
			if(runner!=null){runner.run(report, policies);}
			else{
				final ArrayList<String> command=new ArrayList<String>();
				command.add("env");
				final String prior=System.getenv("JAVA_TOOL_OPTIONS");
				final StringBuilder options=new StringBuilder(prior==null ? "" : prior+" ");
				for(final String option:java.lang.management.ManagementFactory.getRuntimeMXBean().getInputArguments()){
					if(option.startsWith("-Xss")){options.append(option).append(' ');}
				}
				options.append("-Dbbtools.diagnostic.ncrnaCpuOut=").append(new File(directory, "cpu.tsv").getPath());
				command.add("JAVA_TOOL_OPTIONS="+options);
				command.add("bash");
				command.add(script.getPath());
				command.add("simd="+Shared.SIMD);
				command.add("config="+NcrnaCallerProcess.configurationPath(config));
				NcrnaProcess.addJvmOptions(command, heapMiB);
				NcrnaProcess.run(command, null, "Joint shortlist round "+directory.getName());
			}
		}
		final Measurement[] result=readMeasurements(report.getPath()+".summary.tsv", policies);
		if(coverage!=null){
			final String[] names=new String[policies.size()];
			for(int p=0; p<names.length; p++){names[p]=policies.get(p).name;}
			final NcrnaSearchCoverage.Counts[] counts=coverage.read(report.getPath(), names, result);
			for(int p=0; p<result.length; p++){
				final Measurement point=result[p];
				result[p]=new Measurement(point.tp, point.fp, point.cost, point.lost, counts[p]);
			}
		}
		if(!existing){bindText(new File(directory, "complete.txt"), completedMarker(), false);}
		return result;
	}

	private String completedMarker(){
		return backend+(coverage==null ? "\tEVALUATION_RETURNED_AND_SUMMARY_PARSED\n"
			: "\tEVALUATION_RETURNED_AND_COVERAGE_VALIDATED\n");
	}

	/** Disabled protection adds no child flags, preserving the legacy round files. */
	private String childCoverage(){
		if(coverage==null){return "";}
		return "\ncoverage=t"+(members==null ? "" : "\nmembers="+new File(root, "members.tsv").getPath()
			+"\nloci="+new File(loci).getAbsolutePath());
	}

	/** Declared support pins are authenticated by the external recipe, as are the
	 * existing runtime/data pins. The small literal map is also bound on resume. */
	private String coverageContract(){
		if(coverage==null){return "";}
		return "coverage_schema=1\ncoverage_baseline=immutable_shipped_per_input\nmaxlastloss="+maxLastLoss
			+"\nmaxmemberlastloss="+maxMemberLastLoss+"\n"+(members==null ? "" : "members="+new File(members).getAbsolutePath()
			+"\nmembers_sha80="+memberSha80+"\nloci="+new File(loci).getAbsolutePath()+"\nloci_sha80="+lociSha80+"\n");
	}

	private String coverageHeader(){
		return coverage==null ? "" : "\tlast_loss_inputs\tfirst_gain_inputs\tmember_last_loss_inputs\tmember_first_gain_inputs\tconstraints_met\ttarget_met";
	}
	private String coverageText(final Measurement point, final long targetTp){
		if(coverage==null){return "";}
		require(point.coverage!=null, "Enabled coverage reports require validated input/member accounting");
		final StringBuilder text=new StringBuilder();
		text.append('\t').append(point.coverage.lastLoss()).append('\t').append(point.coverage.firstGain()).append('\t');
		for(int m=0; m<coverageMembers.length; m++){
			if(m>0){text.append(',');}
			text.append(point.coverage.memberLastLoss(m));
		}
		if(coverageMembers.length==0){text.append("NA");}
		text.append('\t');
		for(int m=0; m<coverageMembers.length; m++){
			if(m>0){text.append(',');}
			text.append(point.coverage.memberFirstGain(m));
		}
		if(coverageMembers.length==0){text.append("NA");}
		return text.append('\t').append(NcrnaPolicySearch.constraintsMet(point, maxLost, maxLastLoss, maxMemberLastLoss))
			.append('\t').append(point.tp>=targetTp).toString();
	}

	/** Reasons are a small diagnostic table, not an alternate selection objective. */
	private void writeRejections(final File directory, final String centerKey, final Measurement center,
			final ArrayList<int[]> neighbors, final Measurement[] points, final boolean[] fresh,
			final int selected, final long targetTp){
		assert(neighbors.size()==points.length && fresh.length>=points.length)
			: "Every considered neighbor has one cached or newly validated measurement";
		final StringBuilder text=new StringBuilder("point\torigin\tselected_move\treason\tTP\tFP\t"+costColumn()+"\tlost_root"+coverageHeader()+"\n");
		appendRejection(text, centerKey, "center", false, center, center, targetTp);
		for(int p=0; p<points.length; p++){
			appendRejection(text, search.key(neighbors.get(p)), fresh[p] ? "new" : "cached", p==selected, points[p], center, targetTp);
		}
		bindText(new File(directory, "coverage.tsv"), text.toString(), resume);
	}

	/** Recompute admissibility against today's center/TP floor without reevaluating
	 * a cached point. constraints_pass does not mean improving, selected or adopted. */
	private void appendRejection(final StringBuilder text, final String key, final String origin,
			final boolean selected, final Measurement point, final Measurement center, final long targetTp){
		require(point!=null && point.coverage!=null, "Every considered point needs its validated coverage vector");
		final StringBuilder reason=new StringBuilder();
		if(maxLost>=0 && point.lost>maxLost){reason.append("maxlost;");}
		if(maxLastLoss>=0 && point.coverage.lastLoss()>maxLastLoss){reason.append("maxlastloss;");}
		if(maxMemberLastLoss>=0){
			for(int m=0; m<coverageMembers.length; m++){
				if(point.coverage.memberLastLoss(m)>maxMemberLastLoss){reason.append("maxmemberlastloss[").append(m).append("];");}
			}
		}
		if(center.tp>=targetTp && point.tp<targetTp){reason.append("target_TP;");}
		if(reason.length()==0){reason.append("constraints_pass");}
		text.append(key).append('\t').append(origin).append('\t').append(selected).append('\t').append(reason)
			.append('\t').append(point.tp).append('\t').append(point.fp).append('\t').append(point.cost).append('\t').append(point.lost)
			.append(coverageText(point, targetTp)).append('\n');
	}

	/** Exercise the actual writer/reader for every declared control value before
	 * launching a native round; this validates numeric semantics, not file identity.
	 */
	private void validateGridRoundTrip(){
		final ArrayList<NcrnaJointPolicy> expected=search.serializationCases();
		final File path=new File(root, "grid_policies.tsv");
		bindText(path, policyText(expected), resume);
		final NcrnaJointPolicy[] actual=NcrnaJointPolicy.load(path.getPath());
		require(actual.length==expected.size(), "Serialized grid must retain every declared control choice");
		for(int p=0; p<actual.length; p++){
			final NcrnaJointPolicy wanted=expected.get(p), found=actual[p];
			require(wanted.name.equals(found.name) && wanted.size()==found.size(), "Grid serialization changed policy membership");
			for(int m=0; m<wanted.size(); m++){
				final Member a=wanted.member(m), b=found.member(m);
				require(a.name.equals(b.name) && a.adaptive==b.adaptive && a.topN==b.topN && a.minHits==b.minHits
					&& a.floor==b.floor && a.topFrac==b.topFrac && a.qFrac==b.qFrac && a.patience==b.patience,
					"Grid serialization changed an actual caller value: "+wanted.name+" "+a.name);
			}
		}
	}

	private void finish(final String status, final int[] center, final Measurement measured, final long targetTp){
		final ArrayList<NcrnaJointPolicy> result=new ArrayList<NcrnaJointPolicy>();
		result.add(baseline);
		result.add(search.policy("selected", center));
		bindText(new File(root, "selected.tsv"), policyText(result), resume);
		bindText(new File(root, "result.tsv"), "backend\tstatus\tTP\tFP\t"+costColumn()+"\tlost_root\ttarget_TP"+coverageHeader()+"\n"
			+backend+"\t"+status+"\t"+measured.tp+"\t"+measured.fp+"\t"+measured.cost+"\t"+measured.lost+"\t"+targetTp+coverageText(measured, targetTp)+"\n", resume);
		outstream.println("Result is a local measured point. Paired losses and apparent FP counts remain explicit; no production cutoff was adopted.");
	}

	/*--------------------------------------------------------------*/
	/*----------------          Persistence         ----------------*/
	/*--------------------------------------------------------------*/

	/** Binds explicit pins and all small configuration inputs without rereading bulk data. */
	private void bindContract(){
		if(resume){
			for(final String name:new String[]{"contract.txt", "caller.config", "identity.tsv", "initial.tsv", "grid.tsv", "familymap.tsv", "reference.tsv"}){
				require(new File(root, name).isFile(), "Resume requires the original complete contract, not reconstructed provenance: "+name);
			}
		}
		final String contract="schema=2\npatience=latched_eligible_raw_identity\nbackend="+backend+"\ntarget="+target+"\ngrading="+grading+"\ntargetrecall="+targetRecall
			+"\nmaxlost="+maxLost+"\nmaxrounds="+maxRounds+"\nheapm="+heapMiB+"\ncachewindows="+cacheWindows
			+"\ncachelength="+cacheLength+"\nreusecandidates="+reuseCandidates+"\noverlap="+overlap
			+"\nsimd="+Shared.SIMD+"\nlauncher="+script.getPath()+"\nreference="+new File(reference).getAbsolutePath()
			+"\nfamilymap="+new File(familymap).getAbsolutePath()+"\n"+coverageContract()
			+(modelBaseCost ? "cost=modelbases\n" : "");
		bindText(new File(root, "contract.txt"), contract, resume);
		bindText(new File(root, "caller.config"), String.join("\n", callerArguments)+"\n", resume);
		bindText(new File(root, "identity.tsv"), readText(identity), resume);
		bindText(new File(root, "initial.tsv"), readText(initial), resume);
		bindText(new File(root, "grid.tsv"), readText(grid), resume);
		bindText(new File(root, "familymap.tsv"), readText(familymap), resume);
		bindText(new File(root, "reference.tsv"), readText(reference), resume);
		if(members!=null){bindText(new File(root, "members.tsv"), readText(members), resume);}
	}

	/** Compare existing text only on resume; new files close once with native I/O. */
	private static void bindText(final File file, final String text, final boolean mayExist){
		if(file.exists()){
			require(mayExist && text.equals(readText(file.getPath())), "Resume configuration/result differs: "+file);
		}else{
			try(Output output=new Output(file.getPath())){output.writer.print(new ByteBuilder(text.getBytes(StandardCharsets.UTF_8)));}
		}
	}

	/** Small configuration/summary text only; this never reads genomes or GFFs. */
	private static String readText(final String path){
		final ByteBuilder text=new ByteBuilder();
		try(Input input=new Input(path)){
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){text.append(line).nl();}
		}
		return new String(text.toBytes(), StandardCharsets.UTF_8);
	}

	private static void validateIdentity(final String path){
		final String[] names={"source_commit", "build_seal", "resources_sha80", "inputs_sha80", "references_sha80"};
		final boolean[] seen=new boolean[names.length];
		final LineParser1 parser=new LineParser1('\t');
		try(Input input=new Input(path)){
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				parser.set(line);
				require(parser.terms()==2, "Identity manifest requires key and value, without a header");
				int key=-1;
				for(int i=0; i<names.length; i++){if(parser.termEquals(names[i], 0)){key=i;}}
				require(key>=0 && !seen[key], "Unknown or duplicate identity key");
				final String value=parser.parseString(1);
				require(value.matches(key==0 ? "[0-9a-f]{40}" : "[0-9a-f]{20}"), "Invalid source commit or sha80 identity");
				seen[key]=true;
			}
		}
		for(int i=0; i<seen.length; i++){require(seen[i], "Missing resume identity: "+names[i]);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Native Results        ----------------*/
	/*--------------------------------------------------------------*/

	/** The default metric and contract stay byte-compatible with prior searches. */
	private String costColumn(){return modelBaseCost ? "aligned_model_bases" : "logical_alignments";}

	private Measurement[] readMeasurements(final String path, final ArrayList<NcrnaJointPolicy> policies){
		final ObjectIntMap<String> indices=new ObjectIntMap<String>(String.class);
		for(int i=0; i<policies.size(); i++){indices.put(policies.get(i).name, i);}
		final Measurement[] result=new Measurement[policies.size()];
		final String[] fields={"policy", "TP", "FP", "FN", "scored_references", "lost_baseline_TP", "gained_TP",
			modelBaseCost ? "aligned_model_bases" : "kernel_logical_alignments"};
		final int[] columns=new int[fields.length];
		Arrays.fill(columns, -1);
		final LineParser1 parser=new LineParser1('\t');
		long scored=-1;
		try(Input input=new Input(path)){
			boolean header=false;
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				if(line.length>0 && line[0]=='#'){continue;}
				parser.set(line);
				if(!header){
					for(int c=0; c<parser.terms(); c++){
						for(int f=0; f<fields.length; f++){if(parser.termEquals(fields[f], c)){columns[f]=c;}}
					}
					for(int f=0; f<columns.length; f++){require(columns[f]>=0, "Missing native summary field: "+fields[f]);}
					header=true;
					continue;
				}
				final int index=indices.get(parser.parseString(columns[0]));
				require(index>=0 && result[index]==null, "Unexpected or duplicate native policy row");
				final long tp=parser.parseLong(columns[1]), fp=parser.parseLong(columns[2]), fn=parser.parseLong(columns[3]);
				final long denominator=parser.parseLong(columns[4]), lost=parser.parseLong(columns[5]), gained=parser.parseLong(columns[6]);
				require(fn>=0 && denominator>=0 && Math.addExact(tp, fn)==denominator, "Native reference denominator does not conserve TP+FN");
				if(scored<0){scored=denominator;}
				require(scored==denominator, "Native policies must share the same reference denominator");
				require(!parser.termEquals("NA", columns[7]), "Search requires measured "+fields[7]+", not legacy scavenger counts");
				result[index]=new Measurement(tp, fp, parser.parseLong(columns[7]), lost);
				if(index==0){require(lost==0 && gained==0, "Shipped policy must be its own immutable reference baseline");}
				else{
					require(result[0]!=null && tp==result[0].tp-lost+gained, "Paired reference changes must explain each policy's recovery");
				}
			}
		}
		for(final Measurement point:result){require(point!=null, "Every requested policy must have a completed native result");}
		return result;
	}

	private static void same(final Measurement before, final Measurement after, final String message){
		require(before.tp==after.tp && before.fp==after.fp && before.cost==after.cost && before.lost==after.lost,
			message+"; stop for investigation rather than rerunning to obtain matching numbers");
		require(NcrnaSearchCoverage.Counts.same(before.coverage, after.coverage), message+"; coverage counts changed");
	}

	private static NcrnaJointPolicy renamed(final String name, final NcrnaJointPolicy point){
		final Member[] members=new Member[point.size()];
		for(int i=0; i<members.length; i++){members[i]=point.member(i);}
		return new NcrnaJointPolicy(name, members);
	}

	/** Plain decimal preserves float round trips without requiring exponent support
	 * in LineParser1's native numeric parser. Small values must remain reloadable.
	 */
	private static String policyText(final ArrayList<NcrnaJointPolicy> policies){
		final StringBuilder text=new StringBuilder(NcrnaJointPolicy.HEADER).append('\n');
		for(final NcrnaJointPolicy policy:policies){
			for(int i=0; i<policy.size(); i++){
				final Member member=policy.member(i);
				text.append(policy.name).append('\t').append(member.name).append('\t').append(member.adaptive ? "adaptive" : "fixed")
					.append('\t').append(member.topN).append('\t').append(plain(member.floor)).append('\t').append(member.minHits)
					.append('\t').append(plain(member.topFrac)).append('\t').append(plain(member.qFrac))
					.append('\t').append(member.patience).append('\n');
			}
		}
		return text.toString();
	}
	private static String plain(final float value){return new BigDecimal(Float.toString(value)).toPlainString();}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/** Finite test seam. Implementations write the native summary schema only;
	 * the standard CLI always invokes the real whole-input caller instead.
	 */
	public interface RoundRunner{
		void run(File reportPrefix, ArrayList<NcrnaJointPolicy> policies);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String caller, reference, familymap, initial, grid, identity, target, grading;
	private double targetRecall=Double.NaN;
	private long maxLost=-1, heapMiB=4096;
	private long maxLastLoss=-1, maxMemberLastLoss=-1;
	private String members, loci, memberSha80, lociSha80;
	private final String[] coverageMembers;
	private final NcrnaSearchCoverage coverage;
	private int maxRounds=100, cacheWindows=128, cacheLength=4096;
	private float overlap=.5f;
	private boolean resume, reuseCandidates=true, modelBaseCost;
	private final File root, script;
	private final PrintStream outstream;
	private final NcrnaJointPolicy baseline;
	private final NcrnaPolicySearch search;
	private final String[] callerArguments;
	private final RoundRunner runner;
	private final String backend;
}
