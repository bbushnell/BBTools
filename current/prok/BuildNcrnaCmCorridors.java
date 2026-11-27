package prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import aligner.CovarianceModel;
import aligner.CovarianceModelParser;
import aligner.CovarianceModelStateCorridors;
import map.IntHashMap;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaReferenceGeometry.Target;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;
import structures.IntList;

/** Fits state corridors from saved unrestricted paths, without realigning.
 * CLI and serial stream ownership follow template.A_SampleByteFile and the
 * shared native text owners. Model/call/group metadata stay on primitive axes;
 * two streamed path passes fit and evaluate without retaining path-row objects.
 * Five source-group folds precede fitting. Full-fit widths describe VISITED
 * state endpoint intervals, not all allocated DP cells or predicted speed.
 * @author Brian Bushnell, Raiden
 */
public final class BuildNcrnaCmCorridors{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, BuildNcrnaCmCorridors.class, false);
		try{new BuildNcrnaCmCorridors(pp.args).process();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private BuildNcrnaCmCorridors(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("calls")){callPath=b;}
			else if(a.equals("audit")){auditPath=b;}
			else if(a.equals("models")){modelPath=b;}
			else if(a.equals("cms")){cmPaths=b.split(",", -1);}
			else if(a.equals("refs")){
				require(b!=null && !b.isEmpty(), "Context resources need explicit original consensus FASTAs");
				referencePaths=b.split(",", -1);
				for(final String path:referencePaths){require(!path.isEmpty(), "Original consensus paths cannot contain empty entries");}
			}
			else if(a.equals("anchorpadding")){
				anchorPadding=Parse.parseIntKMG(b);
				anchorOption=true;
			}
			else if(a.equals("margin")){margin=Parse.parseIntKMG(b);}
			else if(a.equals("expectedcalls")){expectedCalls=Parse.parseIntKMG(b);}
			else if(a.equals("expectedsearches")){expectedSearches=Parse.parseIntKMG(b);}
			else if(a.equals("expectedsteps")){expectedSteps=Parse.parseKMG(b);}
			else if(a.equals("expectedmodels")){expectedModels=Parse.parseIntKMG(b);}
			else if(a.equals("expectedassemblies")){expectedAssemblies=Parse.parseIntKMG(b);}
			else if(a.equals("mincalls")){minCalls=Parse.parseIntKMG(b);}
			else if(a.equals("mingroups")){minGroups=Parse.parseIntKMG(b);}
			else if(a.equals("folds")){folds=Parse.parseIntKMG(b);}
			else if(a.equals("uniform")){uniform=Parse.parseIntKMG(b);}
			else if(a.equals("trimpercent")){trimPercent=Parse.parseIntKMG(b);}
			else{require(parser.parse(arg, a, b), "Unknown saved-path corridor option: "+arg);}
		}
		partsPath=parser.in1;
		out=parser.out1;
		require(partsPath!=null && callPath!=null && auditPath!=null && modelPath!=null && out!=null
			&& cmPaths!=null && cmPaths.length>0 && expectedCalls>0 && expectedSearches>0
			&& expectedSearches<=expectedCalls && expectedSteps>0 && expectedModels>0 && expectedAssemblies>0
			&& minCalls>0 && minGroups>0 && folds>=2 && folds<=7 && uniform>=0 && margin>=0 && trimPercent>=0 && trimPercent<100
			&& !parser.overwrite && !parser.append && parser.in2==null && parser.out2==null
			&& parser.maxReads<0 && parser.samplerate==1 && !new File(out).exists(),
			"Bind complete parts/call/model/audit axes, counts, CM files, fresh output and trimpercent=0..99; at most seven source folds fit one retained-mask byte");
		require(referencePaths==null ? !anchorOption : anchorOption && anchorPadding>=0,
			"Version2 context export requires refs=original_consensuses and anchorpadding=original_teacher_context together");
	}

	/*--------------------------------------------------------------*/
	/*----------------          Preparation         ----------------*/
	/*--------------------------------------------------------------*/

	private void process(){
		System.err.println("Whole-path trimming percent: "+trimPercent);
		readSources();
		readModels();
		if(referencePaths!=null){readReferences();}
		readOriginalCalls();
		readParts();
		readSummaries();
		require(new File(out).mkdir(), "Cannot create the fresh corridor-result directory");
		writeFolds();
		selectTrajectories();
		readPaths(true);
		freezeProfiles();
		readPaths(false);
		writeReports();
		System.err.println("NCRNA_CM_CORRIDOR_FIT_PASS_SAVED_PATHS_ONLY");
	}

	private void readSources(){
		blockedGroups=new boolean[expectedAssemblies];
		try(Input input=new Input(auditPath)){
			RrnaResourceIO.header(input.reader, AUDIT_HEADER, auditPath);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==8, "The resolved assembly audit must retain its eight declared fields");
				final int group=p.parseInt(4);
				require(group>=0 && group<expectedAssemblies && (p.termEquals("true", 5) || p.termEquals("false", 5))
					&& (p.termEquals("true", 6) || p.termEquals("false", 6)),
					"Source groups and exclusion flags must be actual resolved-audit values");
				probe.set(p, 0);
				final Target assembly=probe.copy();
				probe.set(p, 1);
				require(!sources.contains(probe), "One saved input ID cannot name two audit sources");
				blockedGroups[group]|=p.termEquals("true", 5) || p.termEquals("true", 6);
				sources.put(probe.copy(), new Source(assembly, group));
			}
		}
		require(sources.size()==expectedAssemblies, "The source-group audit must conserve the complete formal roster");
	}

	private void readModels(){
		cms=new CovarianceModel[cmPaths.length];
		uniformProfiles=new CovarianceModelStateCorridors[cms.length];
		for(int i=0; i<cms.length; i++){
			cms[i]=CovarianceModelParser.read(cmPaths[i]);
			require(cms[i].accession!=null && cms[i].accession.matches("[A-Za-z0-9][A-Za-z0-9_.-]*"),
				"Each refined CM needs a stable accession usable as its member-resource filename");
			for(int j=0; j<i; j++){require(!cms[i].accession.equals(cms[j].accession), "Refined CM accessions must be unique");}
			uniformProfiles[i]=new CovarianceModelStateCorridors(cms[i], new short[Math.multiplyExact(cms[i].states(), 4)], new boolean[cms[i].states()]);
		}
		try(Input input=new Input(modelPath)){
			RrnaResourceIO.header(input.reader, MODEL_HEADER, modelPath);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==6 && p.parseInt(0)==models.size(), "The census must retain its complete dense consensus order");
				final String member=p.parseString(1), name=p.parseString(2);
				int cm=-1;
				for(int i=0; i<cms.length; i++){if(member.equals(cms[i].accession)){cm=i;}}
				require(p.parseInt(3)>=0 && p.parseInt(4)>=0 && (p.termEquals("true", 5) || p.termEquals("false", 5)),
					"Census counts and pre-anchor support flags must remain explicit");
				probe.set(p, 2);
				require(!modelNames.contains(probe), "Consensus names must be unique before binding saved winner ordinals");
				final Target nameKey=probe.copy();
				final Model model=new Model(models.size(), name, member, nameKey, cm, p.parseInt(3), p.termEquals("true", 5), folds+1);
				modelNames.put(nameKey, model);
				models.add(model);
			}
		}
		require(models.size()==expectedModels, "A missing consensus must not silently change the resource model axis");
	}

	/** Bind sequences by exact original name; keep FASTA order for member exports. */
	private void readReferences(){
		require(referencePaths!=null && anchorPadding>=0, "Context fitting needs its declared original FASTAs and teacher padding");
		for(final String path:referencePaths){
			final NcrnaSupport.Library library=NcrnaSupport.loadLibrary(path);
			for(int m=0; m<library.names.length; m++){
				final byte[] name=library.names[m].getBytes(StandardCharsets.UTF_8);
				probe.set(name, 0, name.length);
				final Model model=modelNames.get(probe);
				require(model!=null && model.sequence==null,
					"Original consensus FASTAs must cover the frozen census once, without unknown names or duplicate assignments");
				model.sequence=library.bases[m];
				referenceOrder.add(model.ordinal);
			}
		}
		require(referenceOrder.size==models.size(), "Context export cannot omit even an unrefined or sparse census model's original sequence");
	}

	private void readOriginalCalls(){
		callModel=new int[expectedCalls];
		callGroup=new int[expectedCalls];
		starts=new int[expectedCalls];
		stops=new int[expectedCalls];
		minus=new boolean[expectedCalls];
		callSources=new Source[expectedCalls];
		seqids=new byte[expectedCalls][];
		int count=0, refined=0;
		try(Input input=new Input(callPath)){
			RrnaResourceIO.header(input.reader, PrepareNcrnaCmTraining.CALL_HEADER, callPath);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(count<expectedCalls && p.terms()==13 && p.parseInt(0)==count, "Original training calls must retain their full dense ordinal axis");
				final int m=p.parseInt(2);
				require(m>=0 && m<models.size(), "Every training winner must index the frozen consensus census");
				final Model model=models.get(m);
				probe.set(p, 3);
				final Source source=sources.get(probe);
				probe.set(p, 1);
				require(source!=null && !blockedGroups[source.group] && source.assembly.equals(probe) && p.termEquals(model.member, 8)
					&& wordEquals(model.nameKey, 9), "Training calls must bind their exact admitted source and winning consensus/member");
				starts[count]=p.parseInt(5);
				stops[count]=p.parseInt(6);
				require(starts[count]>0 && stops[count]>=starts[count] && (p.termEquals("+", 7) || p.termEquals("-", 7)),
					"Original training intervals require positive ordered genomic coordinates and an explicit strand");
				minus[count]=p.termEquals("-", 7);
				seqids[count]=p.parseByteArray(4);
				callModel[count]=m;
				callGroup[count]=source.group;
				callSources[count]=source;
				model.calls.add(count++);
				if(model.cm>=0){refined++;}
			}
		}
		require(count==expectedCalls && refined==expectedSearches, "Original all-member and refined-member call totals must match the completed training pass");
		for(final Model model:models){require(model.calls.size==model.censusCalls, "Per-model accepted-call counts must agree with the unchanged census");}
	}

	private void readParts(){
		try(Input input=new Input(partsPath)){
			RrnaResourceIO.header(input.reader, PART_HEADER, partsPath);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==2, "Each completed part requires one call-summary path and one aggregate state-path file");
				parts.add(new String[]{p.parseString(0), p.parseString(1)});
			}
		}
		require(!parts.isEmpty(), "At least one completed saved-path part is required");
	}

	private void readSummaries(){
		seen=new boolean[expectedCalls];
		finite=new boolean[expectedCalls];
		lengths=new int[expectedCalls];
		departures=new int[expectedCalls];
		Arrays.fill(departures, -1);
		int count=0;
		for(final String[] part:parts){
			try(Input input=new Input(part[0])){
				RrnaResourceIO.header(input.reader, TraceNcrnaCmTraining.CALL_HEADER, part[0]);
				for(byte[] line; (line=input.reader.nextLine())!=null;){
					p.set(line);
					require(p.terms()==23, "Saved path summaries must retain all23fields including anchor/cost outcomes");
					final int call=p.parseInt(0);
					require(call>=0 && call<expectedCalls && !seen[call] && models.get(callModel[call]).cm>=0
						&& p.parseInt(2)==callModel[call], "Each refined call must occur once, on its original model ordinal");
					probe.set(p, 3);
					require(sources.get(probe)==callSources[call], "A saved path summary cannot substitute another source input");
					probe.set(p, 1);
					require(callSources[call].assembly.equals(probe) && fieldEquals(seqids[call], 4)
						&& p.parseInt(5)==starts[call] && p.parseInt(6)==stops[call]
						&& p.termEquals(minus[call] ? "-" : "+", 7), "Summary source and original genomic span/strand must match the complete call axis");
					final int left=p.parseInt(8), right=p.parseInt(9), length=p.parseInt(10);
					require(left>0 && left<=starts[call] && right>=stops[call] && (long)right-left+1==length,
						"Saved padded windows must contain their original call and retain their actual physical length");
					if(referencePaths!=null){
						require(left==Math.max(1L, (long)starts[call]-anchorPadding) && right<=(long)stops[call]+anchorPadding,
							"Context metadata must match the saved teacher's actual padded geometry, allowing physical contig-edge clipping");
					}
					finite[call]=p.termEquals("OK", 11);
					require(finite[call] || p.termEquals("IMPOSSIBLE", 11) || p.termEquals("UNSUPPORTED_SEQUENCE", 11), "Unknown exact-search status cannot become a fitted example");
					require(p.termEquals("true", 15) || p.termEquals("false", 15), "Anchor availability must be an explicit outcome");
					final boolean anchored=p.termEquals("true", 15);
					final int departure=p.parseInt(16);
					require(anchored ? finite[call] && departure>=0 : departure==-1, "Only finite independently anchored paths have a nonnegative departure");
					if(finite[call]){
						require(Float.isFinite(p.parseFloat(12)) && p.parseInt(13)>0 && p.parseInt(14)>=p.parseInt(13)
							&& p.parseInt(14)<=length, "A finite teacher result needs actual ordered endpoints inside its saved window");
						finiteCalls++;
					}
					if(anchored){usableCalls++;}
					lengths[call]=length;
					departures[call]=departure;
					maxLength=Math.max(maxLength, length);
					seen[call]=true;
					count++;
				}
			}
		}
		require(count==expectedSearches, "Saved parts must conserve the complete searched-call population");
		for(int call=0; call<expectedCalls; call++){require(seen[call]==(models.get(callModel[call]).cm>=0), "Every refined call, and no unrefined call, needs exactly one search summary");}
		retained=new byte[expectedCalls];
		pathSteps=new int[expectedCalls];
		widths=new long[WIDTH_NAMES.length][2][Math.addExact(maxLength, 2)];
		// Exact locus framing is needed only for this metadata join, not the
		// large state streams. Release it before allocating fitted profiles.
		seqids=null;
		starts=stops=null;
		minus=null;
		callSources=null;
		sources.clear();
	}

	/*--------------------------------------------------------------*/
	/*----------------       Fit and evaluate       ----------------*/
	/*--------------------------------------------------------------*/

	private void writeFolds(){
		final boolean[] written=new boolean[expectedAssemblies];
		try(Output output=new Output(out+"/source_folds.tsv")){
			output.writer.println("resolved_group\tfold");
			for(final int group:callGroup){
				if(written[group]){continue;}
				written[group]=true;
				output.writer.print(row.clear().append(group).tab().append(group%folds).nl());
			}
		}
	}

	private void selectTrajectories(){
		final NcrnaCmCorridorSelection selector=new NcrnaCmCorridorSelection(minCalls, minGroups, folds, trimPercent);
		final IntHashMap before=new IntHashMap(), after=new IntHashMap();
		for(final Model model:models){
			if(model.cm<0){continue;}
			for(int fit=0; fit<=folds; fit++){
				final int excluded=fit-1, bit=1<<fit;
				before.clear();
				for(int i=0; i<model.calls.size; i++){
					final int call=model.calls.get(i), group=callGroup[call];
					if(departures[call]>=0 && (excluded<0 || group%folds!=excluded)){
						model.eligible[fit]++;
						before.increment(group);
					}
				}
				model.eligibleGroups[fit]=before.size();
				final int kept=selector.select(model.calls, callGroup, departures, excluded, retained, bit);
				if(kept==0){continue;}
				after.clear();
				int count=0;
				for(int i=0; i<model.calls.size; i++){
					final int call=model.calls.get(i);
					if((retained[call]&bit)!=0){
						after.increment(callGroup[call]);
						count++;
					}
				}
				require(count==kept && kept>=minCalls && after.size()>=minGroups, "Selected whole trajectories must retain both declared support minima");
				model.fits[fit]=new Bounds(cms[model.cm], kept, after.size());
			}
		}
	}

	private void readPaths(final boolean fitting){
		final boolean[] pathsSeen=new boolean[expectedCalls];
		long rows=0;
		for(final String[] part:parts){
			int current=-1, steps=0, maxChild=-1;
			boolean fullContains=true, foldContains=true;
			try(Input input=new Input(part[1])){
				RrnaResourceIO.header(input.reader, TraceNcrnaCmTraining.PATH_HEADER, part[1]);
				for(byte[] line; (line=input.reader.nextLine())!=null;){
					p.set(line);
					require(p.terms()==13, "Every saved state visit needs its original13path/anchor fields");
					final int call=p.parseInt(0);
					require(call>=0 && call<expectedCalls && finite[call] && seen[call] && p.parseInt(1)==callModel[call], "Only a finite teacher call may own a path on this model axis");
					if(call!=current){
						if(current>=0){finishPath(current, steps, maxChild, fitting, fullContains, foldContains);}
						require(!pathsSeen[call], "A teacher call's state visits must form one complete, unrepeated block");
						pathsSeen[call]=true;
						current=call;
						steps=0;
						maxChild=-1;
						fullContains=foldContains=true;
					}
					final Model model=models.get(callModel[call]);
					final CovarianceModel cm=cms[model.cm];
					final int state=p.parseInt(3), start=p.parseInt(4), end=p.parseInt(5);
					require(p.parseInt(2)==steps && state>=0 && state<cm.states() && start>=0 && start<=end
						&& end<=lengths[call] && (steps>0 || state==0) && (cm.type[state]!=CovarianceModel.E || start==end),
						"Teacher steps must retain root/order, actual consumed-prefix intervals and empty END semantics");
					final int left=p.parseInt(7), right=p.parseInt(8);
					require((left==-1 || left>steps) && (right==-1 || right>steps), "Saved child steps must point forward inside this finite trajectory");
					maxChild=Math.max(maxChild, Math.max(left, right));
					if(departures[call]>=0){
						for(int bound=0; bound<4; bound++){raw[bound]=p.parseInt(bound+9);}
						require(raw[0]>=0 && raw[0]<=raw[1] && raw[1]<=lengths[call] && raw[2]>=0 && raw[2]<=raw[3]
							&& raw[3]<=lengths[call], "Independent raw anchors must stay ordered inside the actual window");
						if(fitting){
							for(int fit=0; fit<=folds; fit++){
								if((retained[call]&(1<<fit))!=0){model.fits[fit].add(state, start, end, raw);}
							}
						}else{
							fullContains&=contains(model, 0, call, state, start, end, true);
							foldContains&=contains(model, 1+callGroup[call]%folds, call, state, start, end, false);
							project(uniformProfiles[model.cm], state, lengths[call], raw);
							width(6);
							if((retained[call]&1)!=0){require(fullContains, "Full-fit minima/maxima must contain every retained trajectory through the serving projection");}
						}
						// The summary departure is also an independent conservation
						// check of which whole trajectories the selector ranked.
						if(fitting){
							final int departure=Math.max(0, Math.max(Math.max(raw[0]-start, start-raw[1]), Math.max(raw[2]-end, end-raw[3])));
							observedDeparture=Math.max(observedDeparture, departure);
						}
					}else{
						for(int field=9; field<13; field++){p.setBounds(field);require(p.a()==p.b(), "Unavailable anchors must remain empty, never invented zero-coordinate bounds");}
					}
					steps++;
					rows++;
				}
			}
			if(current>=0){finishPath(current, steps, maxChild, fitting, fullContains, foldContains);}
		}
		require(rows==expectedSteps, "Both streamed passes must conserve the completed teacher's exact state-visit count");
		for(int call=0; call<expectedCalls; call++){require(pathsSeen[call]==finite[call], "Every finite search, and no impossible/unsupported search, needs one complete saved path");}
	}

	private void finishPath(final int call, final int steps, final int maxChild, final boolean fitting,
			final boolean fullContains, final boolean foldContains){
		require(steps>0 && maxChild<steps, "A completed finite trajectory must contain all its referenced child steps");
		if(fitting){
			pathSteps[call]=steps;
			require(departures[call]<0 || observedDeparture==departures[call], "Saved and reconstructed whole-path departures must agree before fitting");
			observedDeparture=0;
		}else{
			require(pathSteps[call]==steps, "The evaluation stream must retain the fitted teacher's per-call step counts");
			if(departures[call]<0){return;}
			final Model model=models.get(callModel[call]);
			final int full=model.fits[0]==null ? 2 : (retained[call]&1)!=0 ? 0 : 1;
			final int cv=model.fits[1+callGroup[call]%folds]==null ? 4 : 3;
			model.trajectories[full]++;
			model.contained[full]+=fullContains ? 1 : 0;
			model.trajectories[cv]++;
			model.contained[cv]+=foldContains ? 1 : 0;
		}
	}

	private boolean contains(final Model model, final int fit, final int call, final int state,
			final int start, final int end, final boolean full){
		final CovarianceModelStateCorridors profile=model.profiles[fit];
		final boolean learned=profile!=null && profile.observed(state);
		final CovarianceModelStateCorridors selected=profile==null ? uniformProfiles[model.cm] : profile;
		project(selected, state, lengths[call], raw);
		width((full ? 0 : 3)+(profile==null ? 2 : learned ? 0 : 1));
		return selected.contains(state, lengths[call], projected, start, end);
	}

	private void project(final CovarianceModelStateCorridors profile, final int state, final int length, final int[] rawBounds){
		assert(rawBounds.length==4) : "All fit/evaluation projections reuse the same four-bound consumed-prefix schema";
		System.arraycopy(rawBounds, 0, projected, 0, 4);
		profile.project(state, length, projected, uniform);
	}
	private void width(final int kind){
		assert(kind>=0 && kind<widths.length) : "Width observations must retain their learned/unseen/sparse and full/CV category";
		widths[kind][0][Math.max(0, projected[1]-projected[0]+1)]++;
		widths[kind][1][Math.max(0, projected[3]-projected[2]+1)]++;
	}

	private void freezeProfiles(){
		for(final Model model:models){
			for(int fit=0; fit<=folds; fit++){
				if(model.fits[fit]!=null){model.profiles[fit]=model.fits[fit].freeze(margin);}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	private void writeReports(){
		final String[] names=new String[models.size()], members=new String[models.size()];
		final CovarianceModel[] axis=new CovarianceModel[models.size()];
		final CovarianceModelStateCorridors[] profiles=new CovarianceModelStateCorridors[models.size()];
		final int[] calls=new int[models.size()], groups=new int[models.size()];
		long supportedCensusCalls=0;
		int supportedCensusModels=0, fittedModels=0;
		try(Output support=new Output(out+"/support.tsv"); Output containment=new Output(out+"/containment.tsv")){
			support.writer.println("model_ordinal0\tmodel\tmember\texcluded_fold\taccepted_calls\tusable_calls\tusable_groups\tretained_calls\tretained_groups\ttrimmed_calls\tstatus");
			containment.writer.println("model_ordinal0\tmodel\tmember\tstratum\tanchored_trajectories\twhole_paths_contained");
			for(int m=0; m<models.size(); m++){
				final Model model=models.get(m);
				names[m]=model.name;
				members[m]=model.member;
				axis[m]=model.cm<0 ? null : cms[model.cm];
				profiles[m]=model.profiles[0];
				if(model.cm>=0 && model.censusSupported){
					supportedCensusCalls+=model.censusCalls;
					supportedCensusModels++;
				}
				if(model.fits[0]!=null){
					calls[m]=model.fits[0].calls;
					groups[m]=model.fits[0].groups;
					fittedModels++;
				}
				for(int fit=0; fit<=folds; fit++){
					final Bounds bounds=model.fits[fit];
					support.writer.print(row.clear().append(m).tab().append(model.name).tab().append(model.member).tab().append(fit-1)
						.tab().append(model.censusCalls).tab().append(model.eligible[fit]).tab().append(model.eligibleGroups[fit])
						.tab().append(bounds==null ? 0 : bounds.calls).tab().append(bounds==null ? 0 : bounds.groups)
						.tab().append(bounds==null ? 0 : model.eligible[fit]-bounds.calls).tab()
						.append(model.cm<0 ? "UNREFINED_MEMBER" : bounds==null ? "UNIFORM_SPARSE" : "LEARNED").nl());
				}
				for(int kind=0; kind<STRATA.length; kind++){
					containment.writer.print(row.clear().append(m).tab().append(model.name).tab().append(model.member).tab().append(STRATA[kind])
						.tab().append(model.trajectories[kind]).tab().append(model.contained[kind]).nl());
				}
			}
		}
		final long resourceBytes=writeResources(names, members, axis, profiles, calls, groups);
		try(Output output=new Output(out+"/widths.tsv")){
			output.writer.println("population\tendpoint\tinterval_width_bases\tvisited_states\tqualification");
			for(int kind=0; kind<widths.length; kind++){
				for(int end=0; end<2; end++){
					for(int width=0; width<widths[kind][end].length; width++){
						if(widths[kind][end][width]==0){continue;}
						output.writer.print(row.clear().append(WIDTH_NAMES[kind]).tab().append(end==0 ? "start" : "end").tab().append(width)
							.tab().append(widths[kind][end][width]).append("\tVISITED_STATE_ENDPOINT_INTERVALS_NOT_TOTAL_DP_CELLS\n"));
					}
				}
			}
		}
		try(Output output=new Output(out+"/summary.tsv")){
			output.writer.println("all_training_calls\tsearched_calls\tfinite_paths\tusable_anchors\tunavailable_anchors\tstate_visits\tpre_anchor_supported_models\tpre_anchor_supported_calls\tfull_fit_models\tcorridor_gzip_bytes\tqualification");
			output.writer.print(row.clear().append(expectedCalls).tab().append(expectedSearches).tab().append(finiteCalls).tab().append(usableCalls)
				.tab().append(finiteCalls-usableCalls).tab().append(expectedSteps).tab().append(supportedCensusModels).tab().append(supportedCensusCalls)
				.tab().append(fittedModels).tab().append(resourceBytes)
				.append("\tTRAINING_AND_SOURCE_GROUP_CV_NOT_CALLER_ACCURACY\n"));
		}
	}

	/** The caller loads one member at a time. Give it exactly that consensus
	 * axis, while support/containment tables keep their full original ordinals.
	 * No duplicate combined resource or permissive partial reader is needed.
	 */
	private long writeResources(final String[] names, final String[] members, final CovarianceModel[] axis,
			final CovarianceModelStateCorridors[] profiles, final int[] calls, final int[] groups){
		assert(names.length==models.size()) : "Member exports must be subsets of the complete, already reported consensus axis";
		long total=0;
		final IntList selected=new IntList();
		try(Output output=new Output(out+"/resources.tsv")){
			output.writer.println("member\tmodels\tlearned_models\tfile\tstored_gzip_bytes");
			for(final CovarianceModel cm:cms){
				selected.clear();
				if(referencePaths==null){
					for(int m=0; m<axis.length; m++){if(axis[m]==cm){selected.add(m);}}
				}else{
					for(int r=0; r<referenceOrder.size; r++){
						final int m=referenceOrder.get(r);
						if(axis[m]==cm){selected.add(m);}
					}
				}
				require(selected.size>0, "Every supplied refined CM must have at least one consensus on the bound census axis");
				final String[] memberNames=new String[selected.size], memberIds=new String[selected.size];
				final CovarianceModel[] memberModels=new CovarianceModel[selected.size];
				final CovarianceModelStateCorridors[] memberProfiles=new CovarianceModelStateCorridors[selected.size];
				final int[] memberCalls=new int[selected.size], memberGroups=new int[selected.size];
				final byte[][] memberSequences=referencePaths==null ? null : new byte[selected.size][];
				int learned=0;
				for(int i=0; i<selected.size; i++){
					final int m=selected.get(i);
					memberNames[i]=names[m];
					memberIds[i]=members[m];
					memberModels[i]=cm;
					memberProfiles[i]=profiles[m];
					memberCalls[i]=calls[m];
					memberGroups[i]=groups[m];
					if(memberSequences!=null){memberSequences[i]=models.get(m).sequence;}
					if(profiles[m]!=null){learned++;}
				}
				final String file=cm.accession+".corridors.tsv.gz", path=out+'/'+file;
				if(memberSequences==null){
					NcrnaCmCorridorResource.write(path, memberNames, memberIds, memberModels, memberProfiles, memberCalls, memberGroups);
				}else{
					NcrnaCmCorridorResource.writeContext(path, memberNames, memberIds, memberModels, memberProfiles,
						memberCalls, memberGroups, memberSequences, anchorPadding);
				}
				final NcrnaCmCorridorResource.Data loaded=NcrnaCmCorridorResource.read(path, memberNames, memberIds, memberModels, memberSequences);
				require(memberSequences==null || loaded.unclippedAnchors && loaded.anchorPadding==anchorPadding,
					"The serving resource must retain the fitted anchor context and residual-before-clipping contract");
				for(int i=0; i<selected.size; i++){
					final CovarianceModelStateCorridors original=memberProfiles[i], copy=loaded.profiles[i];
					require(loaded.retainedCalls[i]==memberCalls[i] && loaded.resolvedGroups[i]==memberGroups[i]
						&& (copy==null)==(original==null), "The actual resource loader must preserve fitted support and fallback descriptors");
					if(original==null){continue;}
					for(int state=0; state<cm.states(); state++){
						require(copy.observed(state)==original.observed(state), "Resource loading must preserve unseen-state uniform fallback");
						if(original.observed(state)){
							for(int bound=0; bound<4; bound++){require(copy.offset(state, bound)==original.offset(state, bound), "Signed residual bounds must round-trip exactly through the serving reader");}
						}
					}
				}
				final long bytes=new File(path).length();
				require(bytes>0, "A completed per-member resource must have actual stored bytes");
				total=Math.addExact(total, bytes);
				output.writer.print(row.clear().append(cm.accession).tab().append(selected.size).tab().append(learned).tab().append(file).tab().append(bytes).nl());
			}
		}
		return total;
	}

	private boolean wordEquals(final Target expected, final int field){
		p.setBounds(field);
		final byte[] line=p.line();
		int stop=p.a();
		while(stop<p.b() && line[stop]>32){stop++;}
		probe.set(line, p.a(), stop);
		return expected.equals(probe);
	}
	private boolean fieldEquals(final byte[] expected, final int field){
		p.setBounds(field);
		assert(expected!=null) : "The original call's retained seqid bytes must survive until its exact summary join";
		if(p.b()-p.a()!=expected.length){return false;}
		for(int i=0; i<expected.length; i++){if(p.line()[p.a()+i]!=expected[i]){return false;}}
		return true;
	}
	private static void require(final boolean ok, final String why){if(!ok){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------        Inner classes         ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Source{
		Source(final Target assembly_, final int group_){
			assembly=assembly_;
			group=group_;
		}
		final Target assembly;
		final int group;
	}
	private static final class Model{
		Model(final int ordinal_, final String name_, final String member_, final Target nameKey_, final int cm_, final int calls_, final boolean supported_, final int fits_){
			assert(ordinal_>=0) : "Reference order joins existing dense census ordinals, never a FASTA position guessed as a model index";
			ordinal=ordinal_;
			name=name_;
			member=member_;
			nameKey=nameKey_;
			cm=cm_;
			censusCalls=calls_;
			censusSupported=supported_;
			fits=new Bounds[fits_];
			profiles=new CovarianceModelStateCorridors[fits_];
			eligible=new int[fits_];
			eligibleGroups=new int[fits_];
		}
		final String name, member;
		final int ordinal;
		byte[] sequence;
		final Target nameKey;
		final int cm, censusCalls;
		final boolean censusSupported;
		final IntList calls=new IntList();
		final Bounds[] fits;
		final CovarianceModelStateCorridors[] profiles;
		final int[] eligible, eligibleGroups;
		final long[] trajectories=new long[STRATA.length], contained=new long[STRATA.length];
	}
	private static final class Bounds{
		Bounds(final CovarianceModel model_, final int calls_, final int groups_){
			model=model_;
			calls=calls_;
			groups=groups_;
			offsets=new short[Math.multiplyExact(model.states(), 4)];
			observed=new boolean[model.states()];
		}
		void add(final int state, final int start, final int end, final int[] raw){
			assert(!frozen && state>=0 && state<observed.length && raw.length==4) : "Only training may change this CM's unexpanded four-anchor state fit before freezing";
			final int base=state*4;
			for(int bound=0; bound<4; bound++){
				final int value=(bound<2 ? start : end)-raw[bound];
				require(value>=Short.MIN_VALUE && value<=Short.MAX_VALUE, "Version1 learned residuals cannot silently truncate outside signed16-bit storage");
				if(!observed[state] || ((bound&1)==0 ? value<offsets[base+bound] : value>offsets[base+bound])){offsets[base+bound]=(short)value;}
			}
			observed[state]=true;
		}
		CovarianceModelStateCorridors freeze(final int margin){
			require(!frozen && margin>=0, "A fit freezes once, with a nonnegative selected margin applied exactly once");
			for(int state=0; margin>0 && state<observed.length; state++){
				if(!observed[state]){continue;}
				for(int bound=0; bound<4; bound++){
					final int index=state*4+bound;
					final long value=(long)offsets[index]+((bound&1)==0 ? -margin : margin);
					require(value>=Short.MIN_VALUE && value<=Short.MAX_VALUE, "The selected margin must fit signed16-bit residual storage without truncation");
					offsets[index]=(short)value;
				}
			}
			frozen=true;
			return new CovarianceModelStateCorridors(model, offsets, observed);
		}
		final CovarianceModel model;
		final int calls, groups;
		final short[] offsets;
		final boolean[] observed;
		boolean frozen;
	}

	private String partsPath, callPath, auditPath, modelPath, out;
	private String[] cmPaths, referencePaths;
	private int anchorPadding=-1, margin;
	private boolean anchorOption;
	private int expectedCalls, expectedSearches, expectedModels, expectedAssemblies, minCalls=100, minGroups=10, folds=5, uniform=32, trimPercent=1;
	private long expectedSteps, finiteCalls, usableCalls;
	private int maxLength, observedDeparture;
	private int[] callModel, callGroup, starts, stops, lengths, departures, pathSteps;
	private boolean[] minus, seen, finite, blockedGroups;
	private byte[] retained;
	private Source[] callSources;
	private byte[][] seqids;
	private CovarianceModel[] cms;
	private CovarianceModelStateCorridors[] uniformProfiles;
	private long[][][] widths;
	private final ArrayList<Model> models=new ArrayList<Model>();
	private final ArrayList<String[]> parts=new ArrayList<String[]>();
	private final ObjectMap<Target, Source> sources=new ObjectMap<Target, Source>(Target.class, Source.class);
	private final ObjectMap<Target, Model> modelNames=new ObjectMap<Target, Model>(Target.class, Model.class);
	private final IntList referenceOrder=new IntList();
	private final Target probe=new Target();
	private final LineParser1 p=new LineParser1('\t');
	private final ByteBuilder row=new ByteBuilder();
	private final int[] raw=new int[4], projected=new int[4];
	private static final String AUDIT_HEADER="assembly_id\tinput_id\tcanonical_fasta\thistorical_alias_group\tresolved_group\theld_out\tduplicate_original_group\tsampled";
	private static final String MODEL_HEADER="model_ordinal0\tmember\tmodel\taccepted_calls\tresolved_assembly_groups\tlearned_corridor_supported";
	private static final String PART_HEADER="calls\tpaths";
	private static final String[] STRATA={"FULL_LEARNED_RETAINED", "FULL_LEARNED_TRIMMED", "FULL_UNIFORM_SPARSE", "CV_LEARNED", "CV_UNIFORM_SPARSE"};
	private static final String[] WIDTH_NAMES={"FULL_LEARNED_STATE", "FULL_UNSEEN_STATE", "FULL_SPARSE_MODEL", "CV_LEARNED_STATE", "CV_UNSEEN_STATE", "CV_SPARSE_MODEL", "UNIFORM_CONTROL"};
}
