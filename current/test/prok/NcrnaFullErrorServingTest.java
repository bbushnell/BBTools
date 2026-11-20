package test.prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import consensus.BaseGraph;
import fileIO.ReadWrite;
import ml.CellNet;
import ml.CellNetParser;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import parse.LineParser1;
import prok.FamilyResourceResolver;
import prok.GeneCaller;
import prok.NcrnaBoundaryTrace;
import prok.NcrnaFamily;
import prok.NcrnaFamilyLoader;
import prok.NcrnaFamilyRegistry;
import prok.NcrnaEndpointGateFeatures;
import prok.NcrnaEndpointModelOffsets;
import prok.NcrnaFixedEndpointTrials;
import prok.NcrnaPositionalBoundaryFeatures;
import prok.NcrnaPositionalBoundaryRefiner;
import prok.NcrnaSignedEndpointErrors;
import prok.NcrnaSignedEndpointErrors.Distribution;
import prok.NcrnaScavenger;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.NcrnaWeightedEndpointTable;
import prok.Orf;
import prok.PackagedPositionalBoundaryConfig;
import prok.ProkObject;
import prok.RrnaPositionalKmerTable.Table;
import shared.Shared;
import structures.ByteBuilder;

/** Fixed-vector parity through the actual family package loader and refiner.
 * Borrows two existing 69/22 exports; it trains nothing. Synthetic counted tables
 * and the installed SAM-I core exercise framing, raw features/scores, strict
 * decisions and physical holds. Scalar scores must match exactly; SIMD scores
 * use the reviewed absolute tolerance, with identical decisions in both modes.
 * These are finite plumbing checks, not accuracy.
 * Optional realroot= loads a group package; realfamily= and realcount= bind its
 * registered family and retained teaching count. Defaults retain the main case.
 * CLI follows PackagedPositionalBoundaryConfigTest's native fixture template.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaFullErrorServingTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaFullErrorServingTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			String five=null, three=null, realRoot=null, fiveGate=null, threeGate=null;
			String realFamily=DEFAULT_REAL_FAMILY;
			long realCount=DEFAULT_REAL_COUNT;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				if(a.equals("net5")){five=b;}
				else if(a.equals("net3")){three=b;}
				else if(a.equals("realroot")){realRoot=b;}
				else if(a.equals("realfamily")){realFamily=b;}
				else if(a.equals("realcount")){realCount=Parse.parseKMG(b);}
				else if(a.equals("gate5")){fiveGate=b;}
				else if(a.equals("gate3")){threeGate=b;}
				else{check(parser.parse(arg, a, b), "Unknown full-error serving fixture option: "+arg);}
			}
			check(realRoot!=null || (DEFAULT_REAL_FAMILY.equals(realFamily) && realCount==DEFAULT_REAL_COUNT),
				"realfamily= and realcount= require realroot= so requested package bindings are exercised");
			check(realFamily!=null && realCount>0 && five!=null && three!=null && parser.out1!=null && !parser.overwrite && !parser.append
				&& parser.in1==null && parser.in2==null && parser.maxReads<0 && parser.samplerate==1
				&& new File(parser.out1).mkdir(), "Supply net5=existing_export net3=existing_export out=fresh_directory; no training");
			check((fiveGate==null)==(threeGate==null), "Gate parity requires both existing32/1exports");
			run(new File(parser.out1), five, three, fiveGate, threeGate);
			if(fiveGate!=null){
				pp.outstream.println("NCRNA_GATE_SERVING_PARITY_PASS scalarCases="+gateScalarCases+" simdCases="+gateSimdCases
					+" startGateChecks="+gateChecks[0]+" stopGateChecks="+gateChecks[1]
					+" maxGatePackagedDirectDifference="+maxGatePairDifference+" rawMetadataExact=true decisionsIdentical=true");
			}
			if(realRoot!=null){
				realCases(realRoot, realFamily, realCount, five, three);
				pp.outstream.println("REAL_SAVED_TABLES_PAIR_PACKAGE_CHECK_PASS_SYNTHETIC_SPANS");
			}
			pp.outstream.println("SERVING_ARITHMETIC_PARITY_PASS scalarCases="+scalarCases+" simdCases="+simdCases
				+" maxSimdPackagedDirectDifference="+maxSimdPairDifference+" maxScalarSimdDifference="+maxScalarSimdDifference
				+" scalarExact=true decisionsIdentical=true simdAbsoluteTolerance="+SIMD_TOLERANCE);
			pp.outstream.println("NCRNA_FULL_ERROR_SERVING_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Checks actual released payload bindings on bounded synthetic candidate spans.
	 * This does not replay retained cobalamin evaluation rows or measure accuracy.
	 */
	private static void realCases(final String root, final String family, final long teachingCount,
			final String fivePath, final String threePath){
		final NcrnaFamilyRegistry.Entry entry=NcrnaFamilyRegistry.get(family);
		final String member=entry.member==null || entry.member.length()==0 ? entry.id : entry.member;
		final NcrnaFamilyLoader.Options options=new NcrnaFamilyLoader.Options();
		options.resources.root=root;
		options.endpoints=false;
		final NcrnaFamily loaded=NcrnaFamilyLoader.load(family, options, new FamilyResourceResolver());
		check(loaded.boundaryRefinerFactory instanceof NcrnaPositionalBoundaryRefiner.Resources,
			"Real package must select the actual positional resource factory for "+family);
		final NcrnaPositionalBoundaryRefiner.Resources resources=(NcrnaPositionalBoundaryRefiner.Resources)loaded.boundaryRefinerFactory;
		check(resources.fullErrors() && resources.margin()==0 && resources.radius()==10,
			"Real payload check requires the explicit full-error radius-ten zero-margin package");
		final Table[][] geometry=new Table[2][3];
		final NcrnaWeightedEndpointTable[][] weighted=new NcrnaWeightedEndpointTable[2][3];
		for(int end=0; end<2; end++){
			for(int probe=0; probe<3; probe++){
				weighted[end][probe]=NcrnaWeightedEndpointTable.loadScores(
					new File(root, member+"/"+(end==0 ? "five" : "three")+probe+".tsv").getPath(), member);
				geometry[end][probe]=weighted[end][probe].geometry;
				check(geometry[end][probe].counts==null && geometry[end][probe].records==teachingCount,
					"Real serving scores must bind the declared teaching population without training rows: family="+family
					+", expected="+teachingCount+", loaded="+geometry[end][probe].records);
			}
		}
		final NcrnaPositionalBoundaryFeatures expected=new NcrnaPositionalBoundaryFeatures(
			geometry[0], geometry[1], weighted[0], weighted[1], false);
		final CellNet[] nets={CellNetParser.load(fivePath), CellNetParser.load(threePath)};
		final NcrnaPositionalBoundaryRefiner worker=(NcrnaPositionalBoundaryRefiner)resources.create(loaded);
		for(final int model:new int[]{0, loaded.library.length-1}){
			final byte[] sequence=bases(loaded.library[model]);
			for(final int shift:new int[]{-3, 0, 4}){
				for(int fixed=0; fixed<=3; fixed++){
					parity(worker, loaded, expected, nets, sequence, model, PAD+shift,
						sequence.length-PAD-1-shift, fixed, 0, resources.oorMargin(), resources.usesOorHead(),
						resources.modelOffsets()==null ? 0 : resources.modelOffsets().median(0, model),
						resources.modelOffsets()==null ? 0 : resources.modelOffsets().median(1, model));
				}
			}
		}
	}

	private static void run(final File root, final String fivePath, final String threePath,
			final String fiveGatePath, final String threeGatePath){
		capturedScoreCases();
		readerCases(root);
		gateConfigurationCases(root);
		final FamilyResourceResolver installed=new FamilyResourceResolver();
		final String[] core={installed.resolve("sam1", "consensus", null).path,
			installed.resolve("sam1", "models", null).path, installed.resolve("sam1", "seeds", null).path};
		final NcrnaFamily baseline=NcrnaFamilyLoader.load("sam1",
			PackagedPositionalBoundaryConfigTest.options(core, null, false), new FamilyResourceResolver());
		final File payload=new File(root, "payload");
		check(payload.mkdir(), "Serving fixture payload must be fresh");
		ReadWrite.copyFile(fivePath, path(payload, "five.bbnet"));
		ReadWrite.copyFile(threePath, path(payload, "three.bbnet"));
		final String hbmName=new File(core[1]).getName();
		ReadWrite.copyFile(core[1], path(payload, hbmName));
		PackagedPositionalBoundaryConfigTest.witness(new File(payload, "training.fa"), baseline, false, false);
		final Table[][] geometry=new Table[2][3];
		final NcrnaWeightedEndpointTable[][] weighted=new NcrnaWeightedEndpointTable[2][3];
		final byte[] teaching=bases(baseline.library[0]);
		final Distribution distribution=new Distribution(new int[]{-17, 0, 4, 22}, new long[]{2, 4, 1, 3});
		final NcrnaSignedEndpointErrors errors=new NcrnaSignedEndpointErrors("sam1", distribution, distribution);
		final String[] files={"five.bbnet", "three.bbnet", hbmName, "training.fa",
			"five0.tsv", "five1.tsv", "five2.tsv", "three0.tsv", "three1.tsv", "three2.tsv"};
		for(int end=0; end<2; end++){
			for(int probe=0; probe<3; probe++){
				final int offset=probe==0 ? (end==0 ? -2 : -7) : probe==1 ? -9 : 1;
				final Table table=new Table("sam1", end==0 ? "5prime" : "3prime", 9, offset, .5, 10);
				final NcrnaWeightedEndpointTable background=new NcrnaWeightedEndpointTable(table, errors);
				for(int locus=0; locus<3; locus++){
					final byte[] sequence=teaching.clone();
					sequence[PAD+locus]=(byte)(locus==0 ? 'A' : locus==1 ? 'C' : 'G');
					final int truth=end==0 ? PAD : sequence.length-PAD-1;
					table.add(sequence, truth);
					background.add(sequence, truth);
				}
				table.finish();
				background.finish();
				geometry[end][probe]=table;
				weighted[end][probe]=background;
				final String saved=path(payload, (end==0 ? "five" : "three")+probe+".tsv");
				background.writeScores(saved);
				final NcrnaWeightedEndpointTable loaded=NcrnaWeightedEndpointTable.loadScores(saved, "sam1");
				check(loaded.geometry.counts==null, "Serving must not allocate positional training counts");
				fails(()->loaded.leaveOneOut());
				for(int word=0; word<table.space; word++){
					check(Double.doubleToRawLongBits(loaded.scoreKey(word))==Double.doubleToRawLongBits(background.scoreKey(word)),
						"Every saved weighted double must survive the score-only round trip exactly");
				}
			}
		}
		final NcrnaPositionalBoundaryFeatures expected=new NcrnaPositionalBoundaryFeatures(
			geometry[0], geometry[1], weighted[0], weighted[1], false);
		final CellNet[] nets={CellNetParser.load(fivePath), CellNetParser.load(threePath)};
		for(int policy=0; policy<4; policy++){
			final float margin=policy==1 ? .1f : 0;
			final float oorMargin=policy<2 ? margin : policy==2 ? -.25f : 0;
			final boolean useOorHead=policy!=3;
			final String name=policy==0 ? "margin0" : policy==1 ? "margin01" : policy==2 ? "signed" : "ignore";
			final String config=config(margin, hbmName)+(policy<2 ? "" : "oormargin="+(policy==2 ? "-0.25" : "ignore")+"\n");
			final File packaged=PackagedPositionalBoundaryConfigTest.tripleRoot(root, payload,
				name, config, files);
			final NcrnaFamily loaded=NcrnaFamilyLoader.load("sam1",
				PackagedPositionalBoundaryConfigTest.options(core, packaged, null), new FamilyResourceResolver());
			final NcrnaPositionalBoundaryRefiner.Resources resources=(NcrnaPositionalBoundaryRefiner.Resources)loaded.boundaryRefinerFactory;
			check(resources.fullErrors() && resources.radius()==10 && resources.margin()==margin
				&& resources.oorMargin()==oorMargin && resources.usesOorHead()==useOorHead,
				"Actual package loading must install the explicit weighted 22-output route");
			final ArrayList<NcrnaFamily> families=new ArrayList<NcrnaFamily>();
			families.add(loaded);
			fails(()->NcrnaBoundaryTrace.validateTarget("sam1", families));
			NcrnaBoundaryTrace.validateTarget("sam1", families, true);
			check(!loaded.boundaryOnRawEndpoints, "Omitting rawendpoints preserves the historical family default");
			if(policy==0){
				for(final boolean raw:new boolean[]{false, true}){
					final File rawRoot=PackagedPositionalBoundaryConfigTest.tripleRoot(root, payload,
						"raw_"+raw, config+"rawendpoints="+raw+"\n", files);
					final NcrnaFamily rawFamily=NcrnaFamilyLoader.load("sam1",
						PackagedPositionalBoundaryConfigTest.options(core, rawRoot, null), new FamilyResourceResolver());
					check(rawFamily.boundaryOnRawEndpoints==raw && !loaded.boundaryOnRawEndpoints,
						"The explicit raw-end choice applies only to its selected loaded family");
					invocationCases(root, rawFamily, raw);
				}
			}
			final NcrnaPositionalBoundaryRefiner worker=(NcrnaPositionalBoundaryRefiner)resources.create(loaded);
			check(worker.afterAlignmentTrim(), "Full-error inference belongs to the real post-trim seam");
			for(final int model:new int[]{0, loaded.library.length-1}){
				final byte[] sequence=bases(loaded.library[model]);
				for(final int shift:new int[]{-3, 0, 4}){
					for(int fixed=0; fixed<=3; fixed++){
						parity(worker, loaded, expected, nets, sequence, model, PAD+shift,
							sequence.length-PAD-1-shift, fixed, margin, oorMargin, useOorHead);
					}
				}
				sequence[PAD]='N';
				parity(worker, loaded, expected, nets, sequence, model, PAD, sequence.length-PAD-1, 0, margin, oorMargin, useOorHead);
			}
			final NcrnaFamily off=NcrnaFamilyLoader.load("sam1",
				PackagedPositionalBoundaryConfigTest.options(core, packaged, false), new FamilyResourceResolver());
			check(off.boundaryRefinerFactory==null && off.boundary5NetTemplate==null, "Explicit OFF installs no learned route");
			final File wrong=PackagedPositionalBoundaryConfigTest.tripleRoot(root, payload,
				"wrong_"+name, config.replace("table5left=five1.tsv", "table5left=five2.tsv"), files);
			fails(()->NcrnaFamilyLoader.load("sam1", PackagedPositionalBoundaryConfigTest.options(core, wrong, null), new FamilyResourceResolver()));
			final NcrnaFamilyLoader.Options overridden=PackagedPositionalBoundaryConfigTest.options(core, packaged, null);
			overridden.boundaryMarginStart=overridden.boundaryMarginStop=.125f;
			overridden.boundaryMarginStartSet=overridden.boundaryMarginStopSet=true;
			final NcrnaPositionalBoundaryRefiner.Resources changed=(NcrnaPositionalBoundaryRefiner.Resources)
				NcrnaFamilyLoader.load("sam1", overridden, new FamilyResourceResolver()).boundaryRefinerFactory;
			check(changed.margin()==.125f && changed.oorMargin()==(policy<2 ? .125f : oorMargin)
				&& changed.usesOorHead()==useOorHead,
				"A CLI positional-margin override carries the old head default but preserves an explicit signed head margin");
		}
		for(final String value:new String[]{"NaN", "Infinity", "-Infinity", "wrong"}){
			final File bad=new File(root, "bad_head_"+value+".config");
			PackagedPositionalBoundaryConfigTest.write(bad, config(0, hbmName)+"oormargin="+value+"\n");
			fails(()->PackagedPositionalBoundaryConfig.read(bad.getPath()));
		}
		final NcrnaPositionalBoundaryRefiner forced=(NcrnaPositionalBoundaryRefiner)new NcrnaPositionalBoundaryRefiner.Resources(
			baseline, geometry[0], geometry[1], weighted[0], weighted[1], new FixedNet(20, 2), new FixedNet(0, 2),
			0, 0, false).create(baseline);
		final int[] forcedMask={0};
		forced.setScoreObserver((end, input, scores, selected)->{
			check(selected==(end==0 ? 20 : 0), "The forced joint test must propose both literal inward shifts");
			forcedMask[0]|=1<<end;
		});
		final Orf crossed=candidate(teaching, baseline, 0, PAD, teaching.length-PAD-1);
		final int length=crossed.stop-crossed.start+1;
		check(!forced.refine(crossed, teaching, 0, .37f, .8f, true, 0, length, length)
			&& crossed.start==PAD && crossed.stop==teaching.length-PAD-1 && forcedMask[0]==3,
			"A literal shrink of20bases violates the allowed joint length and must retain both raw ends");
		if(fiveGatePath!=null){
			gateCases(root, payload, core, files, baseline, geometry, weighted, fivePath, threePath, fiveGatePath, threeGatePath);
		}
		centeredCases(root, payload, core, files, baseline, geometry, weighted, fivePath, threePath);
	}

	/** Real caller seam with unusable sequence: invoked refiners make no neural move.
	 * A failed trim, both fixed ends and candidate/path rejection have distinct counts. */
	private static void invocationCases(final File root, final NcrnaFamily family, final boolean raw){
		final String path=new File(root, "invocations_"+raw+".tsv").getPath();
		final byte[] sequence=new byte[240];
		Arrays.fill(sequence, (byte)'N');
		sequence[200]='A';// Actual contig GC is defined; both endpoint word windows remain unusable.
		try(Output output=new Output(path)){
			output.writer.print(NcrnaBoundaryTrace.INVOCATION_HEADER);
			final NcrnaBoundaryTrace.Writer writer=new NcrnaBoundaryTrace.Writer(output.writer, "sam1", true);
			writer.input(0, "literal_unusable.fa");
			final NcrnaBoundaryTrace.Worker owner=writer.worker(0, 0);
			final NcrnaBoundaryTrace.Family trace=owner.family(family);
			final NcrnaScavenger caller=new NcrnaScavenger(family.library, family.models, family.modelNames,
				family.kmerSet, family.kLong, family.minLen, family.windowPad, family.indexK, family.indexTopN,
				family.adaptive, family.adaptFloor, family.adaptTopFrac, family.adaptQFrac, family.fixedMinHits,
				family.scoreA, family.scoreB, family.idPass, family.idBorderline){
				@Override
				public boolean trimToAlignmentExtent(final Orf orf, final byte[] bases, final int model,
						final int start, final int stop, final float identity){return false;}
			};
			GeneCaller.applyNcrnaFamilyControls(caller, family);
			caller.setBoundaryTrace(trace);
			owner.begin("fixture", sequence.length, 0);
			try{
				final ArrayList<Orf> selected=new ArrayList<Orf>();
				for(int trial=0; trial<5; trial++){
					final int fixed=trial==3 ? 3 : trial==4 ? 1 : 0;
					final int start=fixed==0 ? 80 : 0, stop=fixed==3 ? 239 : 139;
					final Orf call=candidate(sequence, family, 0, start, stop);
					call.ncrnaPartial=fixed;
					caller.finishAcceptedBoundary(call, sequence, 0, 0, 239, 0, .8f, stop-start+1,
						.8f, true, null, null, fixed, null);
					check(call.start==start && call.stop==stop && call.orfScore==42,
						"Unusable sequence and fixed ends must preserve literal coordinates and detection score");
					if(trial!=2){trace.returned(call);}
					if(trial!=1 && trial!=2){selected.add(call);}
				}
				owner.finish(selected);
			}finally{owner.clear();}
		}
		int candidates=0, invoked=0, retained=0, pathRejected=0, scavengerRejected=0;
		final LineParser1 fields=new LineParser1('\t');
		try(Input input=new Input(path)){
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				if(line.length==0 || line[0]=='#'){continue;}
				fields.set(line);
				check(!fields.termEquals("S", 0), "Invocation-only observation must never export score-replay rows");
				if(!fields.termEquals("C", 0)){continue;}
				check(fields.terms()==37 && fields.parseInt(14)==1 && fields.parseInt(30)==0,
					"Every literal candidate has known invocation provenance without copied score buffers");
				candidates++;
				invoked+=fields.parseInt(16);
				retained+=fields.parseInt(16)*fields.parseInt(32);
				if(fields.termEquals("PATH_REJECTED", 33)){pathRejected++;}
				if(fields.termEquals("SCAVENGER_REJECTED", 33)){scavengerRejected++;}
			}
		}
		check(candidates==5 && invoked==(raw ? 4 : 1) && retained==(raw ? 2 : 1)
			&& pathRejected==1 && scavengerRejected==1,
			"Raw ON invokes four candidates and retains two; OFF invokes/retains only the one-fixed-end case; fully fixed never invokes");
	}

	/** Supported signed offsets precede both fixed feature vectors and every NN hold.
	 * Literal offsets are independent of the production resource's center method.
	 */
	private static void centeredCases(final File root, final File payload, final String[] core, final String[] files,
			final NcrnaFamily baseline, final Table[][] geometry, final NcrnaWeightedEndpointTable[][] weighted,
			final String fivePath, final String threePath){
		check(baseline.library.length>1, "The centering fixture needs distinct supported and unsupported model rows");
		final ByteBuilder rows=new ByteBuilder(NcrnaFixedEndpointTrials.OFFSET_HEADER).nl();
		for(int model=0; model<baseline.modelNames.length; model++){
			rows.append("sam1\t").append(baseline.modelNames[model].split(" ", 2)[0]);
			rows.append(model==0 ? "\t10\t10\t5\t-3\n" : "\t9\t10\t0\t0\n");
		}
		final File offsets=new File(payload, "model_offsets.tsv");
		PackagedPositionalBoundaryConfigTest.write(offsets, rows.toString());
		final String[] payloads=Arrays.copyOf(files, files.length+1);
		payloads[files.length]=offsets.getName();
		final NcrnaPositionalBoundaryFeatures expected=new NcrnaPositionalBoundaryFeatures(
			geometry[0], geometry[1], weighted[0], weighted[1], false);
		final CellNet[] direct={CellNetParser.load(fivePath), CellNetParser.load(threePath)};
		for(final boolean ignored:new boolean[]{false, true}){
			final String recipe=config(0, new File(core[1]).getName())+"modeloffsets=model_offsets.tsv\n"
				+(ignored ? "oormargin=ignore\n" : "");
			final File packaged=PackagedPositionalBoundaryConfigTest.tripleRoot(root, payload,
				ignored ? "centered_ignore" : "centered_margin0", recipe, payloads);
			final NcrnaFamily family=NcrnaFamilyLoader.load("sam1",
				PackagedPositionalBoundaryConfigTest.options(core, packaged, null), new FamilyResourceResolver());
			final NcrnaPositionalBoundaryRefiner worker=(NcrnaPositionalBoundaryRefiner)family.boundaryRefinerFactory.create(family);
			for(final int model:new int[]{0, family.library.length-1}){
				final byte[] sequence=bases(family.library[model]);
				for(final int shift:new int[]{-3, 0, 4}){
					for(int fixed=0; fixed<=3; fixed++){
						parity(worker, family, expected, direct, sequence, model, PAD+shift,
							sequence.length-PAD-1-shift, fixed, 0, 0, !ignored,
							model==0 ? 5 : 0, model==0 ? -3 : 0);
					}
				}
			}
			final File absent=PackagedPositionalBoundaryConfigTest.tripleRoot(root, payload,
				ignored ? "centered_absent_ignore" : "centered_absent_margin0",
				recipe.replace("modeloffsets=model_offsets.tsv", "modeloffsets=absent.tsv"), files);
			check(NcrnaFamilyLoader.load("sam1", PackagedPositionalBoundaryConfigTest.options(core, absent, false),
				new FamilyResourceResolver()).boundaryRefinerFactory==null, "Explicit OFF must not open model-offset payloads");
			fails(()->NcrnaFamilyLoader.load("sam1", PackagedPositionalBoundaryConfigTest.options(core, absent, null), new FamilyResourceResolver()));
		}
		final NcrnaEndpointModelOffsets bound=NcrnaEndpointModelOffsets.load(offsets.getPath(), baseline, "sam1", 10);
		check(bound.trainingSupport(0)==10 && bound.trainingSupport(baseline.library.length-1)==9,
			"The boundary support floor is inclusive and retains unsupported model rows");
		for(int kind=0; kind<4; kind++){
			final String text=rows.toString();
			final int first=text.indexOf('\n')+1, next=text.indexOf('\n', first)+1;
			final String broken=kind==0 ? text.substring(0, first)+text.substring(next)
				: kind==1 ? text.substring(0, next)+text.substring(first, next)+text.substring(next)
				: kind==2 ? text.replace("\t9\t10\t0\t0", "\t9\t10\t1\t0")
				: text.replace("\t10\t10\t5\t-3", "\t10\t20\t5\t-3");
			final File bad=new File(root, "bad_model_offsets_"+kind+".tsv");
			PackagedPositionalBoundaryConfigTest.write(bad, broken);
			fails(()->NcrnaEndpointModelOffsets.load(bad.getPath(), baseline, "sam1", 10));
		}
		final File gated=new File(root, "uncentered_old_gates.config");
		final String gateRecipe=config(0, new File(core[1]).getName())+"oormargin=ignore\n"
			+"gate5=gate5.bbnet\ngate3=gate3.bbnet\ngatecutoff5=0.1\ngatecutoff3=0.1\n"
			+"gatefeatures="+NcrnaEndpointGateFeatures.FEATURE_ORDER+"\n";
		PackagedPositionalBoundaryConfigTest.write(gated, gateRecipe);
		PackagedPositionalBoundaryConfig.read(gated.getPath());
		final File centeredGated=new File(root, "centered_old_gates.config");
		PackagedPositionalBoundaryConfigTest.write(centeredGated, gateRecipe+"modeloffsets=model_offsets.tsv\n");
		fails(()->PackagedPositionalBoundaryConfig.read(centeredGated.getPath()));
		literalCenteredCases(baseline, geometry, weighted, bound, expected);
	}

	/** Gate-field presence is independent of a legitimate numeric zero cutoff. */
	private static void gateConfigurationCases(final File root){
		final String gates="gate5=gate5.bbnet\ngate3=gate3.bbnet\ngatecutoff5=0\ngatecutoff3=0\n"
			+"gatefeatures="+NcrnaEndpointGateFeatures.FEATURE_ORDER+"\n";
		final String recipe=config(0, "training.hbmc")+gates;
		final File complete=new File(root, "zero_gate_cutoffs.config");
		PackagedPositionalBoundaryConfigTest.write(complete, recipe);
		PackagedPositionalBoundaryConfig.read(complete.getPath());
		final String[] fields={"gate5=gate5.bbnet\n", "gate3=gate3.bbnet\n", "gatecutoff5=0\n", "gatecutoff3=0\n",
			"gatefeatures="+NcrnaEndpointGateFeatures.FEATURE_ORDER+"\n"};
		for(int i=0; i<fields.length; i++){
			final File incomplete=new File(root, "incomplete_gate_"+i+".config");
			PackagedPositionalBoundaryConfigTest.write(incomplete, recipe.replace(fields[i], ""));
			fails(()->PackagedPositionalBoundaryConfig.read(incomplete.getPath()));
		}
		final File alone=new File(root, "cutoffs_without_gates.config");
		PackagedPositionalBoundaryConfigTest.write(alone, config(0, "training.hbmc")+"gatecutoff5=0\ngatecutoff3=0\n");
		fails(()->PackagedPositionalBoundaryConfig.read(alone.getPath()));
		final File comparedHead=new File(root, "gate_with_numeric_head.config");
		PackagedPositionalBoundaryConfigTest.write(comparedHead, recipe+"oormargin=0\n");
		fails(()->PackagedPositionalBoundaryConfig.read(comparedHead.getPath()));
	}

	/** Literal geometry distinguishes an invalid center from an invalid NN pair.
	 * The latter holds the valid center, never the earlier uncentered coordinates.
	 */
	private static void literalCenteredCases(final NcrnaFamily family, final Table[][] geometry,
			final NcrnaWeightedEndpointTable[][] weighted, final NcrnaEndpointModelOffsets offsets,
			final NcrnaPositionalBoundaryFeatures expected){
		final byte[] sequence=bases(family.library[0]);
		final int start=PAD, stop=sequence.length-PAD-1, centeredStart=start-5, centeredStop=stop+3;
		final NcrnaPositionalBoundaryRefiner worker=(NcrnaPositionalBoundaryRefiner)new NcrnaPositionalBoundaryRefiner.Resources(
			family, geometry[0], geometry[1], weighted[0], weighted[1], new FixedNet(20, 2), new FixedNet(0, 2),
			0, 0, false, null, null, 0, 0, offsets).create(family);
		final int[] observed={0};
		worker.setScoreObserver((end, input, scores, selected)->observed[0]|=1<<end);
		final Orf outside=candidate(sequence, family, 0, 2, stop);
		check(!worker.refine(outside, sequence, 0, .37f, .8f, true, 0, 16, sequence.length)
			&& outside.start==2 && outside.stop==stop && observed[0]==0,
			"A negative centered start holds both original ends and performs no inference");
		final int centeredLength=centeredStop-centeredStart+1;
		for(final int[] limits:new int[][]{{centeredLength+1, sequence.length}, {16, centeredLength-1}}){
			final Orf held=candidate(sequence, family, 0, start, stop);
			check(!worker.refine(held, sequence, 0, .37f, .8f, true, 0, limits[0], limits[1])
				&& held.start==start && held.stop==stop && observed[0]==0,
				"Centering must use both minimum and maximum lengths supplied by the caller, before inference");
		}
		final Orf joint=candidate(sequence, family, 0, start, stop);
		check(worker.refine(joint, sequence, 0, .37f, .8f, true, 0, centeredLength-15, sequence.length)
			&& joint.start==centeredStart && joint.stop==centeredStop && observed[0]==3 && joint.orfScore==42,
			"Two valid individual ten-base inward moves make a jointly short span and must retain the valid center");
		observed[0]=0;
		final NcrnaPositionalBoundaryRefiner stable=(NcrnaPositionalBoundaryRefiner)new NcrnaPositionalBoundaryRefiner.Resources(
			family, geometry[0], geometry[1], weighted[0], weighted[1], new FixedNet(10, 2), new FixedNet(10, 2),
			0, 0, false, null, null, 0, 0, offsets).create(family);
		stable.setScoreObserver((end, input, scores, selected)->observed[0]|=1<<end);
		final Orf fixed=candidate(sequence, family, 0, 0, stop);
		stable.refine(fixed, sequence, 0, .37f, .8f, true, 1, 16, sequence.length);
		check(fixed.start==0 && fixed.stop==centeredStop && observed[0]==2,
			"A physical fixed start remains zero while the free end is centered and inferred");
		sequence[centeredStart]='N';
		check(expected.fill(sequence, centeredStart, centeredStop, family.library[0], family.models[0],
			plateau(family.models[0]), .37f, 0)==2, "A start-probe ambiguity makes only the start feature unavailable");
		observed[0]=0;
		final Orf missing=candidate(sequence, family, 0, start, stop);
		stable.refine(missing, sequence, 0, .37f, .8f, true, 0, 16, sequence.length);
		check(missing.start==centeredStart && missing.stop==centeredStop && observed[0]==2,
			"A valid center remains in place when one end lacks features; only the available end is inferred");
	}

	/** Actual gate exports traverse the package loader; deterministic nets exercise equality and holds. */
	private static void gateCases(final File root, final File payload, final String[] core, final String[] files,
			final NcrnaFamily baseline, final Table[][] geometry, final NcrnaWeightedEndpointTable[][] weighted,
			final String fivePath, final String threePath, final String fiveGatePath, final String threeGatePath){
		ReadWrite.copyFile(fiveGatePath, path(payload, "gate5.bbnet"));
		ReadWrite.copyFile(threeGatePath, path(payload, "gate3.bbnet"));
		final String[] gateFiles=Arrays.copyOf(files, files.length+2);
		gateFiles[files.length]="gate5.bbnet";
		gateFiles[files.length+1]="gate3.bbnet";
		final CellNet[] direct={CellNetParser.load(fivePath), CellNetParser.load(threePath)};
		final CellNet[] gates={CellNetParser.load(fiveGatePath), CellNetParser.load(threeGatePath)};
		final NcrnaPositionalBoundaryFeatures expected=new NcrnaPositionalBoundaryFeatures(
			geometry[0], geometry[1], weighted[0], weighted[1], false);
		for(int policy=0; policy<2; policy++){
			final float[] cutoffs=policy==0 ? new float[]{.1f, .1f} : new float[]{.63589126f, .4522214f};
			final String config=config(0, new File(core[1]).getName())+"oormargin=ignore\n"
				+"gate5=gate5.bbnet\ngate3=gate3.bbnet\ngatecutoff5="+cutoffs[0]+"\ngatecutoff3="+cutoffs[1]
				+"\ngatefeatures="+NcrnaEndpointGateFeatures.FEATURE_ORDER+"\n";
			final File packaged=PackagedPositionalBoundaryConfigTest.tripleRoot(root, payload, "gate_policy"+policy, config, gateFiles);
			final NcrnaFamily loaded=NcrnaFamilyLoader.load("sam1",
				PackagedPositionalBoundaryConfigTest.options(core, packaged, null), new FamilyResourceResolver());
			final NcrnaPositionalBoundaryRefiner.Resources resources=(NcrnaPositionalBoundaryRefiner.Resources)loaded.boundaryRefinerFactory;
			check(resources.gated() && resources.oorMargin()==0 && !resources.usesOorHead(),
				"The actual package must replace the head veto with both configured movement gates");
			final NcrnaPositionalBoundaryRefiner worker=(NcrnaPositionalBoundaryRefiner)resources.create(loaded);
			for(final int model:new int[]{0, loaded.library.length-1}){
				final byte[] sequence=bases(loaded.library[model]);
				for(final int shift:new int[]{-3, 0, 4}){
					for(int fixed=0; fixed<=3; fixed++){
						gateParity(worker, loaded, expected, direct, gates, cutoffs, sequence, model,
							PAD+shift, sequence.length-PAD-1-shift, fixed);
					}
				}
				sequence[PAD]='N';
				gateParity(worker, loaded, expected, direct, gates, cutoffs, sequence, model,
					PAD, sequence.length-PAD-1, 0);
			}
			// Literal non-current scores guarantee that the REAL gate exports are exercised at both ends.
			final CellNet[] literal={new FixedNet(12, 2), new FixedNet(8, 2)};
			final NcrnaPositionalBoundaryRefiner forced=(NcrnaPositionalBoundaryRefiner)new NcrnaPositionalBoundaryRefiner.Resources(
				loaded, geometry[0], geometry[1], weighted[0], weighted[1], literal[0], literal[1], 0,
				0, false, gates[0], gates[1], cutoffs[0], cutoffs[1]).create(loaded);
			final byte[] literalBases=bases(loaded.library[0]);
			gateParity(forced, loaded, expected, literal, gates, cutoffs, literalBases, 0, PAD, literalBases.length-PAD-1, 0);
			final NcrnaFamily off=NcrnaFamilyLoader.load("sam1",
				PackagedPositionalBoundaryConfigTest.options(core, packaged, false), new FamilyResourceResolver());
			check(off.boundaryRefinerFactory==null, "Explicit OFF bypasses the complete gate route");
			final File wrong=PackagedPositionalBoundaryConfigTest.tripleRoot(root, payload, "gate_wrong_width"+policy,
				config.replace("gate5=gate5.bbnet", "gate5=five.bbnet"), gateFiles);
			fails(()->NcrnaFamilyLoader.load("sam1", PackagedPositionalBoundaryConfigTest.options(core, wrong, null), new FamilyResourceResolver()));
			final NcrnaFamilyLoader.Options override=PackagedPositionalBoundaryConfigTest.options(core, packaged, null);
			override.boundaryMarginStart=override.boundaryMarginStop=.1f;
			override.boundaryMarginStartSet=override.boundaryMarginStopSet=true;
			fails(()->NcrnaFamilyLoader.load("sam1", override, new FamilyResourceResolver()));
			final File missing=PackagedPositionalBoundaryConfigTest.tripleRoot(root, payload, "gate_missing_payload"+policy,
				config.replace("gate5.bbnet", "absent5.bbnet").replace("gate3.bbnet", "absent3.bbnet"), files);
			check(NcrnaFamilyLoader.load("sam1", PackagedPositionalBoundaryConfigTest.options(core, missing, false),
				new FamilyResourceResolver()).boundaryRefinerFactory==null, "OFF must not open missing gate payloads");
			fails(()->NcrnaFamilyLoader.load("sam1", PackagedPositionalBoundaryConfigTest.options(core, missing, null), new FamilyResourceResolver()));
		}
		check(gateChecks[0]>0 && gateChecks[1]>0, "Both real gate exports must have actual observed fixed-vector checks");
		literalGateCases(baseline, geometry, weighted);
	}

	/** Independent raw32input and score oracle; scalar exactness and SIMD decision parity are separate. */
	private static void gateParity(final NcrnaPositionalBoundaryRefiner worker, final NcrnaFamily family,
			final NcrnaPositionalBoundaryFeatures expected, final CellNet[] nets, final CellNet[] gates,
			final float[] cutoffs, final byte[] bases, final int model, final int start, final int stop, final int fixed){
		final boolean saved=Shared.SIMD;
		int scalarStart=-1, scalarStop=-1;
		try{
			for(int mode=0; mode<2; mode++){
				if(mode==1 && (!saved || !simd.Vector.simd128)){continue;}
				Shared.SIMD=mode==1;
				final boolean vector=mode==1;
				final int usable=expected.fill(bases, start, stop, family.library[model], family.models[model], plateau(family.models[model]), .37f, fixed);
				final int[] choices={10, 10}, proposals={-1, -1}, observed={0, 0};
				final float[][] scores=new float[2][], inputs=new float[2][];
				int gateMask=0;
				for(int end=0; end<2; end++){
					if((usable&(1<<end))==0){continue;}
					nets[end].applyInput(expected.values[end]);
					nets[end].feedForward();
					scores[end]=nets[end].getOutput();
					proposals[end]=literalProposal(scores[end]);
					if(proposals[end]<0){continue;}
					gateMask|=1<<end;
					final int shift=proposals[end]-10;
					final long length=(long)stop-start+1+(end==0 ? -shift : shift);
					final float[] raw=expected.values[end], input=new float[32];
					System.arraycopy(scores[end], 0, input, 0, 22);
					input[22]=raw[63];
					input[23]=length/(float)family.library[model].length;
					input[24]=raw[64];
					input[25]=literalGC(bases, start, stop);
					input[26]=literalGC(family.library[model], 0, family.library[model].length-1);
					input[27]=input[25]-input[26];
					input[28]=2f*raw[65]-1f;
					System.arraycopy(raw, 66, input, 29, 3);
					inputs[end]=input;
				}
				worker.setGateObserver((end, input, score, proposed, accepted)->{
					check(proposed==proposals[end], "Gate proposal must match direct first-stage inference");
					for(int i=0; i<32; i++){
						if(!sameScore(input[i], inputs[end][i], vector && i<22)){
							throw new AssertionError("Gate input mismatch: field="+i+", model="+model+", end="+end+", span="+start+".."+stop
								+", fixed="+fixed+", vector="+vector+", actual="+input[i]+", expected="+inputs[end][i]
								+", actualBits="+Float.floatToRawIntBits(input[i])+", expectedBits="+Float.floatToRawIntBits(inputs[end][i])
								+"; metadata and scalar scores must be exact; accelerated score limit="+SIMD_TOLERANCE
								+" because simd/SIMD.java:330 and simd/SIMD128.java:330 use FloatVector.reduceLanes(ADD), whose sum order is not fixed");
						}
					}
					// Feed the SAME actual32inputs to isolate gate-copy arithmetic from first-stage SIMD drift.
					gates[end].applyInput(input);
					gates[end].feedForward();
					final float direct=gates[end].getOutput(0);
					maxGatePairDifference=Math.max(maxGatePairDifference, Math.abs((double)score-direct));
					choices[end]=direct>=cutoffs[end] ? proposed : 10;
					check(sameScore(score, direct, vector) && accepted==(choices[end]!=10), "Gate score and inclusive decision must match same-input direct inference");
					check((observed[1]&(1<<end))==0, "Each proposed end may reach the gate exactly once");
					observed[1]|=1<<end;
					gateChecks[end]++;
					Arrays.fill(input, Float.NaN);
				});
				worker.setScoreObserver((end, input, actual, chosen)->{
					check(Arrays.equals(input, expected.values[end]) && chosen==choices[end], "Gate serving preserves all69first-stage features and its selected offset");
					for(int i=0; i<22; i++){check(sameScore(actual[i], scores[end][i], vector), "Every first-stage score must match direct inference");}
					observed[0]|=1<<end;
				});
				final Orf call=candidate(bases, family, model, start, stop);
				worker.refine(call, bases, model, .37f, .8f, true, fixed, 16, bases.length);
				check(observed[0]==usable && observed[1]==gateMask, "Fixed/unavailable/current-winning/tied ends bypass their gate");
				int five=start+choices[0]-10, three=stop+choices[1]-10;
				if(five<0 || three>=bases.length || three-five+1<16 || three-five+1>bases.length){five=start; three=stop;}
				check(call.start==five && call.stop==three && call.orfScore==42, "Gate decisions precede unchanged physical holds and caller score");
				if(vector){
					check(call.start==scalarStart && call.stop==scalarStop, "Scalar and SIMD must make identical final gate decisions");
					gateSimdCases++;
				}else{
					scalarStart=call.start;
					scalarStop=call.stop;
					gateScalarCases++;
				}
			}
		}finally{
			Shared.SIMD=saved;
			worker.setGateObserver(null);
			worker.setScoreObserver(null);
		}
	}

	/** Fixed observed pairs make the chosen accelerated-score contract reproducible
	 * without waiting for another intermittent vector evaluation. Scalar results
	 * remain exact, including when the same pair is accepted in accelerated mode. */
	private static void capturedScoreCases(){
		final int[][] bits={{1061254580, 1061254574}, {-1168271484, -1168267007},
			{-1168267007, -1168271484}};
		for(final int[] pair:bits){
			final float actual=Float.intBitsToFloat(pair[0]), expected=Float.intBitsToFloat(pair[1]);
			check(sameScore(actual, expected, true),
				"Captured accelerated scores must satisfy the selected absolute limit: actualBits="+pair[0]
				+", expectedBits="+pair[1]+", difference="+Math.abs((double)actual-expected));
			check(!sameScore(actual, expected, false), "The accelerated allowance must not relax exact scalar comparisons");
		}
		check(!sameScore(.5f, .500002f, true), "A two-millionth score error exceeds the selected accelerated limit");
		check(!sameScore(Float.NaN, 0, true) && !sameScore(0, Float.NaN, true)
			&& !sameScore(Float.POSITIVE_INFINITY, 0, true) && !sameScore(Float.POSITIVE_INFINITY, Float.POSITIVE_INFINITY, true),
			"An accelerated allowance cannot admit nonfinite scores");
	}

	/** Exact scalar bits; the established SIMD allowance does not excuse changed decisions. */
	private static boolean sameScore(final float a, final float b, final boolean vector){
		return vector ? Math.abs((double)a-b)<=SIMD_TOLERANCE : Float.floatToRawIntBits(a)==Float.floatToRawIntBits(b);
	}
	private static int literalProposal(final float[] scores){
		int best=0;
		for(int i=1; i<21; i++){if(scores[i]>scores[best]){best=i;}}
		if(best==10){return -1;}
		for(int i=0; i<21; i++){if(i!=best && scores[i]==scores[best]){return -1;}}
		return best;
	}
	/** Explicit alphabet interpretation avoids borrowing the production GC function as its own oracle. */
	private static float literalGC(final byte[] bases, final int start, final int stop){
		int count=0, gc=0;
		for(int i=start; i<=stop; i++){
			final byte base=bases[i];
			if(base=='A' || base=='C' || base=='G' || base=='T'){count++;}
			if(base=='C' || base=='G'){gc++;}
		}
		check(count>0, "The literal GC fixture must contain recognized bases");
		return gc/(float)count;
	}

	/** Literal score equality accepts, while current/tied/fixed ends never evaluate a gate. */
	private static void literalGateCases(final NcrnaFamily family, final Table[][] geometry,
			final NcrnaWeightedEndpointTable[][] weighted){
		final byte[] bases=bases(family.library[0]);
		final int start=PAD, stop=bases.length-PAD-1, length=stop-start+1;
		for(final float cutoff:new float[]{0, .1f}){
			for(final float score:new float[]{Math.nextDown(cutoff), cutoff, Math.nextUp(cutoff), Float.NaN}){
				final int[] calls={0};
				final NcrnaPositionalBoundaryRefiner worker=(NcrnaPositionalBoundaryRefiner)new NcrnaPositionalBoundaryRefiner.Resources(
					family, geometry[0], geometry[1], weighted[0], weighted[1], new FixedNet(12, 2), new FixedNet(8, 2),
					0, 0, false, new FixedGate(score, calls), new FixedGate(score, calls), cutoff, cutoff).create(family);
				final Orf call=candidate(bases, family, 0, start, stop);
				if(Float.isNaN(score)){
					fails(()->worker.refine(call, bases, 0, .37f, .8f, true, 0, 16, bases.length));
					continue;
				}
				worker.refine(call, bases, 0, .37f, .8f, true, 0, 16, bases.length);
				final boolean accepted=score>=cutoff;
				check(calls[0]==2 && call.start==(accepted ? start+2 : start) && call.stop==(accepted ? stop-2 : stop),
					"Inclusive gate equality and adjacent floats must give their literal decisions despite a larger head21");
				final Orf joint=candidate(bases, family, 0, start, stop);
				worker.refine(joint, bases, 0, .37f, .8f, true, 0, length, length);
				check(joint.start==start && joint.stop==stop && joint.orfScore==42,
					"Accepted gate proposals cannot bypass the unchanged joint-length hold or alter detection score");
			}
		}
		for(int kind=0; kind<3; kind++){
			final int[] calls={0};
			final FixedNet five=kind==0 ? new FixedNet(10, 2) : new FixedNet(12, 2, 10);
			final FixedNet three=kind==0 ? new FixedNet(10, 2) : new FixedNet(8, 2, 7);
			final NcrnaPositionalBoundaryRefiner worker=(NcrnaPositionalBoundaryRefiner)new NcrnaPositionalBoundaryRefiner.Resources(
				family, geometry[0], geometry[1], weighted[0], weighted[1], five, three,
				0, 0, false, new FixedGate(1, calls), new FixedGate(1, calls), .1f, .1f).create(family);
			final Orf call=candidate(bases, family, 0, start, stop);
			worker.refine(call, bases, 0, .37f, .8f, true, kind==2 ? 3 : 0, 16, bases.length);
			check(calls[0]==0 && call.start==start && call.stop==stop,
				"Current winners, reachable ties and fixed ends hold without an extra gate evaluation");
		}
	}

	/** Observes the real loaded worker; expected features use pre-serialization tables. */
	private static void parity(final NcrnaPositionalBoundaryRefiner worker, final NcrnaFamily family,
			final NcrnaPositionalBoundaryFeatures expected, final CellNet[] nets, final byte[] bases, final int model,
			final int start, final int stop, final int fixed, final float margin){
		parity(worker, family, expected, nets, bases, model, start, stop, fixed, margin, margin, true);
	}

	private static void parity(final NcrnaPositionalBoundaryRefiner worker, final NcrnaFamily family,
			final NcrnaPositionalBoundaryFeatures expected, final CellNet[] nets, final byte[] bases, final int model,
			final int start, final int stop, final int fixed, final float margin, final float oorMargin, final boolean useOorHead){
		parity(worker, family, expected, nets, bases, model, start, stop, fixed, margin, oorMargin, useOorHead, 0, 0);
	}

	private static void parity(final NcrnaPositionalBoundaryRefiner worker, final NcrnaFamily family,
			final NcrnaPositionalBoundaryFeatures expected, final CellNet[] nets, final byte[] bases, final int model,
			final int start, final int stop, final int fixed, final float margin, final float oorMargin, final boolean useOorHead,
			final int median5, final int median3){
		final boolean savedSimd=Shared.SIMD;
		final float[][] scalarInputs=new float[2][], scalarScores=new float[2][];
		final int[] scalarChoices={10, 10};
		try{
			Shared.SIMD=false;
			final int usable=parityMode(worker, family, expected, nets, bases, model, start, stop, fixed, margin, oorMargin, useOorHead,
				false, scalarInputs, scalarScores, scalarChoices, median5, median3);
			if(savedSimd && simd.Vector.simd128){
				Shared.SIMD=true;
				check(parityMode(worker, family, expected, nets, bases, model, start, stop, fixed, margin, oorMargin, useOorHead,
					true, scalarInputs, scalarScores, scalarChoices, median5, median3)==usable,
					"Scalar and SIMD must infer the same usable physical ends");
			}
		}finally{
			Shared.SIMD=savedSimd;
			worker.setScoreObserver(null);
		}
	}

	/** Compares the packaged worker with direct evaluation in one arithmetic mode.
	 * Scalar results also anchor exact cross-mode choices; SIMD noise cannot hide
	 * a different offset, move/hold decision, feature vector or usable-end mask.
	 */
	private static int parityMode(final NcrnaPositionalBoundaryRefiner worker, final NcrnaFamily family,
			final NcrnaPositionalBoundaryFeatures expected, final CellNet[] nets, final byte[] bases, final int model,
			final int start, final int stop, final int fixed, final float margin, final float oorMargin, final boolean useOorHead,
			final boolean vector, final float[][] scalarInputs, final float[][] scalarScores, final int[] scalarChoices,
			final int median5, final int median3){
		final long centeredStart=(long)start-((fixed&1)==0 ? median5 : 0);
		final long centeredStop=(long)stop-((fixed&2)==0 ? median3 : 0);
		if(centeredStart<0 || centeredStop>=bases.length || centeredStop-centeredStart+1<16
				|| centeredStop-centeredStart+1>bases.length){
			worker.setScoreObserver((end, input, scores, chosen)->{
				throw new AssertionError("An invalid center must not run either network");
			});
			final Orf held=candidate(bases, family, model, start, stop);
			check(!worker.refine(held, bases, model, .37f, .8f, true, fixed, 16, bases.length)
				&& held.start==start && held.stop==stop && held.orfScore==42,
				"An invalid literal centered pair retains the original coordinates and detection score");
			if(vector){simdCases++;}else{scalarCases++;}
			return 0;
		}
		final int baseStart=(int)centeredStart, baseStop=(int)centeredStop;
		final int usable=expected.fill(bases, baseStart, baseStop, family.library[model], family.models[model],
			plateau(family.models[model]), .37f, fixed);
		final int[] selected={10, 10}, observed={0};
		worker.setScoreObserver((end, input, scores, chosen)->{
			check((usable&(1<<end))!=0 && (observed[0]&(1<<end))==0, "Only usable, unfixed ends may be observed once");
			check(Arrays.equals(input, expected.values[end]), "All69raw features must match before and after package serialization");
			if(vector){
				check(Arrays.equals(input, scalarInputs[end]), "Changing NN arithmetic must preserve all raw features");
			}else{scalarInputs[end]=input.clone();}
			final CellNet net=nets[end];
			net.applyInput(expected.values[end]);
			net.feedForward();
			final float[] direct=net.getOutput();
			check(scores.length==22 && direct.length==22, "Both serving and direct evaluation must retain every full-error output");
			for(int output=0; output<scores.length; output++){
				final double difference=Math.abs((double)scores[output]-direct[output]);
				final boolean same=Float.floatToRawIntBits(scores[output])==Float.floatToRawIntBits(direct[output]);
				if(vector){maxSimdPairDifference=Math.max(maxSimdPairDifference, difference);}
				if(vector ? !(difference<=SIMD_TOLERANCE) : !same){
					throw new AssertionError("Packaged/direct raw-score mismatch: model="+model+", end="+end
						+", span="+start+".."+stop+", fixed="+fixed+", margin="+margin+", oorMargin="+oorMargin
						+", mode="+(vector ? "SIMD" : "scalar")+", allowedDifference="+(vector ? SIMD_TOLERANCE : 0)
						+", output="+output+", packaged="+scores[output]+", direct="+direct[output]
						+", packagedBits="+Float.floatToRawIntBits(scores[output])+", directBits="+Float.floatToRawIntBits(direct[output])
						+(vector ? "; the sum order is not fixed for FloatVector.reduceLanes(ADD) at simd/SIMD.java:330 and simd/SIMD128.java:330"
							: "; scalar scores must match exactly"));
				}
			}
			selected[end]=literal(scores, margin, oorMargin, useOorHead);
			final int evaluated=literal(direct, margin, oorMargin, useOorHead);
			if(chosen!=selected[end] || chosen!=evaluated || vector && chosen!=scalarChoices[end]){
				throw new AssertionError("Serving decision mismatch: model="+model+", end="+end+", span="+start+".."+stop
					+", fixed="+fixed+", margin="+margin+", oorMargin="+oorMargin+", mode="+(vector ? "SIMD" : "scalar")
					+", packaged="+chosen+", oracle="+selected[end]+", direct="+evaluated
					+", scalar="+(vector ? scalarChoices[end] : chosen));
			}
			if(vector){
				for(int output=0; output<scores.length; output++){
					maxScalarSimdDifference=Math.max(maxScalarSimdDifference, Math.abs((double)scores[output]-scalarScores[end][output]));
				}
			}else{
				scalarScores[end]=scores.clone();
				scalarChoices[end]=chosen;
			}
			observed[0]|=1<<end;
			Arrays.fill(input, Float.NaN);
			Arrays.fill(scores, Float.NaN);// Diagnostic copies cannot corrupt the worker.
		});
		final Orf call=candidate(bases, family, model, start, stop);
		worker.refine(call, bases, model, .37f, .8f, true, fixed, 16, bases.length);
		check(observed[0]==usable, "Every and only usable end must traverse actual network inference");
		int five=baseStart+selected[0]-10, three=baseStop+selected[1]-10;
		if(five<0 || three>=bases.length || three-five+1<16 || three-five+1>bases.length){
			five=baseStart;
			three=baseStop;
		}
		check(call.start==five && call.stop==three && call.orfScore==42, "Selected ends, physical holds and caller score must survive the real refiner");
		worker.setScoreObserver(null);
		if(vector){simdCases++;}else{scalarCases++;}
		return usable;
	}

	private static int literal(final float[] scores, final float margin, final float oorMargin, final boolean useOorHead){
		int best=0;
		for(int i=1; i<21; i++){if(scores[i]>scores[best]){best=i;}}
		for(int i=0; i<21; i++){
			if(i!=best && (double)scores[best]<=(double)scores[i]+margin){return 10;}
		}
		if(useOorHead && (double)scores[best]<=(double)scores[21]+oorMargin){return 10;}
		return best;
	}

	/** Independent dense K5 wire fixture has histogram9loci but table4loci. */
	private static void readerCases(final File root){
		final String bins="#-2\t2\t0.4\t4\t0\t0\n#0\t4\t0.0\t4\t0\t0\n#3\t3\t0.6\t3\t1\t0\n";
		final File good=new File(root, "dense.tsv");
		writeDense(good, bins, -1, false, false);
		final NcrnaWeightedEndpointTable table=NcrnaWeightedEndpointTable.loadScores(good.getPath(), "fixture");
		for(int word=0; word<1024; word++){
			check(Double.doubleToRawLongBits(table.scoreKey(word))==Double.doubleToRawLongBits((word+1)/524288.0), "Cold decimal parsing must retain exact saved double bits");
		}
		check(table.geometry.counts==null, "Loaded geometry must be score-only");
		fails(()->table.leaveOneOut());
		fails(()->NcrnaWeightedEndpointTable.loadScores(good.getPath(), "wrong_member"));
		for(int kind=0; kind<4; kind++){
			final File bad=new File(root, "bad"+kind+".tsv");
			writeDense(bad, kind==0 ? bins.replace("0.6", "0.7") : bins, kind==1 ? 17 : -1, kind==2, kind==3);
			fails(()->NcrnaWeightedEndpointTable.loadScores(bad.getPath(), "fixture"));
		}
	}
	private static void writeDense(final File path, final String bins, final int omit, final boolean duplicate, final boolean nonfinite){
		final ByteBuilder text=new ByteBuilder("#NCRNA_SIGNED_BACKGROUND_SCORES_V1\nmember\tend\tk\tradius_bp\tprobe_offset_bp\talpha_per_key\ttable_loci\nfixture\t5prime\t5\t3\t0\t0.5\t4\n#offset_bp\tobserved_errors\tbackground_weight\tvalid\tclipped\tambiguous\n").append(bins).append("kmer\tlog_score\n");
		final byte[] alphabet={'A', 'C', 'G', 'T'};
		for(int word=0; word<1024; word++){
			if(word==omit){continue;}
			final int key=duplicate && word==1 ? 0 : word;
			for(int shift=8; shift>=0; shift-=2){text.append(alphabet[(key>>>shift)&3]);}
			text.tab().append(nonfinite && word==7 ? "NaN" : Double.toString((word+1)/524288.0)).nl();
		}
		PackagedPositionalBoundaryConfigTest.write(path, text.toString());
	}
	private static String config(final float margin, final String hbmName){
		return "schema=ncrna_positional_full_error_v1\nfeatures="+PackagedPositionalBoundaryConfig.TRIPLE_FEATURE_ORDER
			+"\nfamily=sam1\nenabled=t\nradius=10\nmargin="+margin+"\nstartinside=7\nstartoutside=2\nstopinside=8\nstopoutside=1\n"
			+"library=training.fa\ntraininghbm="+hbmName+"\nnet5=five.bbnet\nnet3=three.bbnet\n"
			+"table5=five0.tsv\ntable5left=five1.tsv\ntable5right=five2.tsv\ntable3=three0.tsv\ntable3left=three1.tsv\ntable3right=three2.tsv\n";
	}
	private static byte[] bases(final byte[] model){
		final byte[] result=new byte[model.length+2*PAD];
		final byte[] alphabet={'A', 'C', 'G', 'T'};
		for(int i=0; i<result.length; i++){result[i]=alphabet[(i*17+i/7)&3];}
		for(int i=0; i<model.length; i++){
			final byte base=model[i];
			result[PAD+i]=base=='A' || base=='C' || base=='G' || base=='T' ? base : (byte)'A';
		}
		return result;
	}
	private static Orf candidate(final byte[] bases, final NcrnaFamily family, final int model, final int start, final int stop){
		final Orf call=new Orf("fixture", start, stop, 0, 0, bases, false, ProkObject.RNA);
		call.ncrnaFamily=family.name;
		call.trnaModel=family.modelNames[model];
		call.orfScore=42;
		return call;
	}
	/** Independent coverage oracle: max ref+deletion count over the complete HBM,
	 * matching TrnaBoundaryFeatures' documented plateau rather than its midpoint.
	 */
	private static float plateau(final BaseGraph model){
		check(model!=null && model.ref.length==model.del.length, "Fixture needs the paired HBM coverage axis");
		int maximum=0;
		for(int position=0; position<model.ref.length; position++){
			maximum=Math.max(maximum, model.ref[position].countSum+model.del[position].countSum);
		}
		check(maximum>0, "The loaded fixture HBM must contain observed coverage");
		return maximum;
	}
	private static final class FixedNet extends CellNet{
		FixedNet(final int peak_, final float head_){
			this(peak_, head_, -1);
		}
		FixedNet(final int peak_, final float head_, final int tied_){
			peak=peak_;
			head=head_;
			tied=tied_;
			scores[peak]=1;
			scores[21]=head;
			if(tied>=0){scores[tied]=1;}
		}
		@Override
		public int numInputs(){return 69;}
		@Override
		public int numOutputs(){return 22;}
		@Override
		public CellNet copy(final boolean deltas){return new FixedNet(peak, head, tied);}
		@Override
		public CellNet applyInput(final float[] input){
			NcrnaFullErrorServingTest.check(input.length==69, "Literal full-error net requires69features");
			return this;
		}
		@Override
		public float feedForward(){return 1;}
		@Override
		public float[] getOutput(){return scores;}
		@Override
		public float getOutput(final int index){return scores[index];}
		private final int peak;
		private final int tied;
		private final float head;
		private final float[] scores=new float[22];
	}
	/** Shared invocation counter survives the actual worker's copy(false) boundary. */
	private static final class FixedGate extends CellNet{
		FixedGate(final float value_, final int[] calls_){value=value_; calls=calls_;}
		@Override
		public int numInputs(){return 32;}
		@Override
		public int numOutputs(){return 1;}
		@Override
		public CellNet copy(final boolean deltas){return new FixedGate(value, calls);}
		@Override
		public CellNet applyInput(final float[] input){NcrnaFullErrorServingTest.check(input.length==32, "Literal gate requires32inputs"); return this;}
		@Override
		public float feedForward(){calls[0]++; return value;}
		@Override
		public float getOutput(final int index){NcrnaFullErrorServingTest.check(index==0, "Literal gate has one output"); return value;}
		private final float value;
		private final int[] calls;
	}
	private static String path(final File root, final String file){return new File(root, file).getPath();}
	private static void fails(final Runnable action){
		try{action.run();}catch(IllegalArgumentException | IllegalStateException expected){return;}
		throw new AssertionError("Malformed or unsupported full-error state must be rejected");
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private static final String DEFAULT_REAL_FAMILY="cobalamin_rf00174";
	private static final long DEFAULT_REAL_COUNT=88390;
	private static final int PAD=80;
	/** Test-only absolute score limit selected by G11 on November19. The vector
	 * ADD reduction at simd/SIMD.java:330 and simd/SIMD128.java:330 may associate
	 * sums differently over time (JDK25 FloatVector.reduceLanes contract). This is
	 * an acceptance choice, not a mathematical bound for arbitrary networks or
	 * a production decision epsilon. Inputs, scalar scores and decisions stay exact. */
	private static final double SIMD_TOLERANCE=1e-6;
	private static long scalarCases, simdCases;
	private static long gateScalarCases, gateSimdCases;
	private static final long[] gateChecks=new long[2];
	private static double maxGatePairDifference;
	private static double maxSimdPairDifference, maxScalarSimdDifference;
}
