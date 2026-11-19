package test.prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import fileIO.ReadWrite;
import consensus.BaseGraph;
import consensus.BaseGraphHelper;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.FamilyResourceResolver;
import prok.NcrnaBoundaryRefiner;
import prok.NcrnaFamily;
import prok.NcrnaFamilyLoader;
import prok.NcrnaEndpointGateFeatures;
import prok.NcrnaCenteredGateFeatures;
import prok.NcrnaNativeTrainer;
import prok.NcrnaPositionalBoundaryRefiner;
import prok.NcrnaScavenger;
import prok.NcrnaTextIO.Output;
import prok.Orf;
import prok.PackagedBoundaryConfig;
import prok.PackagedPositionalBoundaryConfig;
import prok.ProkObject;
import prok.RrnaPositionalKmerTable;
import prok.TrnaConsensusBuilder;
import shared.Shared;
import structures.ByteBuilder;

/** Actual family-loader ON/OFF and library/geometry binding checks.
 * A tiny6→3 native fit supplies an exported finite payload at radius1; this is
 * a schema/lifecycle oracle, not another biological fit or a timing experiment.
 * The installed SAM-I core is borrowed read-only and all outputs are fresh.
 * @author Brian Bushnell, Raiden
 */
public final class PackagedPositionalBoundaryConfigTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, PackagedPositionalBoundaryConfigTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			String installation=null;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				if(a.equals("installation")){installation=b;}
				else{check(parser.parse(arg, a, b), "Unknown positional-package fixture argument: "+arg);}
			}
			check(parser.out1!=null && installation!=null && !parser.overwrite && !parser.append
				&& parser.in1==null && parser.in2==null && parser.maxReads<0 && parser.samplerate==1
				&& new File(parser.out1).mkdir(), "Supply installation= and out=fresh_fixture_directory");
			run(new File(parser.out1), new File(installation));
			pp.outstream.println("PACKAGED_POSITIONAL_BOUNDARY_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(final File root, final File installation){
		gateConfigCases(root);
		final FamilyResourceResolver installed=new FamilyResourceResolver();
		final String[] core={installed.resolve("sam1", "consensus", null).path,
			installed.resolve("sam1", "models", null).path, installed.resolve("sam1", "seeds", null).path};
		final NcrnaFamily baseline=NcrnaFamilyLoader.load("sam1", options(core, null, false), new FamilyResourceResolver());
		check(baseline.boundaryRefinerFactory==null && baseline.boundary5NetTemplate==null,
			"Absent positional configuration with explicit OFF preserves ordinary family loading");
		final File payload=new File(root, "payload");
		check(payload.mkdir(), "Finite payload root must be fresh");
		write(new File(payload, "train.tsv"), matrix());
		write(new File(payload, "monitor.tsv"), matrix());
		write(new File(payload, "recipe.config"), "dims=6,4,3\nseed=58\nepochs=2\nvfraction=0\nfinal=sigmoid\nhidden=tanh\nsimd=f\nt=1\nprintevery=1\now=t\n");
		new NcrnaNativeTrainer(installation, 1).trainPositional(path(payload, "recipe.config"),
			new String[]{path(payload, "train.tsv")}, path(payload, "monitor.tsv"), new File(payload, "net.bbnet"), 6, 3);
		final byte[] bases=new byte[240];
		Arrays.fill(bases, (byte)'A');
		for(int radius=1; radius<=2; radius++){
			table(payload, "five"+radius+".tsv", bases, "sam1", "5prime", 0, 80, radius);
			table(payload, "three"+radius+".tsv", bases, "sam1", "3prime", -8, 139, radius);
		}
		table(payload, "placed_five.tsv", bases, "sam1", "5prime", 5, -2, 80, 1);
		table(payload, "placed_three.tsv", bases, "sam1", "3prime", 8, -4, 139, 1);
		witness(new File(payload, "training.fa"), baseline, false, false);
		witness(new File(payload, "stale.fa"), baseline, true, false);
		witness(new File(payload, "reordered.fa"), baseline, false, true);
		final String good=recipe("sam1", true, 1, "training.fa", "five1.tsv", "three1.tsv");
		final String geometry="startinside=3\nstartoutside=2\nstopinside=5\nstopoutside=3\n";
		final String placed=recipe("sam1", true, 1, "training.fa", "placed_five.tsv", "placed_three.tsv")
			.replace("ncrna_positional_v1", "ncrna_positional_v2")+geometry;
		final File placedRoot=caseRoot(root, payload, "placed", placed, false);
		check(NcrnaFamilyLoader.load("sam1", options(core, placedRoot, null), new FamilyResourceResolver()).boundaryRefinerFactory!=null,
			"Explicit v2 packages bind distinct per-end K5/K8 anchors without changing the radius-sized neural interface");
		final File mismatch=caseRoot(root, payload, "placed_mismatch", placed.replace("stopoutside=3", "stopoutside=4"), false);
		fails(()->NcrnaFamilyLoader.load("sam1", options(core, mismatch, null), new FamilyResourceResolver()));
		final File incomplete=caseRoot(root, payload, "placed_incomplete", placed.replace("stopoutside=3\n", ""), false);
		fails(()->NcrnaFamilyLoader.load("sam1", options(core, incomplete, null), new FamilyResourceResolver()));
		final File legacyGeometry=caseRoot(root, payload, "legacy_geometry", good+geometry, false);
		fails(()->NcrnaFamilyLoader.load("sam1", options(core, legacyGeometry, null), new FamilyResourceResolver()));
		final File on=caseRoot(root, payload, "on", good, false);
		final NcrnaFamily loaded=NcrnaFamilyLoader.load("sam1", options(core, on, null), new FamilyResourceResolver());
		check(loaded.boundaryRefinerFactory!=null && loaded.boundary5NetTemplate==null && loaded.boundary3NetTemplate==null,
			"Enabled positional manifest installs only its custom factory through the actual loader");
		final File marginRoot=caseRoot(root, payload, "margin_defaults", good.replace("margin=0\n", "margin=0.25\n"), false);
		final NcrnaFamily inheritedMargin=NcrnaFamilyLoader.load("sam1", options(core, marginRoot, null), new FamilyResourceResolver());
		check(((NcrnaPositionalBoundaryRefiner.Resources)inheritedMargin.boundaryRefinerFactory).margin()==.25f,
			"Omitted finite margin fields must inherit the package's nonzero margin");
		final NcrnaFamilyLoader.Options zeroMargins=options(core, marginRoot, null);
		zeroMargins.boundaryMarginStart=zeroMargins.boundaryMarginStop=0f;
		zeroMargins.boundaryMarginStartSet=zeroMargins.boundaryMarginStopSet=true;
		final NcrnaFamily explicitZero=NcrnaFamilyLoader.load("sam1", zeroMargins, new FamilyResourceResolver());
		check(((NcrnaPositionalBoundaryRefiner.Resources)explicitZero.boundaryRefinerFactory).margin()==0f,
			"Explicit zero margins must override the package instead of being mistaken for omission");
		final NcrnaFamilyLoader.Options unequalMargins=options(core, marginRoot, null);
		unequalMargins.boundaryMarginStart=0f;
		unequalMargins.boundaryMarginStartSet=true;
		fails(()->NcrnaFamilyLoader.load("sam1", unequalMargins, new FamilyResourceResolver()));
		final NcrnaBoundaryRefiner first=loaded.boundaryRefinerFactory.create(loaded);
		check(first.afterAlignmentTrim() && first!=loaded.boundaryRefinerFactory.create(loaded),
			"An ungrouped registration binds its stable ID and supplies independent post-trim workers");
		final TrimSpy caller=new TrimSpy(loaded);
		caller.setBoundaryRefiner(first);
		final Orf candidate=new Orf("fixture", 75, 145, 0, 0, bases, false, ProkObject.RNA);
		candidate.ncrnaFamily=loaded.name;
		candidate.trnaModel=loaded.modelNames[0];
		candidate.orfScore=42;
		check(caller.finishAcceptedBoundary(candidate, bases, 0, 0, 239, 0, .8f, 71, .8f, true, null, null)
			&& caller.trims==1 && Math.abs(candidate.start-80)<=1 && Math.abs(candidate.stop-139)<=1
			&& candidate.orfScore==42, "Actual factory dispatch must use the80..139 trim result and preserve the caller score");
		final String missing=good.replace("net.bbnet", "missing.bbnet").replace("training.fa", "missing.fa");
		final File off=caseRoot(root, payload, "explicit_off", missing, false);
		final NcrnaFamily disabled=NcrnaFamilyLoader.load("sam1", options(core, off, false), new FamilyResourceResolver());
		check(disabled.boundaryRefinerFactory==null && disabled.boundary5NetTemplate==null,
			"Explicit OFF must not open missing witness/network payloads or install scalar fallback");
		fails(()->NcrnaFamilyLoader.load("sam1", options(core, off, null), new FamilyResourceResolver()));
		final File wrongLibrary=caseRoot(root, payload, "wrong_library", good.replace("training.fa", "stale.fa"), false);
		fails(()->NcrnaFamilyLoader.load("sam1", options(core, wrongLibrary, null), new FamilyResourceResolver()));
		final File wrongOrder=caseRoot(root, payload, "wrong_order", good.replace("training.fa", "reordered.fa"), false);
		fails(()->NcrnaFamilyLoader.load("sam1", options(core, wrongOrder, null), new FamilyResourceResolver()));
		final File wrongRadius=caseRoot(root, payload, "wrong_radius", good.replace("radius=1", "radius=2"), false);
		fails(()->NcrnaFamilyLoader.load("sam1", options(core, wrongRadius, null), new FamilyResourceResolver()));
		final File wrongWidth=caseRoot(root, payload, "wrong_width", recipe("sam1", true, 2, "training.fa", "five2.tsv", "three2.tsv"), false);
		fails(()->NcrnaFamilyLoader.load("sam1", options(core, wrongWidth, null), new FamilyResourceResolver()));
		final File wrongMember=caseRoot(root, payload, "wrong_member", good.replace("family=sam1", "family=cobalamin_rf00174"), false);
		fails(()->NcrnaFamilyLoader.load("sam1", options(core, wrongMember, null), new FamilyResourceResolver()));
		final File both=caseRoot(root, payload, "both_kinds", good, true);
		write(new File(both, "sliding.config"), "family=sam1\nenabled=f\nfeatures=v1\nfeatureorder="+PackagedBoundaryConfig.FEATURE_ORDER
			+"\nmeanlen="+PackagedBoundaryConfig.upperMedianLength(baseline.library)+"\nstartradius=1\nstopradius=1\nmarginstart=0\nmarginstop=0\ncutoff=none\n");
		fails(()->NcrnaFamilyLoader.load("sam1", options(core, both, null), new FamilyResourceResolver()));
		final File declaredOff=caseRoot(root, payload, "declared_off", recipe("sam1", false, 1, null, null, null), false);
		check(NcrnaFamilyLoader.load("sam1", options(core, declaredOff, null), new FamilyResourceResolver()).boundaryRefinerFactory==null,
			"A declared disabled positional package remains a legitimate ordinary-end outcome");
		fails(()->NcrnaFamilyLoader.load("sam1", options(core, declaredOff, true), new FamilyResourceResolver()));
		int bad=0;
		for(final String text:new String[]{good+"radius=1\n", good.replace("margin=0", "margin=NaN"),
			good+"oormargin=-0.1\n", good+"oormargin=ignore\n",
			good.replace("radius=1", "radius=0"), good.replace("ncrna_positional_v1", "unknown"),
			good.replace("net5=net.bbnet", "net5=../escape.bbnet"), good.replace("candidateIdentity", "wrongIdentity")}){
			final File path=new File(root, "bad"+(bad++)+".config");
			write(path, text);
			fails(()->PackagedPositionalBoundaryConfig.read(path.getPath()));
		}
		tripleCases(root, baseline, core, installation);
	}

	/** Gate metadata is complete and unambiguous before any payload is opened. */
	private static void gateConfigCases(final File root){
		final String base="schema=ncrna_positional_full_error_v1\nfeatures="+PackagedPositionalBoundaryConfig.TRIPLE_FEATURE_ORDER
			+"\nfamily=sam1\nenabled=t\nradius=10\nmargin=0\nstartinside=7\nstartoutside=2\nstopinside=8\nstopoutside=1\n"
			+"library=training.fa\ntraininghbm=training.hbm\nnet5=five.bbnet\nnet3=three.bbnet\n"
			+"table5=five0.tsv\ntable5left=five1.tsv\ntable5right=five2.tsv\ntable3=three0.tsv\ntable3left=three1.tsv\ntable3right=three2.tsv\n";
		final String[] fields={"gate5=gate5.bbnet\n", "gate3=gate3.bbnet\n", "gatecutoff5=0.1\n", "gatecutoff3=0.4522214\n",
			"gatefeatures="+NcrnaEndpointGateFeatures.FEATURE_ORDER+"\n"};
		final File deferred=new File(root, "deferred_cobalamin.config");
		write(deferred, base.replace("family=sam1", "family=cobalamin_rf00174")+"rawendpoints=t\ndeferreport=t\n");
		PackagedPositionalBoundaryConfig.read(deferred.getPath());
		final File ordinary=new File(root, "deferred_false.config");
		write(ordinary, base+"deferreport=f\n");
		PackagedPositionalBoundaryConfig.read(ordinary.getPath());
		final File wrongDeferred=new File(root, "deferred_wrong_family.config");
		write(wrongDeferred, base+"deferreport=t\n");
		fails(()->PackagedPositionalBoundaryConfig.read(wrongDeferred.getPath()));
		final ByteBuilder text=new ByteBuilder(base);
		for(final String field:fields){text.append(field);}
		final String good=text.toString();
		final File valid=new File(root, "gate_valid.config");
		write(valid, good);
		PackagedPositionalBoundaryConfig.read(valid.getPath());
		final String centered=good.replace("family=sam1", "family=cobalamin_rf00174")
			.replace(NcrnaEndpointGateFeatures.FEATURE_ORDER, NcrnaCenteredGateFeatures.FEATURE_ORDER)
			+"modeloffsets=offsets.tsv\nrawendpoints=t\ndeferreport=t\n";
		for(final String policy:new String[]{"0", "0.41595876", "never", "always"}){
			final File accepted=new File(root, "centered_gate_"+policy+".config");
			write(accepted, centered.replace("gatecutoff5=0.1", "gatecutoff5="+policy));
			PackagedPositionalBoundaryConfig.read(accepted.getPath());
		}
		int centeredBad=0;
		for(final String invalid:new String[]{centered.replace("modeloffsets=offsets.tsv\n", ""),
			centered.replace("deferreport=t", "deferreport=f"), centered.replace("rawendpoints=t", "rawendpoints=f"),
			centered.replace("family=cobalamin_rf00174", "family=cobalamin_rf01689"),
			centered.replace(NcrnaCenteredGateFeatures.FEATURE_ORDER, NcrnaEndpointGateFeatures.FEATURE_ORDER),
			centered.replace("gatecutoff5=0.1", "gatecutoff5=Infinity"),
			good.replace("gatecutoff5=0.1", "gatecutoff5=never"), good.replace("gatecutoff5=0.1", "gatecutoff5=always")}){
			final File rejected=new File(root, "centered_gate_bad"+(centeredBad++)+".config");
			write(rejected, invalid);
			fails(()->PackagedPositionalBoundaryConfig.read(rejected.getPath()));
		}
		for(int i=0; i<fields.length; i++){
			final File bad=new File(root, "gate_missing"+i+".config");
			write(bad, good.replace(fields[i], ""));
			fails(()->PackagedPositionalBoundaryConfig.read(bad.getPath()));
		}
		int count=0;
		for(final String invalid:new String[]{good.replace("gatecutoff5=0.1", "gatecutoff5=NaN"),
			good.replace("gatecutoff3=0.4522214", "gatecutoff3=Infinity"),
			good.replace("margin=0\n", "margin=0.1\n"), good+"oormargin=0\n",
			good.replace("rescaledRawIdentity", "wrongIdentity"),
			good.replace("ncrna_positional_full_error_v1", "ncrna_positional_triple_v1")}){
			final File bad=new File(root, "gate_invalid"+(count++)+".config");
			write(bad, invalid);
			fails(()->PackagedPositionalBoundaryConfig.read(bad.getPath()));
		}
	}

	/** A real 69-input export exercises the optional loader, with an HBM witness
	 * that can change coverage without changing any consensus sequence.
	 */
	private static void tripleCases(final File root, final NcrnaFamily baseline, final String[] core, final File installation){
		final File payload=new File(root, "triple_payload");
		check(payload.mkdir(), "Optional loader payload must be fresh");
		write(new File(payload, "train.tsv"), matrix(69, 21));
		write(new File(payload, "monitor.tsv"), matrix(69, 21));
		write(new File(payload, "recipe.config"), "dims=69,4,21\nseed=58\nepochs=2\nvfraction=0\nfinal=rslog\nhidden=tanh\nsimd=f\nt=1\nprintevery=1\now=t\n");
		new NcrnaNativeTrainer(installation, 1).trainPositional(path(payload, "recipe.config"),
			new String[]{path(payload, "train.tsv")}, path(payload, "monitor.tsv"), new File(payload, "net.bbnet"), 69, 21);
		final byte[] bases=new byte[240];
		Arrays.fill(bases, (byte)'A');
		final String[] tables={"five.tsv", "five_left.tsv", "five_right.tsv", "three.tsv", "three_left.tsv", "three_right.tsv"};
		final int[] offsets={-2, -9, 1, -7, -9, 1};
		for(int i=0; i<tables.length; i++){
			table(payload, tables[i], bases, "sam1", i<3 ? "5prime" : "3prime", offsets[i], i<3 ? 80 : 139, 10);
		}
		witness(new File(payload, "training.fa"), baseline, false, false);
		final String hbmName=new File(core[1]).getName();
		ReadWrite.copyFile(core[1], path(payload, hbmName));
		final BaseGraph[] changed=TrnaConsensusBuilder.loadModels(core[1]);
		final int[] counts=BaseGraphHelper.getCounts(changed[0], 0);
		BaseGraphHelper.setCounts(changed[0], 0, counts[0]+1, counts[1], counts[2], counts[3], counts[4]);
		TrnaConsensusBuilder.writeTextModels(new ArrayList<BaseGraph>(Arrays.asList(changed)), path(payload, "changed.hbm"));
		final ArrayList<BaseGraph> reordered=new ArrayList<BaseGraph>(Arrays.asList(baseline.models));
		check(reordered.size()>1, "HBM order fixture needs two models");
		final BaseGraph first=reordered.get(0);
		reordered.set(0, reordered.get(1));
		reordered.set(1, first);
		TrnaConsensusBuilder.writeTextModels(reordered, path(payload, "reordered.hbm"));
		reordered.remove(reordered.size()-1);
		TrnaConsensusBuilder.writeTextModels(reordered, path(payload, "short.hbm"));
		final String config="schema=ncrna_positional_triple_v1\nfeatures="+PackagedPositionalBoundaryConfig.TRIPLE_FEATURE_ORDER
			+"\nfamily=sam1\nenabled=t\nradius=10\nmargin=0\nstartinside=7\nstartoutside=2\nstopinside=8\nstopoutside=1\n"
			+"library=training.fa\ntraininghbm="+hbmName+"\nnet5=net.bbnet\nnet3=net.bbnet\n"
			+"table5=five.tsv\ntable3=three.tsv\ntable5left=five_left.tsv\ntable5right=five_right.tsv\n"
			+"table3left=three_left.tsv\ntable3right=three_right.tsv\n";
		final String[] files={"net.bbnet", "training.fa", hbmName, "changed.hbm", "reordered.hbm", "short.hbm",
			"five.tsv", "five_left.tsv", "five_right.tsv", "three.tsv", "three_left.tsv", "three_right.tsv"};
		final File good=tripleRoot(root, payload, "triple_on", config, files);
		check(NcrnaFamilyLoader.load("sam1", options(core, good, null), new FamilyResourceResolver()).boundaryRefinerFactory!=null,
			"Actual optional loader must accept six correctly ordered K9 tables, a real69-input export and matching HBM witness");
		for(final String bad:new String[]{"changed.hbm", "reordered.hbm", "short.hbm"}){
			final File invalid=tripleRoot(root, payload, "triple_"+bad, config.replace("traininghbm="+hbmName, "traininghbm="+bad), files);
			fails(()->NcrnaFamilyLoader.load("sam1", options(core, invalid, null), new FamilyResourceResolver()));
		}
		final File wrongOffset=tripleRoot(root, payload, "triple_wrong_offset", config.replace("table5left=five_left.tsv", "table5left=five_right.tsv"), files);
		fails(()->NcrnaFamilyLoader.load("sam1", options(core, wrongOffset, null), new FamilyResourceResolver()));
		String missing=config;
		for(final String file:files){missing=missing.replace("="+file+"\n", "=missing_"+file+"\n");}
		final File off=tripleRoot(root, payload, "triple_explicit_off", missing, new String[0]);
		check(NcrnaFamilyLoader.load("sam1", options(core, off, false), new FamilyResourceResolver()).boundaryRefinerFactory==null,
			"Explicit optional OFF cannot open any absent HBM witness, table, library or neural payload");
		fails(()->NcrnaFamilyLoader.load("sam1", options(core, off, null), new FamilyResourceResolver()));
	}

	static File tripleRoot(final File root, final File payload, final String name, final String config, final String[] files){
		final File out=new File(root, name);
		check(out.mkdir(), "Every optional loader fixture needs a fresh root");
		for(final String file:files){ReadWrite.copyFile(path(payload, file), path(out, file));}
		write(new File(out, "positional.config"), config);
		write(new File(out, "manifest.tsv"), "#ncrna_family\t1\tsam1\nrole\tmodel\tend\tpath\tformat\tbytes\tsha80\n"
			+"positional_config\t.\t.\tpositional.config\tconfig\t-\t-\n");
		return out;
	}

	static NcrnaFamilyLoader.Options options(final String[] core, final File root, final Boolean override){
		final NcrnaFamilyLoader.Options result=new NcrnaFamilyLoader.Options();
		result.resources.consensus=core[0];
		result.resources.models=core[1];
		result.resources.seeds=core[2];
		result.resources.root=root==null ? null : root.getPath();
		result.boundaryOverride=override;
		result.endpoints=false;
		return result;
	}
	private static File caseRoot(final File root, final File payload, final String name, final String recipe, final boolean sliding){
		final File result=new File(root, name);
		check(result.mkdir(), "Every loader case needs a fresh resource root");
		for(final String file:new String[]{"net.bbnet", "training.fa", "stale.fa", "reordered.fa", "five1.tsv", "three1.tsv", "five2.tsv", "three2.tsv", "placed_five.tsv", "placed_three.tsv"}){
			ReadWrite.copyFile(path(payload, file), path(result, file));
		}
		write(new File(result, "positional.config"), recipe);
		write(new File(result, "manifest.tsv"), "#ncrna_family\t1\tsam1\nrole\tmodel\tend\tpath\tformat\tbytes\tsha80\n"
			+"positional_config\t.\t.\tpositional.config\tconfig\t-\t-\n"
			+(sliding ? "boundary_config\t.\t.\tsliding.config\tconfig\t-\t-\n" : ""));
		return result;
	}
	private static void table(final File root, final String name, final byte[] bases, final String member,
			final String end, final int offset, final int truth, final int radius){
		table(root, name, bases, member, end, 9, offset, truth, radius);
	}
	private static void table(final File root, final String name, final byte[] bases, final String member,
			final String end, final int k, final int offset, final int truth, final int radius){
		final RrnaPositionalKmerTable.Table table=new RrnaPositionalKmerTable.Table(member, end, k, offset, 1.0/(1<<(2*k)), radius);
		table.add(bases, truth);
		table.finish();
		table.write(path(root, name));
	}
	/** Reordering preserves all records but breaks the model-index contract checked by the loader. */
	static void witness(final File out, final NcrnaFamily family, final boolean stale, final boolean reordered){
		check(!reordered || family.library.length>1, "The ordered-library rejection fixture needs two distinct model indices");
		try(Output output=new Output(out.getPath())){
			final ByteBuilder buffer=new ByteBuilder();
			for(int model=0; model<family.library.length; model++){
				final int index=reordered && model<2 ? 1-model : model;
				final byte[] bases=stale && model==0 ? family.library[index].clone() : family.library[index];
				if(stale && model==0){bases[0]=bases[0]=='A' ? (byte)'C' : (byte)'A';}
				output.writer.print(buffer.clear().append('>').append(family.modelNames[index]).nl().append(bases).nl());
			}
		}
	}
	private static String recipe(final String family, final boolean enabled, final int radius, final String library,
			final String five, final String three){
		return "schema=ncrna_positional_v1\nfeatures="+PackagedPositionalBoundaryConfig.FEATURE_ORDER+"\nfamily="+family
			+"\nenabled="+enabled+"\nradius="+radius+"\nmargin=0\n"+(enabled ? "library="+library+"\nnet5=net.bbnet\nnet3=net.bbnet\ntable5="+five+"\ntable3="+three+"\n" : "");
	}
	private static String matrix(){
		return matrix(6, 3);
	}
	private static String matrix(final int inputs, final int outputs){
		final ByteBuilder text=new ByteBuilder("#dims\t").append(inputs).tab().append(outputs).nl();
		for(int row=0; row<9; row++){
			for(int feature=0; feature<inputs; feature++){text.append((row+feature)%3).tab();}
			for(int label=0; label<outputs; label++){
				if(label>0){text.tab();}
				text.append(label==row%outputs ? 1 : 0);
			}
			text.nl();
		}
		return text.toString();
	}
	private static final class TrimSpy extends NcrnaScavenger{
		TrimSpy(final NcrnaFamily family){super(family.library, null, family.modelNames, family.kmerSet, 17, 16, 64);}
		@Override
		public boolean trimToAlignmentExtent(final Orf candidate, final byte[] bases, final int model,
				final int start, final int stop, final float identity){
			trims++;
			candidate.start=80;
			candidate.stop=139;
			return true;
		}
		int trims;
	}
	static String path(final File root, final String name){return new File(root, name).getPath();}
	static void write(final File path, final String text){try(Output out=new Output(path.getPath())){out.writer.print(text);}}
	private static void fails(final Runnable action){
		try{action.run();}catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Invalid positional package must be rejected before use");
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
}
