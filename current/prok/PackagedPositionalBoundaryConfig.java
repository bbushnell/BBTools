package prok;

import java.util.Arrays;
import java.util.Locale;

import consensus.BaseGraph;
import map.ObjectSet;
import ml.CellNetParser;
import parse.Parse;
import prok.NcrnaTextIO.ConfigInput;
import prok.NcrnaPositionalBoundaryRefiner.GatePolicy;

/** Explicit member-pooled positional resources for the ordinary family loader.
 * Configuration owns relative payload names, while FamilyResourceResolver owns
 * path construction. Disabled selection never opens tables or neural payloads.
 * An enabled factory copies its networks once per caller worker; tables remain
 * immutable and shared. This does not change verification, trim or path scores.
 * Full-error packages may set oormargin to a signed score gap or ignore; when
 * omitted, the extra head uses the effective ordinary margin, including CLI overrides.
 * Optional modeloffsets names the ordered six-field training export with support
 * floor ten; it centers both ends before features. Only the separately declared
 * 35-input whole-move gate can accompany offsets; old 32-input gates cannot.
 * Optional rawendpoints enables refinement even when trimming does not succeed.
 * Omission preserves the family's historical eligibility; physical/feature/length
 * guards still apply. This is an experimental policy change, not a default.
 * Optional deferreport retains ordinary post-trim bounds for cobalamin selection
 * and applies the effective refined report only after selection.
 * @author Brian Bushnell, Raiden
 */
public final class PackagedPositionalBoundaryConfig{

	/** Reads one complete configuration without opening any learned payload. */
	public static PackagedPositionalBoundaryConfig read(final String path){return new PackagedPositionalBoundaryConfig(path);}

	private PackagedPositionalBoundaryConfig(final String path){
		if(path==null){throw new IllegalArgumentException("A positional configuration needs its source path");}
		source=path;
		String schema=null;
		String featureOrder=null;
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		try(ConfigInput input=new ConfigInput(path)){
			while(input.next()){
				final String a=input.key.toLowerCase(Locale.ROOT), b=input.value;
				require(!a.isEmpty() && !b.isEmpty() && seen.add(a), "Configuration keys/values must be nonempty and unique: "+a);
				switch(a){
					case "schema":
						require(b.equals("ncrna_positional_v1") || b.equals("ncrna_positional_v2")
							|| b.equals("ncrna_positional_triple_v1") || b.equals("ncrna_positional_full_error_v1"),
							"Unsupported positional configuration schema");
						schema=b;
						break;
					case "features":
						featureOrder=b;
						break;
					case "family":
						family=NcrnaFamilyRegistry.sharedBoundaryMember(b).id;
						break;
					case "enabled":
						require(Parse.isBoolean(b) || b.equals("0") || b.equals("1"), "Enabled must be a boolean");
						enabled=Parse.parseBoolean(b);
						break;
					case "rawendpoints":
						require(Parse.isBoolean(b) || b.equals("0") || b.equals("1"), "Raw endpoints must be a boolean");
						rawEndpoints=Parse.parseBoolean(b);
						break;
					case "deferreport":
						require(Parse.isBoolean(b) || b.equals("0") || b.equals("1"), "Deferred reporting must be a boolean");
						deferReport=Parse.parseBoolean(b);
						break;
					case "radius":
						radius=NcrnaNumbers.intKMG(a, b);
						break;
					case "margin":
						margin=NcrnaNumbers.floatValue(a, b);
						break;
					case "oormargin":
						oorMarginSet=true;
						useOorHead=!b.equals("ignore");
						oorMargin=useOorHead ? NcrnaNumbers.floatValue(a, b) : 0;
						require(Float.isFinite(oorMargin), "Out-of-range margin must be finite and signed, or ignore");
						break;
					case "startinside":
						startInside=NcrnaNumbers.intValue(a, b);
						break;
					case "startoutside":
						startOutside=NcrnaNumbers.intValue(a, b);
						break;
					case "stopinside":
						stopInside=NcrnaNumbers.intValue(a, b);
						break;
					case "stopoutside":
						stopOutside=NcrnaNumbers.intValue(a, b);
						break;
					case "library":
						trainingLibrary=relative(b);
						break;
					case "net5":
						fiveNet=relative(b);
						break;
					case "net3":
						threeNet=relative(b);
						break;
					case "gate5": fiveGate=relative(b); break;
					case "gate3": threeGate=relative(b); break;
					case "gatecutoff5":
						fiveGatePolicy=gatePolicy(b);
						fiveGateCutoff=gateCutoff(a, b, fiveGatePolicy);
						break;
					case "gatecutoff3":
						threeGatePolicy=gatePolicy(b);
						threeGateCutoff=gateCutoff(a, b, threeGatePolicy);
						break;
					case "gatefeatures": gateFeatureOrder=b; break;
					case "table5":
						fiveTable=relative(b);
						break;
					case "table3":
						threeTable=relative(b);
						break;
					case "table5left": fiveLeft=relative(b); break;
					case "table5right": fiveRight=relative(b); break;
					case "table3left": threeLeft=relative(b); break;
					case "table3right": threeRight=relative(b); break;
					case "traininghbm": trainingModels=relative(b); break;
					case "modeloffsets": modelOffsets=relative(b); break;
					default: throw new IllegalArgumentException("Unknown positional configuration key "+a+" in "+source);
				}
			}
		}
		for(final String key:REQUIRED){require(seen.contains(key), "Missing positional configuration key: "+key);}
		fullErrors="ncrna_positional_full_error_v1".equals(schema);
		require(!deferReport || enabled && fullErrors && family!=null && family.startsWith("cobalamin_"),
			"Deferred reporting is an explicit enabled cobalamin full-error experiment");
		gated=false;
		for(final String key:GATE_FIELDS){gated|=seen.contains(key);}
		final boolean centeredGate=gated && NcrnaCenteredGateFeatures.FEATURE_ORDER.equals(gateFeatureOrder);
		require(modelOffsets==null || enabled && fullErrors && (!gated || centeredGate),
			"Model offsets require an enabled full-error package; gates trained on uncentered scores cannot be reused");
		for(final String key:GATE_FIELDS){
			require(seen.contains(key)==gated, "A movement gate declares both nets, both cutoffs and its exact feature order: "+key);
		}
		if(gated){
			require(enabled && fullErrors && margin==0 && (!oorMarginSet || !useOorHead)
				&& Float.isFinite(fiveGateCutoff) && Float.isFinite(threeGateCutoff),
				"Movement gates require enabled full-error zero-margin proposals, an ignored head and explicit per-end policies");
			if(centeredGate){
				require(modelOffsets!=null && deferReport && Boolean.TRUE.equals(rawEndpoints)
					&& (family.equals("cobalamin_rf00174") || family.equals("cobalamin_rf01482")),
					"The 35-input whole-move experiment requires bound offsets and explicit raw/deferred main or second cobalamin reporting");
			}else{
				require(NcrnaEndpointGateFeatures.FEATURE_ORDER.equals(gateFeatureOrder)
					&& fiveGatePolicy==GatePolicy.THRESHOLD && threeGatePolicy==GatePolicy.THRESHOLD,
					"Legacy movement gates require the exact32-input layout and finite cutoffs; never/always belong only to the centered35 route");
			}
		}
		require(!seen.contains("oormargin") || fullErrors, "Only the full-error schema has an out-of-range margin");
		triple=fullErrors || "ncrna_positional_triple_v1".equals(schema);
		require((triple ? TRIPLE_FEATURE_ORDER : FEATURE_ORDER).equals(featureOrder),
			"Feature order must match the selected positional training and inference schema");
		final boolean geometry=triple || "ncrna_positional_v2".equals(schema);
		for(final String key:GEOMETRY){
			require(geometry==seen.contains(key), "Positional v2 declares all four geometry fields; v1 retains its fixed K9 geometry: "+key);
		}
		NcrnaPositionalTables.wordLength(startInside, startOutside);
		NcrnaPositionalTables.wordLength(stopInside, stopOutside);
		if(triple){NcrnaPositionalTables.tripleWordLength(startInside, startOutside, stopInside, stopOutside);}
		for(final String key:TRIPLE_PAYLOAD){
			require(seen.contains(key)==(triple && enabled), "Only an enabled triple schema declares each extra table and its training HBM witness: "+key);
		}
		require(radius>0 && radius<=(Integer.MAX_VALUE-4)/2 && Float.isFinite(margin) && margin>=0,
			"Positional radius must fit feature widths and confidence margin must be finite/nonnegative");
		require(!fullErrors || radius==10, "Full-error v1 binds the reviewed 69-input/22-output radius-ten geometry");
		final boolean complete=fiveNet!=null && threeNet!=null && fiveTable!=null && threeTable!=null && trainingLibrary!=null;
		final boolean absent=fiveNet==null && threeNet==null && fiveTable==null && threeTable==null && trainingLibrary==null;
		require(enabled ? complete : absent, "Enabled positional configuration requires both nets/tables and the training-library witness; disabled package omits all five");
	}

	/** Explicit family/global ON/OFF overrides the package; scalar resources cannot mix.
	 * Explicit CLI margins, including zero, override the package as one common end margin.
	 * Radius overrides, if present, must describe this exact resource geometry.
	 */
	public void applyDefaults(final String name, final NcrnaFamilyLoader.Options options, final byte[][] library){
		require(options!=null && family.equals(name), "Positional configuration must bind the selected caller member");
		selected=options.boundaryOverride==null ? enabled : options.boundaryOverride;
		if(selected){
			require(fiveNet!=null && threeNet!=null && fiveTable!=null && threeTable!=null && trainingLibrary!=null,
				"Cannot enable a positional package with absent payload bindings");
			require(radius<=PackagedBoundaryConfig.upperMedianLength(library), "Positional radius exceeds the loaded library's upper-median length");
			require(options.resources.boundary==null && options.resources.boundary3==null
				&& options.resources.boundaryStart==null && options.resources.boundaryStop==null,
				"Positional configuration cannot also select scalar network/table overrides");
			final float five=options.boundaryMarginStartSet ? options.boundaryMarginStart : margin;
			final float three=options.boundaryMarginStopSet ? options.boundaryMarginStop : margin;
			require(Float.isFinite(five) && five>=0 && five==three,
				"This positional package requires one finite nonnegative confidence margin at both ends");
			margin=five;
			require(!gated || margin==0, "A movement gate replaces margin selection; a nonzero CLI margin cannot be combined with it");
			require((options.boundaryStartOffsets==null || NcrnaBoundaryScorer.centeredRadius(options.boundaryStartOffsets)==radius)
				&& (options.boundaryStopOffsets==null || NcrnaBoundaryScorer.centeredRadius(options.boundaryStopOffsets)==radius),
				"Explicit offsets cannot reinterpret the positional table/network radius");
		}
		// The custom post-trim factory replaces scalar refinement, including when OFF.
		options.boundary=false;
		configured=true;
	}

	/** Validates real table geometry and network widths before publishing a factory. */
	public void loaded(final NcrnaFamily loaded, final FamilyResourceResolver resolver){
		require(configured && loaded!=null && resolver!=null && loaded.name.equals(family),
			"Resolve the explicit ON/OFF selection before loading positional resources for their declared family");
		if(!selected){return;}
		if(rawEndpoints!=null){loaded.boundaryOnRawEndpoints=rawEndpoints;}
		loaded.deferBoundaryReport=deferReport;
		require(loaded.boundary5NetTemplate==null && loaded.boundary3NetTemplate==null
			&& loaded.boundaryRefinerFactory==null && loaded.rrnaEndpointFeatures==null,
			"Only one endpoint inference route may be installed on a caller family");
		final NcrnaFamilyRegistry.Entry entry=NcrnaFamilyRegistry.get(family);
		final String member=entry.group==null ? entry.id : entry.member;
		final NcrnaSupport.Library trained=NcrnaSupport.loadLibrary(
			resolver.configPayload(family, source, "positional_training_library", trainingLibrary).path);
		require(trained.names.length==loaded.library.length, "Positional network training library differs from the actual loaded library");
		for(int model=0; model<trained.names.length; model++){
			require(trained.names[model].equals(NcrnaSupport.recordName(loaded.modelNames[model]))
				&& Arrays.equals(trained.bases[model], loaded.library[model]),
				"Positional network training library differs in ordered model ID or sequence at index "+model);
		}
		if(fullErrors){
			loadFull(loaded, resolver, member);
			return;
		}
		final RrnaPositionalKmerTable.Table five=RrnaPositionalKmerTable.loadScores(
			resolver.configPayload(family, source, "positional_table5", fiveTable).path, member);
		final RrnaPositionalKmerTable.Table three=RrnaPositionalKmerTable.loadScores(
			resolver.configPayload(family, source, "positional_table3", threeTable).path, member);
		require(five.radius==radius && three.radius==radius && five.k==startInside+startOutside
			&& three.k==stopInside+stopOutside && five.offset==-startOutside && three.offset==1-stopInside,
			"Positional tables must match the package's declared radius, per-end word length and anchor placement");
		if(triple){
			final BaseGraph[] trainedModels=TrnaConsensusBuilder.loadModels(
				resolver.configPayload(family, source, "positional_training_hbm", trainingModels).path);
			verifyModels(trainedModels, loaded);
			final RrnaPositionalKmerTable.Table[] fiveProbes={five,
				RrnaPositionalKmerTable.loadScores(resolver.configPayload(family, source, "positional_table5left", fiveLeft).path, member),
				RrnaPositionalKmerTable.loadScores(resolver.configPayload(family, source, "positional_table5right", fiveRight).path, member)};
			final RrnaPositionalKmerTable.Table[] threeProbes={three,
				RrnaPositionalKmerTable.loadScores(resolver.configPayload(family, source, "positional_table3left", threeLeft).path, member),
				RrnaPositionalKmerTable.loadScores(resolver.configPayload(family, source, "positional_table3right", threeRight).path, member)};
			loaded.boundaryRefinerFactory=new NcrnaPositionalBoundaryRefiner.Resources(loaded, fiveProbes, threeProbes,
				CellNetParser.load(resolver.configPayload(family, source, "positional_net5", fiveNet).path),
				CellNetParser.load(resolver.configPayload(family, source, "positional_net3", threeNet).path), margin);
		}else{
			loaded.boundaryRefinerFactory=new NcrnaPositionalBoundaryRefiner.Resources(loaded, five, three,
				CellNetParser.load(resolver.configPayload(family, source, "positional_net5", fiveNet).path),
				CellNetParser.load(resolver.configPayload(family, source, "positional_net3", threeNet).path), margin);
		}
	}

	/** Signed tables are selected only by the explicit full-error package schema. */
	private void loadFull(final NcrnaFamily loaded, final FamilyResourceResolver resolver, final String member){
		assert(fullErrors && triple && radius==10) : "The full-error schema fixes both feature and output geometry before payload loading";
		final String[] fivePaths={fiveTable, fiveLeft, fiveRight}, threePaths={threeTable, threeLeft, threeRight};
		final NcrnaWeightedEndpointTable[] five=new NcrnaWeightedEndpointTable[3], three=new NcrnaWeightedEndpointTable[3];
		final RrnaPositionalKmerTable.Table[] fiveGeometry=new RrnaPositionalKmerTable.Table[3], threeGeometry=new RrnaPositionalKmerTable.Table[3];
		for(int probe=0; probe<3; probe++){
			five[probe]=NcrnaWeightedEndpointTable.loadScores(
				resolver.configPayload(family, source, "positional_signed_table5_"+probe, fivePaths[probe]).path, member);
			three[probe]=NcrnaWeightedEndpointTable.loadScores(
				resolver.configPayload(family, source, "positional_signed_table3_"+probe, threePaths[probe]).path, member);
			fiveGeometry[probe]=five[probe].geometry;
			threeGeometry[probe]=three[probe].geometry;
		}
		require(fiveGeometry[0].radius==radius && threeGeometry[0].radius==radius
			&& fiveGeometry[0].k==startInside+startOutside && threeGeometry[0].k==stopInside+stopOutside
			&& fiveGeometry[0].offset==-startOutside && threeGeometry[0].offset==1-stopInside,
			"Signed primary probes must match the declared package word length, placement and correction radius");
		verifyModels(TrnaConsensusBuilder.loadModels(
			resolver.configPayload(family, source, "positional_training_hbm", trainingModels).path), loaded);
		loaded.boundaryRefinerFactory=new NcrnaPositionalBoundaryRefiner.Resources(loaded, fiveGeometry, threeGeometry, five, three,
			CellNetParser.load(resolver.configPayload(family, source, "positional_net5", fiveNet).path),
			CellNetParser.load(resolver.configPayload(family, source, "positional_net3", threeNet).path), margin,
			gated || !useOorHead ? 0 : oorMarginSet ? oorMargin : margin, !gated && useOorHead,
			gated ? CellNetParser.load(resolver.configPayload(family, source, "positional_gate5", fiveGate).path) : null,
			gated ? CellNetParser.load(resolver.configPayload(family, source, "positional_gate3", threeGate).path) : null,
			fiveGateCutoff, threeGateCutoff, modelOffsets==null ? null : NcrnaEndpointModelOffsets.load(
				resolver.configPayload(family, source, "positional_model_offsets", modelOffsets).path, loaded, member, 10),
			fiveGatePolicy, threeGatePolicy);
	}

	/** HBM coverage is an input now: identical consensus bases alone are not a
	 * training/serve witness. Compare the actual ordered graphs and precisely the
	 * coordinate layout and coverage sums consumed by the shared fuzziness helper.
	 */
	private void verifyModels(final BaseGraph[] trained, final NcrnaFamily loaded){
		require(trained!=null && loaded.models!=null && trained.length==loaded.models.length,
			"Optional positional training HBM count differs from the actual loaded family");
		for(int model=0; model<trained.length; model++){
			final BaseGraph expected=trained[model], actual=loaded.models[model];
			require(expected!=null && actual!=null && expected.pad==actual.pad
				&& expected.ref.length==actual.ref.length && expected.del.length==actual.del.length
				&& NcrnaSupport.recordName(expected.name).equals(NcrnaSupport.recordName(actual.name))
				&& Arrays.equals(expected.original, actual.original),
				"Optional positional training HBM order, coordinates or original sequence differ at model "+model);
			for(int position=0; position<expected.ref.length; position++){
				require(expected.ref[position].countSum+expected.del[position].countSum
					==actual.ref[position].countSum+actual.del[position].countSum,
					"Optional positional training HBM coverage differs at model "+model+", position "+position);
			}
		}
	}

	private static String relative(final String path){
		FamilyResourceManifest.validateRelative(path);
		return path;
	}
	/** Explicit endpoints cannot be confused with arbitrary finite score thresholds. */
	private static GatePolicy gatePolicy(final String value){
		return value.equals("never") ? GatePolicy.NEVER : value.equals("always") ? GatePolicy.ALWAYS : GatePolicy.THRESHOLD;
	}
	private float gateCutoff(final String name, final String value, final GatePolicy policy){
		if(policy!=GatePolicy.THRESHOLD){return 0;}
		final float cutoff=NcrnaNumbers.floatValue(name, value);
		require(Float.isFinite(cutoff), "A numeric gate cutoff must be finite; use explicit never/always for centered controls");
		return cutoff;
	}
	private void require(final boolean condition, final String why){
		if(!condition){throw new IllegalArgumentException(why+" in "+source);}
	}

	private final String source;
	private String family, fiveNet, threeNet, fiveTable, threeTable, trainingLibrary;
	private String fiveLeft, fiveRight, threeLeft, threeRight, trainingModels;
	private String fiveGate, threeGate, gateFeatureOrder;
	private String modelOffsets;
	private boolean enabled, selected, configured, triple, fullErrors, gated;
	/** Explicit experimental eligibility override; omission preserves the family default. */
	private Boolean rawEndpoints;
	/** False preserves every family's ordinary acceptance and reporting geometry. */
	private boolean deferReport=false;
	/** Used only with the complete gate field set; a selected gate accepts score>=its cutoff. */
	private float fiveGateCutoff, threeGateCutoff;
	private GatePolicy fiveGatePolicy=GatePolicy.THRESHOLD, threeGatePolicy=GatePolicy.THRESHOLD;
	private int radius;
	/** Legacy v1 defaults; v2 must explicitly bind all four semantic counts. */
	private int startInside=9, startOutside=0, stopInside=9, stopOutside=0;
	private float margin;
	/** Omission follows the effective positional margin; ignoring the head is separate state. */
	private float oorMargin;
	private boolean oorMarginSet;
	private boolean useOorHead=true;
	public static final String FEATURE_ORDER="positionalScores,lengthRatio,contigGC,candidateIdentity";
	public static final String TRIPLE_FEATURE_ORDER="bestScores,leftScores,rightScores,lengthRatio,contigGC,candidateIdentity,hbmCoverageMinus,hbmCoverageAt,hbmCoveragePlus";
	private static final String[] REQUIRED={"schema", "features", "family", "enabled", "radius", "margin"};
	private static final String[] GEOMETRY={"startinside", "startoutside", "stopinside", "stopoutside"};
	private static final String[] TRIPLE_PAYLOAD={"table5left", "table5right", "table3left", "table3right", "traininghbm"};
	private static final String[] GATE_FIELDS={"gate5", "gate3", "gatecutoff5", "gatecutoff3", "gatefeatures"};
}
