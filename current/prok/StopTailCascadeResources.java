package prok;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;

import dna.Data;
import ddl.DDLFamilySearch.StructuralInputs;
import ddl.ProteinBlacklistPolicy;
import ddl.ProteinKmerGeometry;
import fileIO.ByteFile;
import parse.LineParser1;
import prot.DigestSuffix;
import prot.FamilyHeuristicRanking;
import prot.FamilyHeuristicScore;
import prot.FamilyGcInputs;
import prot.FamilyMemberScalarMeans;
import prot.FamilyTailCascade;
import prot.FamilyAlignmentMetric;
import prot.FamilyBlacklistBinding;
import prot.FamilyHitFloorBinding;
import prot.FamilyScoringPackage;
import prot.FamilyCompactCodec;
import prot.FamilyLengthPolicy;
import prot.ProteinSubstitutionMatrix;
import prot.FamilyCalibration.VerificationMode;
import shared.Tools;

/** Versioned resource binding for the callable stop-tail classifier.
 * Relative payload paths belong to the manifest directory; explicit absolute
 * paths support development resources. Version 2 explicitly binds the acceptance
 * margin scale and early-index rejection, and permits omission of both blacklist
 * bindings for models built without one. Version 3 additionally pins optional
 * 100–124 and 150–199 classifiers, each with its own explicit version-2 policy.
 * Version 4 instead requires an independent100–149 family package plus the
 * existing150–199 child. Version 1 preserves legacy defaults; version3 still
 * excludes125–149. Version5 adds an independent40–99 child to the version4 policy;
 * its optional short40minaa scanner minimum defaults to40 and may be raised to99.
 * Version5 defaults to the earliest configured-code start within the upstream
 * interval; explicit firststart=f preserves whole-interval queries.
 * Its40–99 version2 child may select admission=short_scaled or family_mean;
 * omitted admission or admission=fixed retains the fixed raw-hit floor.
 * An independent40–99 version2 child may explicitly select plain profiles with
 * positional=f and no partitions/overlap; positional children retain both fields.
 * Version5 defaults to searching the short set as well as the applicable long
 * view, preserving one binary contribution and separate per-set evidence.
 * Explicit searchall=f retains single-set length routing for comparisons.
 * Version6 explicitly binds an exact150 root,100..149 child and independent40..99
 * child for capped150 all-family-set search; it cannot inherit the200-aa root gate.
 * Version7 explicitly binds an exact100 root and independent40..99 child for
 * capped100 all-family-set search.
 * Version8 binds an exact120 root,100..119 child and independent40..99 child.
 * Earlier manifest versions retain their original query-boundary policy.
 * Version2+ may additionally bind child-local heuristic medians/KD and independent
 * retained-member scalar means. Omitted optional bindings preserve the old policy.
 * Version2+ can explicitly bind alignmentmode/alignmentminspan to the package
 * header; the clipped metric requires the released common pivot and old arrays.
 * Version2+ may bind aabits or aapartition for reduced-alphabet suffix residues.
 * Version2+ may bind positionmode explicitly, including fixed-quadrant. Omission
 * retains the classifier's legacy positional default; supplied modes must match
 * both the positional flag and stored profile layout when the classifier loads.
 * These options remain local to each classifier; children never inherit a
 * parent's geometry. The shared classifier validates the geometry and profiles.
 * Scientific validation stays in the classifier.
 * A #bbgap_tail_fitted version1 child instead pins its original fitting parent,
 * absolute margin, cohort manifest and saved floors. It is allowed only in the
 * independent40..99 slot and loads the same policy through FamilyHitFloorBinding.
 * @author Keqing, Fischl
 */
public final class StopTailCascadeResources{

	/** Read bindings without loading models or changing global nucleotide input mode. */
	public static StopTailCascadeResources read(final String name) throws IOException{
		final StopTailCascadeResources result=read(name, false, false, false, false);
		return requireRoot(result);
	}
	/** Explicit declared-identity binding; neither root, children nor payloads are hashed.
	 * Newly emitted resources use literal unmeasured; inherited identities retain their labels.
	 * Structural mode supports unfitted packages, including official and explicitly bound
	 * supported blacklist geometries, without scalar-mean replacement.
	 * Loading checks immutable filter keys.
	 * Schema, roster, policy and observable file stability remain checked.
	 * A compact package-root directory instead binds the approved SHORT and TAILS views.
	 */
	public static StopTailCascadeResources readStructural(final String name, final String inheritedPin) throws IOException{
		if(name!=null){
			final String resolved=Data.findPath(name,false);
			if(resolved!=null && Files.isDirectory(Paths.get(resolved))){
				if(!"unmeasured".equals(inheritedPin)){
					throw new IllegalArgumentException("Compact caller roots use the structural identity unmeasured");
				}
				return readCompact(Paths.get(resolved));
			}
		}
		return requireRoot(read(name,false,false,false,false,149,new StructuralInputs(),inheritedPin));
	}
	/** Bind the approved SHORT40..99 and TAILS100 views without inventing a legacy manifest.
	 * Both views must resolve the same alphabet, blacklist and matrix inside this root.
	 * Loading remains deferred until loadSupport; the accepted capped100/first-start
	 * routing and each header's absolute margin are retained. */
	public static StopTailCascadeResources readCompact(final Path packageRoot)throws IOException{
		if(packageRoot==null || !Files.isDirectory(packageRoot)){
			throw new IllegalArgumentException("Compact caller loading requires a package-root directory");
		}
		final Path root=packageRoot.toRealPath();
		try{
			final FamilyScoringPackage shortPack=new FamilyScoringPackage(
				compactFamily(root,"short"), VerificationMode.STRUCTURAL);
			final FamilyScoringPackage tailPack=new FamilyScoringPackage(
				compactFamily(root,"tails"), VerificationMode.STRUCTURAL);
			final Path[] shortInputs=shortPack.compactInputs(),tailInputs=tailPack.compactInputs();
			assert(shortInputs.length==4 && tailInputs.length==4)
				: "FamilyScoringPackage.compactInputs binds family, alphabet, blacklist and matrix in that order";
			for(int i=1;i<shortInputs.length;i++){
				if(!Files.isSameFile(shortInputs[i],tailInputs[i])){
					throw new IllegalArgumentException(
						"Compact caller views must share their alphabet, blacklist and matrix resources");
				}
			}
			final String provenance=compactValue(shortPack.compactHeader(),"packageprovenance");
			if(!provenance.equals(compactValue(tailPack.compactHeader(),"packageprovenance"))){
				throw new IllegalArgumentException("Compact caller views have different package provenance");
			}
			final StructuralInputs tracker=new StructuralInputs();
			final StopTailCascadeResources shortView=compactView(shortPack,tracker,null,40,99,.985);
			final StopTailCascadeResources result=compactView(tailPack,tracker,shortView,100,100,1);
			result.verifyManifest();return result;
		}catch(RuntimeException e){throw e;}catch(IOException e){throw e;}catch(Exception e){
			throw new IOException("Could not bind compact caller package: "+root,e);
		}
	}
	private static Path compactFamily(final Path root,final String view)throws IOException{
		assert(root!=null && (view.equals("short") || view.equals("tails")))
			: "The composite caller binds only the approved short and tails directories";
		final Path plain=root.resolve(view).resolve("families.tsv");
		final Path compressed=root.resolve(view).resolve("families.tsv.gz");
		if(Files.exists(plain)==Files.exists(compressed)){
			throw new IllegalArgumentException(
				"Compact caller view needs exactly one families.tsv or families.tsv.gz: "+view);
		}
		final Path selected=Files.exists(compressed) ? compressed : plain;
		if(Files.isSymbolicLink(selected) || !Files.isRegularFile(selected) || !selected.toRealPath().startsWith(root)){
			throw new IllegalArgumentException("Compact caller family file must be a regular in-package input: "+selected);
		}
		return selected;
	}
	private static StopTailCascadeResources compactView(final FamilyScoringPackage pack,final StructuralInputs tracker,
			final StopTailCascadeResources shortView,final int minimum,final int maximum,final double margin)throws Exception{
		final FamilyCompactCodec.Header header=pack.compactHeader();
		final HashMap<String,String> options=pack.resolvedOptions();
		if(Integer.parseInt(compactValue(header,"minqueryaa"))!=minimum ||
			Integer.parseInt(compactValue(header,"maxqueryaa"))!=maximum ||
			Double.compare(parseMarginScale(compactValue(header,"marginscale")),margin)!=0 ||
			!"hitfourth_h020_c150".equals(options.get("heuristic"))){
			throw new IllegalArgumentException(
				"Compact caller views require SHORT40..99/margin0.985 and TAILS100/margin1 with the released ranking");
		}
		// The released caller uses this measured-tail policy; omission must not
		// silently select the generic compact loader's optional no-GC path.
		final String gcFilter=compactValue(header,"gcfilter");
		if(Double.parseDouble(compactValue(header,"familygclambda"))!=24 ||
			!(gcFilter.equals("true") || gcFilter.equals("t")) ||
			Double.parseDouble(compactValue(header,"gc_z"))!=8 ||
			Integer.parseInt(compactValue(header,"gc_minimum_members"))!=30 ||
			Integer.parseInt(compactValue(header,"taillow"))!=minimum ||
			Integer.parseInt(compactValue(header,"tailhigh"))!=maximum){
			throw new IllegalArgumentException(
				"Compact caller GC policy requires coefficient24, filter=true, z8, "
				+ "minimum30 members and the matching tail interval");
		}
		final Path[] paths=pack.compactInputs();
		final ArrayList<Path> inputs=new ArrayList<Path>(paths.length);
		for(Path input:paths){tracker.track(input);inputs.add(input);}
		return new StopTailCascadeResources(paths[0],"unmeasured",options,inputs,margin,
			shortView,null,null,minimum,maximum,40,options.get("heuristic"),null,null,null,
			true,false,false,shortView!=null,null,tracker,null,null,pack);
	}
	private static String compactValue(final FamilyCompactCodec.Header header,final String key){
		final String value=header==null ? null : header.get(key);
		if(value==null || value.isEmpty()){throw new IllegalArgumentException("Compact caller header requires "+key);}
		return value;
	}
	private static StopTailCascadeResources requireRoot(final StopTailCascadeResources result){
		if(result.minQueryAa!=200 && !(result.minQueryAa==150 && result.maxQueryAa==150) &&
			!(result.minQueryAa==120 && result.maxQueryAa==120) &&
			!(result.minQueryAa==100 && result.maxQueryAa==100)){
			throw new IllegalArgumentException("A short-bin manifest must be bound as a child of a composite caller manifest");
		}
		return result;
	}
	/** Explicit research entry point; ordinary callers cannot load this fixed100-only manifest. */
	public static StopTailCascadeResources readFixed100(final String name) throws IOException{
		return read(name, false, true, false, false);
	}
	public static StopTailCascadeResources readFixed100Structural(final String name,
			final String inheritedPin)throws IOException{
		return read(name,false,true,false,false,149,new StructuralInputs(),inheritedPin);
	}
	private static StopTailCascadeResources read(final String name, final boolean child, final boolean fixed100,
			final boolean independent100, final boolean independent40) throws IOException{
		return read(name, child, fixed100, independent100, independent40, 149);
	}
	/** The composite schema supplies the lower child's ceiling; the child cannot widen it. */
	private static StopTailCascadeResources read(final String name, final boolean child, final boolean fixed100,
			final boolean independent100, final boolean independent40, final int lowerMaximum) throws IOException{
		return read(name,child,fixed100,independent100,independent40,lowerMaximum,null,null);
	}
	private static StopTailCascadeResources read(final String name, final boolean child, final boolean fixed100,
			final boolean independent100, final boolean independent40, final int lowerMaximum,
			final StructuralInputs structural,final String inheritedPin) throws IOException{
		assert(lowerMaximum==119 || lowerMaximum==149) :
			"Only cap120 and the retained cap150/200 lower-child contracts are released";
		if(name==null){throw new IllegalArgumentException("A stop-tail classifier manifest is required");}
		final String resolved=Data.findPath(name, false);
		if(resolved==null){throw new IllegalArgumentException("Missing stop-tail classifier manifest: "+name);}
		final Path manifest=Paths.get(resolved).toRealPath();
		if(structural!=null){
			if(!declaredIdentity(inheritedPin,true)){
				throw new IllegalArgumentException("Structural caller manifests require an explicit declared identity "
					+"(inherited sha80 or unmeasured)");
			}
			structural.track(manifest);
		}
		final String pin=structural==null ? DigestSuffix.file(manifest.toString()) : inheritedPin;
		final HashMap<String,String> options=new HashMap<String,String>();
		final ByteFile in=ByteFile.makeByteFile(manifest.toString(), false);
		final LineParser1 parser=new LineParser1('\t');
		final int version;
		try{
			final byte[] header=in.nextLine();
			final String schema=header==null ? "" : new String(header, StandardCharsets.UTF_8);
			if(FITTED_HEADER.equals(schema)){
				if(structural!=null){
					throw new IllegalArgumentException("Structural caller loading does not support fitted-floor children");
				}
				if(!child || !independent40 || fixed100){
					throw new IllegalArgumentException("Fitted floors require the independent40..99 child slot");
				}
				return readFitted(manifest,pin,in,parser);
			}
			version=HEADER.equals(schema) ? 1 :
				HEADER_V2.equals(schema) ? 2 :
				HEADER_V3.equals(schema) ? 3 :
				HEADER_V4.equals(schema) ? 4 :
				HEADER_V5.equals(schema) ? 5 :
				HEADER_V6.equals(schema) ? 6 :
				HEADER_V7.equals(schema) ? 7 :
				HEADER_V8.equals(schema) ? 8 : 0;
			if(version==0){
				throw new IllegalArgumentException("Expected a stop-tail manifest version1..8 in "+manifest);
			}
			if(fixed100 && version!=2){throw new IllegalArgumentException("Fixed100 requires a version-2 manifest");}
			if(child && version!=2){
				throw new IllegalArgumentException("A short-bin child must be version 2; "+
					"nested composite manifests are forbidden");
			}
			for(byte[] line;(line=in.nextLine())!=null;){
				if(line.length==0 || line[0]=='#'){continue;}
				parser.set(line);
				if(parser.terms()!=2){throw new IllegalArgumentException("Classifier manifest requires key and value: "+manifest);}
				final String key=parser.parseString(0).toLowerCase(Locale.ROOT), value=parser.parseString(1);
				if((!contains(KEYS, key) && !(version>=2 && (key.equals("marginscale") || key.equals("indexearlyreject")))
						&& !(version>=2 && contains(RANKING_KEYS,key))
						&& !(version>=2 && contains(GC_KEYS,key))
						&& !(version>=2 && contains(METRIC_KEYS,key))
						&& !(version>=2 && contains(DIAGNOSTIC_MATRIX_KEYS,key))
						&& !(version>=2 && contains(GEOMETRY_KEYS,key))
						&& !(version>=2 && contains(BLACKLIST_POLICY_KEYS,key))
						&& !(version>=3 && contains(SHORT_KEYS, key))
						&& !(version>=5 && (key.equals("short40manifest") || key.equals("short40sha80") ||
							key.equals("short40minaa") || key.equals("firststart") || key.equals("searchall")))
						&& !((version==2 || version==6 || version==7 || version==8) &&
							(key.equals("minqueryaa") || key.equals("maxqueryaa")))
						&& !(version==2 && key.equals("admission"))) ||
						value.isEmpty() || options.put(key, value)!=null){
					throw new IllegalArgumentException("Unknown, empty or duplicate classifier binding: "+key);
				}
			}
		}finally{if(in.close()){throw new IllegalStateException("Classifier manifest read failed: "+manifest);}}
		if(options.containsKey("alignmentminpaired") && !options.containsKey("alignmentmode")){
			throw new IllegalArgumentException("alignmentminpaired requires alignmentmode and alignmentminspan");
		}
		if(options.containsKey("alignmentmode")!=options.containsKey("alignmentminspan")){
			throw new IllegalArgumentException("alignmentmode and alignmentminspan must be supplied together");
		}
		// Historical manifests without a selector retain BLOSUM scores. A caller
		// consumes calibrated bindings; it never selects today's construction default.
		if(options.containsKey("matrix") && !"blosum62".equals(options.get("matrix")) &&
				!options.containsKey("matrixbinding")){
			throw new IllegalArgumentException("Custom caller matrices require their calibrated portable binding");
		}
		final prot.ProteinSubstitutionMatrix matrix=FamilyScoringPackage.resolveMatrix(options);
		if(!matrix.same(prot.ProteinSubstitutionMatrix.BLOSUM62) && !options.containsKey("alignmentmode")){
			throw new IllegalArgumentException("Custom caller matrices require a pinned alignment mode and span");
		}
		FamilyAlignmentMetric requestedMetric=FamilyAlignmentMetric.LEGACY;
		if(options.containsKey("alignmentmode")){
			FamilyAlignmentMetric.Mode mode=null;
			for(FamilyAlignmentMetric.Mode candidate:FamilyAlignmentMetric.Mode.values()){
				if(candidate.label.equals(options.get("alignmentmode"))){mode=candidate;}
			}
			if(mode==null){throw new IllegalArgumentException("Unknown caller alignmentmode");}
			requestedMetric=new FamilyAlignmentMetric(mode,Integer.parseInt(options.get("alignmentminspan")),matrix,
				options.containsKey("alignmentminpaired") ? Integer.parseInt(options.get("alignmentminpaired")) : 0);
		}
		final boolean noBlacklist=version>=2 && !options.containsKey("blacklist") && !options.containsKey("blacklistsha80");
		if(structural!=null && (options.containsKey("familymeans") || options.containsKey("familymeanssha80"))){
			throw new IllegalArgumentException("Structural caller loading does not support scalar-mean replacement");
		}
		if(version>=2 && options.containsKey("blacklist")!=options.containsKey("blacklistsha80")){
			throw new IllegalArgumentException("blacklist and blacklistsha80 must be supplied together or both omitted");
		}
		if(structural!=null && !noBlacklist){
			if(!FamilyBlacklistBinding.structuralFileGeometry(options)){
				throw new IllegalArgumentException("Structural caller filtering requires an official or explicit "+
					"supported blacklist geometry");
			}
			FamilyBlacklistBinding.validateOptions(options,VerificationMode.STRUCTURAL);
		}else{
			for(String key:BLACKLIST_POLICY_KEYS){
				if(options.containsKey(key)){
					FamilyBlacklistBinding.validateOptions(options,structural==null ? VerificationMode.SHA80 :
						VerificationMode.STRUCTURAL);
					break;
				}
			}
		}
		final String positionalText=options.get("positional");
		final boolean plainShort=version==2 && child && independent40 &&
			("f".equalsIgnoreCase(positionalText) || "false".equalsIgnoreCase(positionalText));
		// FamilyTailCascade.loadCallable forbids positional layout keys for plain
		// replacement profiles. Only the explicitly selected SHORT child gains this form.
		if(plainShort && (options.containsKey("partitions") || options.containsKey("overlap"))){
			throw new IllegalArgumentException("Plain SHORT caller child must omit partitions and overlap");
		}
		for(String key:KEYS){
			if(noBlacklist && (key.equals("blacklist") || key.equals("blacklistsha80"))){continue;}
			if(plainShort && (key.equals("partitions") || key.equals("overlap"))){continue;}
			if(!options.containsKey(key)){throw new IllegalArgumentException("Missing explicit classifier setting: "+key);}
		}
		// This caller policy is not a loadCallable option; apply its immutable view after loading.
		final double marginScale=version==1 ? 1 : parseMarginScale(options.remove("marginscale"));
		final String minText=options.remove("minqueryaa"), maxText=options.remove("maxqueryaa");
		if((minText==null)!=(maxText==null) || child && minText==null){
			throw new IllegalArgumentException("Short child manifests require explicit minqueryaa and maxqueryaa together");
		}
		final int minAa=minText==null ? 200 : Integer.parseInt(minText),
			maxAa=maxText==null ? 200 : Integer.parseInt(maxText);
		final String admission=options.remove("admission");
		if(admission!=null && !admission.equals("fixed") && !admission.equals("short_scaled") &&
			!admission.equals("family_mean")){
			throw new IllegalArgumentException("admission must be fixed, short_scaled or family_mean");
		}
		final boolean scaledAdmission="short_scaled".equals(admission);
		final boolean familyMeanAdmission="family_mean".equals(admission);
		if((scaledAdmission || familyMeanAdmission) && !(independent40 && child && minAa==40 && maxAa==99)){
			throw new IllegalArgumentException("Length-based admission is supported only by the independent40..99 child");
		}
		final String minimumText=options.remove("short40minaa");
		final String firstStartText=options.remove("firststart");
		if(firstStartText!=null && !firstStartText.equalsIgnoreCase("t") && !firstStartText.equalsIgnoreCase("true") &&
			!firstStartText.equalsIgnoreCase("f") && !firstStartText.equalsIgnoreCase("false")){
			throw new IllegalArgumentException("firststart must be an explicit boolean");
		}
		final boolean firstStart=firstStartText==null ? version>=5 :
			firstStartText.equalsIgnoreCase("t") || firstStartText.equalsIgnoreCase("true");
		final String allText=options.remove("searchall");
		if(allText!=null && !allText.equalsIgnoreCase("t") && !allText.equalsIgnoreCase("true") &&
			!allText.equalsIgnoreCase("f") && !allText.equalsIgnoreCase("false")){
			throw new IllegalArgumentException("searchall must be an explicit boolean");
		}
		final boolean allSets=allText==null ? version>=5 : allText.equalsIgnoreCase("t") || allText.equalsIgnoreCase("true");
		final int shortMinimum=minimumText==null ? 40 : Integer.parseInt(minimumText);
		if(shortMinimum<40 || shortMinimum>99){throw new IllegalArgumentException("short40minaa must be in [40,99]");}
		if(version==6 ? !(minText!=null && minAa==150 && maxAa==150) :
			version==7 ? !(minText!=null && minAa==100 && maxAa==100) :
			version==8 ? !(minText!=null && minAa==120 && maxAa==120) :
			fixed100 ? !(minText!=null && minAa==100 && maxAa==100) : independent40 ? !(child && minAa==40 && maxAa==99) :
			independent100 ? !(child && minAa==100 && maxAa==lowerMaximum) :
			!((minAa==100 && maxAa==124) || (minAa==150 && maxAa==199) || (!child && minAa==200 && maxAa==200))){
			throw new IllegalArgumentException("Unsupported caller query interval "+minAa+".."+maxAa+" for this manifest role");
		}
		if(version>=2){
			final String early=options.get("indexearlyreject");
			if(early==null){throw new IllegalArgumentException("Missing explicit classifier setting: indexearlyreject");}
			if(!early.equalsIgnoreCase("t") && !early.equalsIgnoreCase("true") &&
					!early.equalsIgnoreCase("f") && !early.equalsIgnoreCase("false")){
				throw new IllegalArgumentException("indexearlyreject must be an explicit boolean: "+early);
			}
			options.put("indexearlyreject", early.equalsIgnoreCase("t") || early.equalsIgnoreCase("true") ? "t" : "f");
		}
		final ArrayList<Path> inputs=new ArrayList<Path>();
		inputs.add(manifest);
		// This small diagnostic input always has a newly checked byte pin, even when
		// large inherited references use structural verification. Ordinary metric
		// options remain unchanged and must pass the original reference-header checks.
		final Payload evaluationMatrix=Payload.read(options,manifest.getParent(),inputs,"uncalibratedmatrixfile");
		if(evaluationMatrix!=null &&
			(scaledAdmission || familyMeanAdmission || !"pairwise40".equals(matrix.name))){
			throw new IllegalArgumentException("Uncalibrated matrix comparison requires fixed admission and "+
				"original pairwise40 calibration");
		}
		final String heuristicArm=options.remove("heuristicarm");
		final Payload familyMedians=Payload.read(options,manifest.getParent(),inputs,"familymedians",structural);
		final Payload familyHydro=Payload.read(options,manifest.getParent(),inputs,"familyhydro",structural);
		final Payload familyMeans=Payload.read(options,manifest.getParent(),inputs,"familymeans",structural);
		final GcPolicy gc=GcPolicy.read(options,manifest.getParent(),inputs,structural,minAa,maxAa);
		if(!requestedMetric.same(FamilyAlignmentMetric.LEGACY) &&
			(!"hitfourth_h020_c150".equals(heuristicArm) || familyMeans!=null)){
			throw new IllegalArgumentException("New-metric caller packages require the common pivot and "+
				"unchanged cosine arrays");
		}
		if((scaledAdmission || familyMeanAdmission) &&
			(!"hitfourth_h020_c150".equals(heuristicArm) || familyMeans!=null)){
			throw new IllegalArgumentException("The length-based admission comparison requires the "+
				"new1.5cosine/0.2KD pivot and unchanged family scalar arrays");
		}
		if(heuristicArm==null){
			if(familyMedians!=null || familyHydro!=null){
				throw new IllegalArgumentException("Heuristic medians/KD require an explicit heuristicarm");
			}
		}else{
			boolean supported=false;
			for(int arm=6;arm<FamilyHeuristicRanking.ALL_ARMS;arm++){
				if(arm!=7 && FamilyHeuristicRanking.armName(arm).equals(heuristicArm)){supported=true;break;}
			}
			if(!supported){throw new IllegalArgumentException("Unsupported caller heuristicarm: "+heuristicArm);}
			if(familyMedians==null || (!heuristicArm.equals("hitfourth"))!=(familyHydro!=null)){
				throw new IllegalArgumentException("Heuristic arms require pinned medians and H arms require "+
					"pinned member KD; H0 must omit KD");
			}
			if(minAa<40 || maxAa>200 || minAa>maxAa){
				throw new IllegalArgumentException("Caller heuristic bindings require one validated query interval within40..200");
			}
		}
		for(String key:PATH_KEYS){
			if(noBlacklist && key.equals("blacklist")){continue;}
			final Path path=resolve(manifest.getParent(), options.get(key));
			if(Files.isSymbolicLink(path)){
				throw new IllegalArgumentException("Callable classifier requires a regular non-symlink input: "+path);
			}
			final String hash=options.get(key+"sha80");
			if(!declaredIdentity(hash,structural!=null)){
				throw new IllegalArgumentException(structural==null ? "Expected sha80 for "+key :
					"Invalid declared structural identity for "+key);
			}
			options.put(key, path.toString());
			inputs.add(path);
			if(structural!=null){structural.track(path);}
		}
		final Path profiles=Paths.get(options.get("ddlprofiles"));
		if(structural==null && !DigestSuffix.file(profiles.toString()).equals(options.get("ddlprofilessha80"))){
			throw new IllegalArgumentException("Profile list digest differs from classifier manifest: "+profiles);
		}
		addProfileInputs(profiles,options.get("ddlprofilessha80"),inputs,structural);
		final int lowerCap=version==8 ? 119 : 149;
		final StopTailCascadeResources short40=version>=5 ?
			child(manifest.getParent(), options, "short40", false, lowerCap,structural) : null;
		final StopTailCascadeResources short100=child(manifest.getParent(), options, "short100",
			version>=4, lowerCap,structural);
		final StopTailCascadeResources short150=child(manifest.getParent(), options, "short150", false, lowerCap,structural);
		if(version==3 && short100==null && short150==null){
			throw new IllegalArgumentException("Version 3 requires at least one explicitly bound short bin");
		}
		if(version==4 && (short100==null || short150==null)){
			throw new IllegalArgumentException("Version4 requires both100–149 and150–199 children");
		}
		if(version==5 && (short40==null || short100==null || short150==null)){
			throw new IllegalArgumentException("Version5 requires40–99,100–149 and150–199 children");
		}
		if(version==6 && (short40==null || short100==null || short150!=null || !allSets)){
			throw new IllegalArgumentException("Version6 requires40–99 and100–149 children, no150–199 child and searchall=t");
		}
		if(version==8 && (short40==null || short100==null || short150!=null || !allSets)){
			throw new IllegalArgumentException("Version8 requires40–99 and100–119 children, no150–199 child and searchall=t");
		}
		if(version==7 && (short40==null || short100!=null || short150!=null || !allSets)){
			throw new IllegalArgumentException("Version7 requires a40–99 child, no100–149/150–199 children and searchall=t");
		}
		if(short40!=null && (short40.lengthScaledAdmission || short40.familyMeanAdmission) && shortMinimum!=40){
			throw new IllegalArgumentException("The length-based admission caller comparison requires minimum40");
		}
		if(allSets){
			if(shortMinimum!=40 || short40==null || version!=7 && short100==null ||
				version!=6 && version!=7 && version!=8 && short150==null ||
				short40.lengthScaledAdmission || short40.familyMeanAdmission){
				throw new IllegalArgumentException("searchall requires explicitly calibrated bins, minimum40 and "+
					"fixed short admission");
			}
			if(!"hitfourth_h020_c150".equals(heuristicArm) || !"hitfourth_h020_c150".equals(short40.heuristicArm) ||
				short100!=null && !"hitfourth_h020_c150".equals(short100.heuristicArm) ||
				short150!=null && !"hitfourth_h020_c150".equals(short150.heuristicArm) ||
				familyMeans!=null || short40.familyMeans!=null ||
				short100!=null && short100.familyMeans!=null || short150!=null && short150.familyMeans!=null){
				throw new IllegalArgumentException("searchall requires the standardized1.5cosine/0.2KD "+
					"pivot with unchanged cosine arrays in every view");
			}
			if(short100!=null && !sameProfiles(options,short100,structural) ||
				short150!=null && !sameProfiles(options,short150,structural) || sameProfiles(options,short40,structural)){
				throw new IllegalArgumentException("searchall requires a distinct short-family index and "+
					"one shared index across any long views");
			}
		}
		final StopTailCascadeResources result=new StopTailCascadeResources(manifest, pin, options, inputs, marginScale,
			short40, short100, short150, minAa, maxAa, shortMinimum,heuristicArm,familyMedians,
			familyHydro,familyMeans,firstStart,scaledAdmission,familyMeanAdmission,allSets,null,
			structural,evaluationMatrix,gc,null);
		result.validateCompositeEvaluationMatrix();
		result.verifyManifest();
		return result;
	}
	/** Every explicit view participates in the same diagnostic; a root-only swap is ambiguous. */
	private void validateCompositeEvaluationMatrix()throws IOException{
		for(StopTailCascadeResources child:new StopTailCascadeResources[]{short40,short100,short150}){
			if(child==null){continue;}
			if((evaluationMatrix==null)!=(child.evaluationMatrix==null)){
				throw new IllegalArgumentException("Composite matrix diagnostic requires the root and every "+
					"child to bind the same evaluation matrix");
			}
			if(evaluationMatrix!=null && (!evaluationMatrix.pin.equals(child.evaluationMatrix.pin) ||
				!Files.isSameFile(evaluationMatrix.path,child.evaluationMatrix.path))){
				throw new IllegalArgumentException("Composite matrix diagnostic has mixed evaluation-matrix files or pins");
			}
		}
	}

	private StopTailCascadeResources(final Path manifest_, final String pin_, final HashMap<String,String> options_,
			final ArrayList<Path> inputs_, final double marginScale_,
			final StopTailCascadeResources short40_, final StopTailCascadeResources short100_,
			final StopTailCascadeResources short150_,
			final int minAa_, final int maxAa_, final int shortMinimum_,final String heuristicArm_,
			final Payload familyMedians_,final Payload familyHydro_,final Payload familyMeans_,
			final boolean firstStart_,final boolean scaledAdmission_,
			final boolean familyMeanAdmission_,final boolean allSets_,final FittedPolicy fitted_,
			final StructuralInputs structural_,
			final Payload evaluationMatrix_,final GcPolicy gc_,final FamilyScoringPackage compactPackage_){
		assert(Double.isFinite(marginScale_) && marginScale_>=0) :
			"Resource load must preserve a validated immutable margin policy";
		manifest=manifest_; sha80=pin_; options=options_; inputs=inputs_; marginScale=marginScale_;
		short40=short40_; short100=short100_; short150=short150_;
		minQueryAa=minAa_; maxQueryAa=maxAa_;
		shortMinimum=shortMinimum_;
		firstStart=firstStart_;
		lengthScaledAdmission=scaledAdmission_;
		familyMeanAdmission=familyMeanAdmission_;
		allFamilySets=allSets_;
		fitted=fitted_;
		structuralInputs=structural_;
		evaluationMatrix=evaluationMatrix_;
		gc=gc_;
		compactPackage=compactPackage_;
		heuristicArm=heuristicArm_;familyMedians=familyMedians_;familyHydro=familyHydro_;familyMeans=familyMeans_;
	}

	/** A fitted child names one immutable parent policy instead of restating its
	 * options. This prevents a caller-side override from changing a fitted table's
	 * gates or ranking, and delays the sole index load until loadSupport. */
	private static StopTailCascadeResources readFitted(final Path manifest,final String pin,
			final ByteFile input,final LineParser1 parser)throws IOException{
		final HashMap<String,String> bindings=new HashMap<String,String>();
		for(byte[] line;(line=input.nextLine())!=null;){
			if(line.length==0 || line[0]=='#'){continue;}
			parser.set(line);
			if(parser.terms()!=2){throw new IllegalArgumentException("Fitted child requires key and value: "+manifest);}
			final String key=parser.parseString(0).toLowerCase(Locale.ROOT),value=parser.parseString(1);
			if(!contains(FITTED_KEYS,key) || value.isEmpty() || bindings.put(key,value)!=null){
				throw new IllegalArgumentException("Unknown, empty or duplicate fitted-child binding: "+key);
			}
		}
		for(String key:FITTED_KEYS){
			if(!bindings.containsKey(key)){throw new IllegalArgumentException("Missing fitted-child binding: "+key);}
		}
		final ArrayList<Path> inputs=new ArrayList<Path>();inputs.add(manifest);
		final Payload parent=Payload.read(bindings,manifest.getParent(),inputs,"parent");
		final Payload cohorts=Payload.read(bindings,manifest.getParent(),inputs,"cohorts");
		final Payload floors=Payload.read(bindings,manifest.getParent(),inputs,"floors");
		final double margin=parseMarginScale(bindings.remove("marginscale"));
		final HashMap<String,String> options;
		try{
			final FamilyScoringPackage pack=new FamilyScoringPackage(parent.path);
			options=pack.resolvedOptions();pack.verify();
		}catch(RuntimeException e){throw e;}catch(Exception e){
			throw new IOException("Could not resolve the pinned fitted-child parent: "+parent.path,e);
		}
		if(!parent.pin.equals(DigestSuffix.file(parent.path.toString())) ||
			!"hitfourth_h020_c150".equals(options.get("heuristic")) ||
			!"40".equals(options.get("minqueryaa")) || !"99".equals(options.get("maxqueryaa")) ||
			!options.containsKey("ddlprofiles") || !options.containsKey("ddlprofilessha80")){
			throw new IllegalArgumentException("Fitted child requires its unchanged common-pivot40..99 parent");
		}
		// Keep all underlying files in the same output-protection roster as ordinary children.
		for(String key:PATH_KEYS){
			if(!options.containsKey(key)){continue;}
			final Path path=resolve(parent.path.getParent(),options.get(key));
			if(Files.isSymbolicLink(path)){throw new IllegalArgumentException("Fitted child input cannot be a symlink: "+path);}
			options.put(key,path.toString());inputs.add(path);
		}
		options.put("familymedians",options.remove("lengthmedians"));
		options.put("familymedianssha80",options.remove("lengthmedianssha80"));
		final Payload medians=Payload.read(options,parent.path.getParent(),inputs,"familymedians");
		final Payload hydro=Payload.read(options,parent.path.getParent(),inputs,"familyhydro");
		addProfileInputs(Paths.get(options.get("ddlprofiles")),options.get("ddlprofilessha80"),inputs);
		final StopTailCascadeResources result=new StopTailCascadeResources(manifest,pin,options,inputs,margin,
			null,null,null,40,99,40,"hitfourth_h020_c150",medians,hydro,null,false,false,false,false,
			new FittedPolicy(parent,cohorts,floors),null,null,null,null);
		result.verifyManifest();return result;
	}

	/** Output protection includes profile members and sibling summaries, not just their list. */
	private static void addProfileInputs(final Path profiles,final String pin,
			final ArrayList<Path> inputs)throws IOException{
		addProfileInputs(profiles,pin,inputs,null);
	}
	private static void addProfileInputs(final Path profiles,final String pin,final ArrayList<Path> inputs,
			final StructuralInputs structural)throws IOException{
		assert(profiles!=null && pin!=null && inputs!=null) :
			"Resource protection needs the complete pinned profile-list roster";
		if(structural!=null){structural.track(profiles);}
		final ByteFile list=ByteFile.makeByteFile(profiles.toString(),false);
		final LineParser1 parser=new LineParser1('\t');int entries=0;
		try{
			for(byte[] line;(line=list.nextLine())!=null;){
				parser.set(line);
				if(parser.terms()!=3){
					throw new IllegalArgumentException("Profile list requires path and two digests: "+profiles);
				}
				final Path member=resolve(profiles.getParent(),parser.parseString(0));
				final Path summary=structural==null ? member.getParent().resolve("summary.tsv").toRealPath() :
					member.getParent().resolve("summary.tsv").toAbsolutePath().normalize();
				inputs.add(member);inputs.add(summary);entries++;
				if(structural!=null){structural.track(member);structural.track(summary);}
			}
		}finally{if(list.close()){throw new IllegalStateException("Profile list read failed: "+profiles);}}
		if(entries==0 || structural==null && !pin.equals(DigestSuffix.file(profiles.toString()))){
			throw new IllegalStateException("Profile list is empty or changed during resource binding: "+profiles);
		}
		if(structural!=null){structural.verify();}
	}

	/** Only this complete binding can select the fitted loader branch. */
	private static final class FittedPolicy{
		FittedPolicy(Payload parent_,Payload cohorts_,Payload floors_){
			assert(parent_!=null && cohorts_!=null && floors_!=null) :
				"Fitted loading requires all three pinned artifacts";
			parent=parent_;
			cohorts=cohorts_;
			floors=floors_;
		}
		final Payload parent,cohorts,floors;
	}

	/** Each child is independently pinned and explicitly scaled; a parent's margin is never inherited. */
	private static StopTailCascadeResources child(final Path parent, final HashMap<String,String> options,
			final String prefix, final boolean independent100, final int lowerMaximum,
			final StructuralInputs structural) throws IOException{
		final String name=options.remove(prefix+"manifest"), expected=options.remove(prefix+"sha80");
		if(name==null && expected==null){return null;}
		if(name==null || !declaredIdentity(expected,structural!=null)){
			throw new IllegalArgumentException(prefix+"manifest and "+prefix+"sha80 must be supplied "+
				"together with a sha80 digest");
		}
		final Path path=resolve(parent, name);
		if(Files.isSymbolicLink(path)){
			throw new IllegalArgumentException("Short classifier manifest cannot be a symlink: "+path);
		}
		final StopTailCascadeResources result=read(path.toString(), true, false, independent100,
			prefix.equals("short40"), lowerMaximum,structural,expected);
		if(!result.sha80.equals(expected)){
			throw new IllegalArgumentException("Short classifier manifest digest mismatch: "+path);
		}
		final int minimum=prefix.equals("short40") ? 40 : prefix.equals("short100") ? 100 : 150;
		if(result.minQueryAa!=minimum){
			throw new IllegalArgumentException(prefix+" binding names a different query-length bin: "+path);
		}
		return result;
	}

	/** Load the root view only (200-aa historically, exact100/120/150 for capped roots). Use loadSupport for routing. */
	public FamilyTailCascade.LoadedClassifier load(final int threads) throws Exception{
		return loadWithOverrides(threads,0,-1,0);
	}
	/** Explicit research overrides apply to a detached policy map, never the frozen manifest.
	 * Zero cap, floor -1 and zero reduction retain their respective original settings.
	 */
	public FamilyTailCascade.LoadedClassifier loadWithOverrides(final int threads,final int top,final int minHits,
			final float identityLowerReduction)throws Exception{
		validateOverrides(threads,top,minHits,identityLowerReduction);
		if((top>0 || minHits>=0 || identityLowerReduction!=0) &&
			(fitted!=null || lengthScaledAdmission || familyMeanAdmission)){
			throw new IllegalArgumentException("Caller overrides require an unfitted fixed-admission classifier");
		}
		verifyManifest();
		if(compactPackage!=null){
			FamilyTailCascade.LoadedClassifier loaded=compactPackage.loadWithOverrides(threads,top,minHits)
				.withEvaluationMarginScale(marginScale);
			if(identityLowerReduction!=0){loaded=loaded.withIdentityLowerReduction(identityLowerReduction);}
			verifyManifest();return loaded;
		}
		if(fitted!=null){
			final FamilyHitFloorBinding binding=new FamilyHitFloorBinding(fitted.parent.path,fitted.parent.pin,
				marginScale,fitted.cohorts.path,fitted.cohorts.pin,threads);
			final FamilyTailCascade.LoadedClassifier loaded=binding.load(fitted.floors.path,fitted.floors.pin);
			verifyManifest();return loaded;
		}
		final HashMap<String,String> selected=new HashMap<String,String>(options);
		if(top>0){selected.put("top",Integer.toString(top));}
		if(minHits>=0){selected.put("ddlreject",Integer.toString(minHits));}
		FamilyTailCascade.LoadedClassifier loaded=(structuralInputs==null ? FamilyTailCascade.loadCallable(selected,threads) :
			FamilyTailCascade.loadCallableStructural(selected,threads,structuralInputs)).withEvaluationMarginScale(marginScale);
		if(lengthScaledAdmission){loaded=loaded.withLengthScaledShortAdmission();}
		if(familyMeanAdmission){loaded=loaded.withFamilyMeanAdmission();}
		if(familyMeans!=null){
			loaded=loaded.withMemberScalarMeans(new FamilyMemberScalarMeans(loaded,familyMeans.path,
				familyMeans.pin));
		}
		if(heuristicArm!=null){
			loaded=loaded.withHeuristicRanking(new FamilyHeuristicScore(loaded,familyMedians.path,familyMedians.pin,
				familyHydro==null ? null : familyHydro.path,familyHydro==null ? null : familyHydro.pin,
				minQueryAa,maxQueryAa,heuristicArm,
				structuralInputs==null ? VerificationMode.SHA80 : VerificationMode.STRUCTURAL));
		}
		if(identityLowerReduction!=0){loaded=loaded.withIdentityLowerReduction(identityLowerReduction);}
		if(gc!=null){
			loaded=loaded.withFamilyGcPenalty(FamilyGcInputs.loadTailStats(gc.stats.path,gc.map.path,gc.families,gc.members,
				minQueryAa,maxQueryAa,gc.lambda,gc.filter,gc.z,structuralInputs));
		}
		if(evaluationMatrix!=null){
			loaded=loaded.withUncalibratedMatrix(ProteinSubstitutionMatrix.read(evaluationMatrix.path.toString()));
		}
		verifyManifest();
		return loaded;
	}

	/** Load all explicitly bound bins once, retaining separate reference gates and margin scales. */
	public StopTailCascadeSupport loadSupport(final int threads) throws Exception{
		return loadSupportWithOverrides(threads,0,-1,0);
	}
	/** Apply the same explicitly requested override to every loaded LONG/SHORT view. */
	public StopTailCascadeSupport loadSupportWithOverrides(final int threads,final int top,final int minHits,
			final float identityLowerReduction)throws Exception{
		return loadSupportWithOverrides(threads,top,minHits,identityLowerReduction,null);
	}
	public StopTailCascadeSupport loadSupportWithOverrides(final int threads,final int top,final int minHits,
			final float identityLowerReduction,final prot.FamilyMatchInformation information)throws Exception{
		return loadSupportWithOverrides(threads,top,minHits,identityLowerReduction,information,null);
	}
	/** Apply one shared coefficient after loading each view's frozen gates and ranking.
	 * The caller owns policy-file stability and output protection; this boundary
	 * checks its parsed immutable policy against the actual joint resource binding.
	 * A null policy preserves the original loading path and native span minima.
	 */
	public StopTailCascadeSupport loadSupportWithOverrides(final int threads,final int top,final int minHits,
			final float identityLowerReduction,final prot.FamilyMatchInformation information,
			final FamilyLengthPolicy lengthPolicy)throws Exception{
		validateOverrides(threads,top,minHits,identityLowerReduction);
		if(lengthPolicy!=null){validateLengthPolicy(lengthPolicy,top,minHits,identityLowerReduction);}
		final boolean capped100=minQueryAa==100 && maxQueryAa==100 && allFamilySets;
		final boolean capped120=minQueryAa==120 && maxQueryAa==120;
		final boolean capped150=minQueryAa==150 && maxQueryAa==150;
		if(minQueryAa!=200 && !capped100 && !capped120 && !capped150){
			throw new IllegalArgumentException("Caller support requires a full-length or explicit capped root manifest");
		}
		verifyManifest();
		FamilyTailCascade.LoadedClassifier full=loadWithOverrides(threads,top,minHits,identityLowerReduction);
		FamilyTailCascade.LoadedClassifier tiny=short40==null ? null :
			short40.loadWithOverrides(threads,top,minHits,identityLowerReduction);
		final FamilyTailCascade.LoadedClassifier small=short100==null ? null :
			short100.loadWithOverrides(threads,top,minHits,identityLowerReduction);
		final FamilyTailCascade.LoadedClassifier medium=short150==null ? null :
			short150.loadWithOverrides(threads,top,minHits,identityLowerReduction);
		if(lengthPolicy!=null){
			full=full.withLengthPolicy(lengthPolicy);
			tiny=tiny.withLengthPolicy(lengthPolicy);
		}
		verifyManifest();
		if(information!=null){
			return new StopTailCascadeSupport(full,small,medium,threads,false,
				tiny!=null || short100!=null && short100.maxQueryAa==149,tiny,shortMinimum,firstStart,
				allFamilySets,maxQueryAa,information);
		}
		if(capped100){return StopTailCascadeSupport.capped100(full,tiny,threads,firstStart);}
		if(capped120){return StopTailCascadeSupport.capped120(full,tiny,small,threads,firstStart);}
		if(capped150){return StopTailCascadeSupport.capped150(full,tiny,small,threads,firstStart);}
		if(allFamilySets){return StopTailCascadeSupport.allFamilySets(full,tiny,small,medium,threads,firstStart);}
		if(tiny!=null){
			return StopTailCascadeSupport.independent40To99(full, tiny, small, medium, threads,
				shortMinimum, firstStart);
		}
		return short100!=null && short100.maxQueryAa==149 ?
			StopTailCascadeSupport.independent100To149(full, small, medium, threads)
			: new StopTailCascadeSupport(full, small, medium, threads);
	}
	/** Preflight before opening outputs; loading repeats this check before workers exist. */
	public void validateLengthPolicy(final FamilyLengthPolicy policy,final int top,final int minHits,
			final float reduction){
		if(policy==null){throw new IllegalArgumentException("Length-policy preflight requires a parsed policy");}
		if(evaluationMatrix!=null ||
			short40!=null && short40.evaluationMatrix!=null){
			throw new IllegalArgumentException("A calibrated length policy cannot consume an uncalibrated "+
				"evaluation-matrix swap");
		}
		policy.requireReference(manifest.toString(),sha80);
		if(minQueryAa!=100 || maxQueryAa!=100 || !allFamilySets || shortMinimum!=40 ||
			short40==null || short40.minQueryAa!=40 || short40.maxQueryAa!=99 || short100!=null || short150!=null){
			throw new IllegalArgumentException("Shared length policy requires the exact LONG100 and SHORT40..99 "+
				"all-set resource");
		}
		if(top!=0 && top!=50 || minHits!=-1 && minHits!=5 || reduction!=0){
			throw new IllegalArgumentException("Shared length policy preserves top50, raw-hit floor5 and the "+
				"original identity lower gate");
		}
		if(marginScale!=1 || short40.marginScale!=.985 ||
			!"70".equals(options.get("alignmentminspan")) || !"28".equals(short40.options.get("alignmentminspan"))){
			throw new IllegalArgumentException("Shared length policy requires the fitted LONG70/margin1 and "+
				"SHORT28/margin0.985 base rules");
		}
		for(StopTailCascadeResources view:new StopTailCascadeResources[]{this,short40}){
			if(view.fitted!=null || view.lengthScaledAdmission || view.familyMeanAdmission || view.familyMeans!=null ||
				!"50".equals(view.options.get("top")) || !"5".equals(view.options.get("ddlreject")) ||
				!"hitfourth_h020_c150".equals(view.heuristicArm)){
				throw new IllegalArgumentException("Shared length policy requires unchanged fixed admission "+
					"and frozen heuristic ranking in both views");
			}
		}
	}
	public boolean hasShortBins(){return short40!=null || short100!=null || short150!=null;}
	/** Read-only development access to the explicitly bound independent40..99 child.
	 * No child inference, overrides or substitute package are permitted.
	 */
	public FamilyTailCascade.LoadedClassifier loadShort40(final int threads)throws Exception{
		if(threads<1 || short40==null || short40.minQueryAa!=40 || short40.maxQueryAa!=99){
			throw new IllegalArgumentException("Development SHORT loading requires an explicit40..99 child "+
				"and positive threads");
		}
		verifyManifest();
		final FamilyTailCascade.LoadedClassifier loaded=short40.load(threads);
		verifyManifest();return loaded;
	}
	/** Configured query ranges; longer intervals are clipped to the explicitly bound root cap. */
	public String queryLengths(){return minQueryAa==100 && maxQueryAa==100 ?
		(short40==null ? "100" : shortMinimum+"-99,100") :
		(short40==null ? "" : shortMinimum+"-99,")+(short100==null ? "" : "100-"+short100.maxQueryAa+",")+
		(short150==null ? "" : "150-199,")+(minQueryAa==120 && maxQueryAa==120 ? "120" :
			minQueryAa==150 && maxQueryAa==150 ? "150" : "200");}

	/** Call before opening any caller output; existing hard links are checked as well as canonical names. */
	public void protectOutputs(final String[] outputs) throws IOException{
		if(outputs==null){throw new IllegalArgumentException("Caller output roster is required for resource protection");}
		if(short40!=null){short40.protectOutputs(outputs);}
		if(short100!=null){short100.protectOutputs(outputs);}
		if(short150!=null){short150.protectOutputs(outputs);}
		for(String name:outputs){
			if(!Tools.isOutputFileName(name)){continue;}
			final Path out=Paths.get(name).toFile().getCanonicalFile().toPath();
			final Path protectedRoot=compactPackage==null ? manifest.getParent() : manifest.getParent().getParent();
			if(out.startsWith(protectedRoot)){
				throw new IllegalArgumentException("Output would be created inside stop-tail resources: "+out);
			}
			final boolean exists=Files.exists(out);
			for(Path input:inputs){
				if(out.equals(input) || exists && Files.isSameFile(out, input)){
					throw new IllegalArgumentException("Output aliases stop-tail resource: "+out+" -> "+input);
				}
			}
		}
	}

	/** Resolved payload location for reporting, or null for an omitted v2 blacklist. */
	public String path(final String key){
		if("familymedians".equals(key)){return familyMedians==null ? null : familyMedians.path.toString();}
		if("familyhydro".equals(key)){return familyHydro==null ? null : familyHydro.path.toString();}
		if("familymeans".equals(key)){return familyMeans==null ? null : familyMeans.path.toString();}
		if(!contains(PATH_KEYS, key)){throw new IllegalArgumentException("Unknown classifier payload role: "+key);}
		return options.get(key);
	}

	/** Detached resource roles, including child-local policy inputs and profile dependencies.
	 * Uses this reader's existing discovery order; never loads classifier arrays. */
	public Map<String, Path> resolvedResourceBindings()throws IOException{
		verifyManifest();
		final LinkedHashMap<String, Path> result=new LinkedHashMap<String, Path>();
		resourceBindings("", result);
		verifyManifest();
		return result;
	}
	/** Append qualified roles from this reader's own resolved fields and tracked dependency order. */
	private void resourceBindings(String prefix, LinkedHashMap<String, Path> result){
		final LinkedHashMap<String, Path> local=new LinkedHashMap<String, Path>();
		local.put("manifest", manifest);
		for(String key:PATH_KEYS){
			if(options.containsKey(key)){
				local.put(key, Paths.get(options.get(key)));
			}
		}
		if(familyMedians!=null){local.put("familymedians", familyMedians.path);}
		if(familyHydro!=null){local.put("familyhydro", familyHydro.path);}
		if(familyMeans!=null){local.put("familymeans", familyMeans.path);}
		if(evaluationMatrix!=null){local.put("uncalibratedmatrixfile", evaluationMatrix.path);}
		if(gc!=null){
			local.put("familygc", gc.stats.path);
			local.put("familygcmap", gc.map.path);
		}
		// Remaining tracked inputs are profile members/sibling summaries, or compact resources.
		// Their stable ordinal is the reader's own discovery order, not a second path parser.
		int ordinal=0;
		for(Path input:inputs){
			if(!local.containsValue(input)){
				local.put("dependency."+(ordinal++), input);
			}
		}
		for(Map.Entry<String, Path> entry:local.entrySet()){result.put(prefix+entry.getKey(), entry.getValue());}
		if(short40!=null){short40.resourceBindings(prefix+"short40.", result);}
		if(short100!=null){short100.resourceBindings(prefix+"short100.", result);}
		if(short150!=null){short150.resourceBindings(prefix+"short150.", result);}
	}

	/** Recheck the complete observed resource graph in structural mode. */
	public void verify()throws IOException{verifyManifest();}
	private void verifyManifest()throws IOException{
		verifyEvaluationMatrices();
		if(compactPackage!=null){
			try{compactPackage.verify();}
			catch(RuntimeException e){throw e;}catch(IOException e){throw e;}catch(Exception e){
				throw new IOException("Compact caller input verification failed: "+manifest,e);
			}
			if(short40!=null){short40.verifyManifest();}
		}
		if(structuralInputs!=null){structuralInputs.verify();return;}
		if(!sha80.equals(DigestSuffix.file(manifest.toString()))){
			throw new IllegalStateException("Classifier manifest changed after it was read: "+manifest);
		}
	}
	/** Small diagnostic files are checked across every child, including strict root verification. */
	private void verifyEvaluationMatrices(){
		if(evaluationMatrix!=null && (Files.isSymbolicLink(evaluationMatrix.path) ||
			!evaluationMatrix.pin.equals(DigestSuffix.file(evaluationMatrix.path.toString())))){
			throw new IllegalStateException("Uncalibrated evaluation matrix changed after its explicit byte binding");
		}
		if(short40!=null){short40.verifyEvaluationMatrices();}
		if(short100!=null){short100.verifyEvaluationMatrices();}
		if(short150!=null){short150.verifyEvaluationMatrices();}
	}
	private static Path resolve(final Path parent, final String value) throws IOException{
		assert(parent!=null && value!=null) : "Manifest paths require an explicit parent and retained value";
		final Path path=Paths.get(value);
		// Keep the declared parent for profile-list relative entries and sibling summaries.
		final Path resolved=(path.isAbsolute() ? path : parent.resolve(path)).toAbsolutePath().normalize();
		if(!Files.isRegularFile(resolved)){throw new IllegalArgumentException("Missing classifier resource file: "+resolved);}
		return resolved;
	}
	private static boolean contains(final String[] values, final String key){
		for(String value:values){if(value.equals(key)){return true;}}
		return false;
	}
	/** Optional inputs are paired, pinned and protected exactly like mandatory resources. */
	private static final class GcPolicy{
		GcPolicy(Payload stats_,Payload map_,int families_,long members_,double lambda_,boolean filter_,double z_){
			stats=stats_;map=map_;families=families_;members=members_;lambda=lambda_;filter=filter_;z=z_;
		}
		static GcPolicy read(HashMap<String,String> options,Path parent,ArrayList<Path> inputs,
				StructuralInputs structural,int minAa,int maxAa)throws IOException{
			boolean present=false;for(String key:GC_KEYS){present|=options.containsKey(key);}
			if(!present){return null;}
			for(String key:GC_KEYS){
				if(!options.containsKey(key)){
					throw new IllegalArgumentException("Incomplete caller GC binding: "+key);
				}
			}
			if(structural==null || !(minAa==40 && maxAa==99 || minAa==100 && maxAa==100)){
				throw new IllegalArgumentException("Caller GC experiments require structural SHORT40..99 or LONG100 resources");
			}
			final Payload stats=Payload.read(options,parent,inputs,"familygc",structural);
			final Payload map=Payload.read(options,parent,inputs,"familygcmap",structural);
			final int families=Integer.parseInt(options.remove("familygcfamilies"));
			final long members=Long.parseLong(options.remove("familygcmembers"));
			final double lambda=parseMarginScale(options.remove("familygclambda"));
			final double z=parseMarginScale(options.remove("familygcz"));
			final String filter=options.remove("familygcfilter");
			if(families<1 || members<families || !(filter.equals("t") || filter.equals("f"))){
				throw new IllegalArgumentException("GC binding requires positive full-roster counts and filter=t/f");
			}
			return new GcPolicy(stats,map,families,members,lambda,filter.equals("t"),z);
		}
		final Payload stats,map;final int families;final long members;final double lambda,z;final boolean filter;
	}

	/** Optional inputs are paired, pinned and protected exactly like mandatory resources. */
	private static final class Payload{
		private Payload(Path path_,String pin_){path=path_;pin=pin_;}
		static Payload read(HashMap<String,String> options,Path parent,ArrayList<Path> inputs,String key) throws IOException{
			return read(options,parent,inputs,key,null);
		}
		static Payload read(HashMap<String,String> options,Path parent,ArrayList<Path> inputs,
				String key,StructuralInputs structural) throws IOException{
			assert(options!=null && parent!=null && inputs!=null && key!=null) :
				"Optional payloads need the owning manifest and its output-protection roster";
			final String name=options.remove(key),pin=options.remove(key+"sha80");
			if(name==null && pin==null){return null;}
			if(name==null || !declaredIdentity(pin,structural!=null)){
				throw new IllegalArgumentException(key+" and "+key+"sha80 require a paired path and sha80");
			}
			final Path path=resolve(parent,name);
			if(Files.isSymbolicLink(path) || structural==null && !DigestSuffix.file(path.toString()).equals(pin)){
				throw new IllegalArgumentException("Optional classifier payload is a symlink or differs from "+
					"its pin: "+key+" "+path);
			}
			if(structural!=null){structural.track(path);}
			inputs.add(path);return new Payload(path,pin);
		}
		final Path path;final String pin;
	}
	/** Version 2 must bind the selected operating point explicitly; version 1 always uses scale 1. */
	private static double parseMarginScale(final String value){
		if(value==null){throw new IllegalArgumentException("Missing explicit classifier setting: marginscale");}
		final double scale;
		try{scale=Double.parseDouble(value);}catch(NumberFormatException e){
			throw new IllegalArgumentException("marginscale must be a finite nonnegative double: "+value, e);
		}
		if(!Double.isFinite(scale) || scale<0){
			throw new IllegalArgumentException("marginscale must be a finite nonnegative double: "+value);
		}
		return scale;
	}

	private static void validateOverrides(int threads,int top,int minHits,float reduction){
		if(threads<1 || top<0 || minHits<-1 || !Float.isFinite(reduction) ||
			reduction<0 || reduction>100){
			throw new IllegalArgumentException("Caller overrides require positive workers, optional positive cap, "+
				"optional nonnegative floor and a finite 0..100 percentage-point lower reduction");
		}
	}
	private static boolean declaredIdentity(String value,boolean structural){
		return value!=null && (value.matches("[0-9a-f]{20}") || structural && value.equals("unmeasured"));
	}
	/** Unmeasured labels cannot establish index equality; structural views must share the actual file. */
	private static boolean sameProfiles(HashMap<String,String> full,StopTailCascadeResources child,
			StructuralInputs structural)throws IOException{
		return structural==null ? full.get("ddlprofilessha80").equals(child.options.get("ddlprofilessha80")) :
			Files.isSameFile(Paths.get(full.get("ddlprofiles")),Paths.get(child.options.get("ddlprofiles")));
	}

	public final Path manifest;
	public final String sha80;
	/** Absolute acceptance-margin multiplier; never changes member rejection, ranking or reference bytes. */
	public final double marginScale;
	public final int minQueryAa, maxQueryAa;
	/** Immutable scanner boundary; never changes reference arrays or acceptance gates. */
	public final boolean firstStart;
	/** Explicit policy of this child, not an inherited root-wide admission override. */
	public final boolean lengthScaledAdmission;
	/** Family-specific raw-hit floors from the pinned reference member census. */
	public final boolean familyMeanAdmission;
	/** Two-set caller search, default for version5; explicit false retains single-set routing. */
	public final boolean allFamilySets;
	/** Explicit optional arm, or null for the unchanged classifier ranking policy. */
	public final String heuristicArm;
	private final Payload familyMedians,familyHydro,familyMeans;
	/** Separate diagnostic provenance; options still describe the calibrated reference metric. */
	private final Payload evaluationMatrix;
	private final GcPolicy gc;
	private final FittedPolicy fitted;
	private final StructuralInputs structuralInputs;
	private final FamilyScoringPackage compactPackage;
	private final int shortMinimum;
	private final HashMap<String,String> options;
	private final ArrayList<Path> inputs;
	private final StopTailCascadeResources short40, short100, short150;
	private static final String HEADER="#bbgap_tail_cascade\t1";
	private static final String FITTED_HEADER="#bbgap_tail_fitted\t1";
	private static final String[] FITTED_KEYS={"parent","parentsha80","cohorts","cohortssha80",
		"floors","floorssha80","marginscale"};
	private static final String HEADER_V2="#bbgap_tail_cascade\t2";
	private static final String HEADER_V3="#bbgap_tail_cascade\t3";
	private static final String HEADER_V4="#bbgap_tail_cascade\t4";
	private static final String HEADER_V5="#bbgap_tail_cascade\t5";
	private static final String HEADER_V6="#bbgap_tail_cascade\t6";
	private static final String HEADER_V7="#bbgap_tail_cascade\t7";
	private static final String HEADER_V8="#bbgap_tail_cascade\t8";
	private static final String[] SHORT_KEYS={"short100manifest", "short100sha80", "short150manifest", "short150sha80"};
	private static final String[] PATH_KEYS={"familyreferences", "ddlprofiles", "blacklist", "scalarbounds"};
	private static final String[] RANKING_KEYS={"heuristicarm","familymedians","familymedianssha80",
		"familyhydro","familyhydrosha80","familymeans","familymeanssha80"};
	private static final String[] METRIC_KEYS={"alignmentmode","alignmentminspan","alignmentminpaired",
		"matrix","matrixbinding"};
	private static final String[] GC_KEYS={"familygc","familygcsha80","familygcmap","familygcmapsha80",
		"familygcfamilies","familygcmembers","familygclambda","familygcfilter","familygcz"};
	private static final String[] DIAGNOSTIC_MATRIX_KEYS={"uncalibratedmatrixfile","uncalibratedmatrixfilesha80"};
	private static final String[] GEOMETRY_KEYS={"aabits","aapartition","positionmode"};
	private static final String[] BLACKLIST_POLICY_KEYS={"blacklistk","blacklistmode","blacklistwords",
		"blacklistprovenance"};
	private static final String[] KEYS={"familyreferences", "familyreferencessha80", "top", "k", "seed", "width",
		"ddlreject", "ddlaccept", "cosine", "rankingfeatures", "scalarbounds", "scalarboundssha80", "scalarscale",
		"scalargates", "scalarz", "identityacceptfloor", "blosumacceptfloor", "positional", "ddlprofiles", "ddlprofilessha80",
		"partitions", "overlap", "blacklist", "blacklistsha80"};
}
