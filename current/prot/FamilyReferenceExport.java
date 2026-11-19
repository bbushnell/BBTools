package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.LinkOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Locale;

import ddl.DDLPentamerBlacklist;
import ddl.DDLFamilyProfiles;
import ddl.ProteinKmerGeometry;
import ddl.DDLFamilySearch.StructuralInputs;
import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.LongHashSet;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Converts pinned all-member reference/calibration collections without biological rescoring.
 * Keeps raw inputs untouched and publishes a PASS only after every input is reverified.
 * @author Sayu, Keqing
 */
public final class FamilyReferenceExport{
	private FamilyReferenceExport(){}
	public static void main(String[] args){
		try{
			final PreParser pp=new PreParser(args,FamilyReferenceExport.class,false);
			try{run(pp.args);}finally{Shared.closeStream(pp.outstream);}
		}catch(Throwable failure){failure.printStackTrace();System.exit(1);}
	}
	public static void run(String[] args) throws Exception{
		final HashMap<String,String> options=new HashMap<String,String>();
		for(String arg:args){
			final String[] pair=Parse.splitOnFirst(arg,'='); require(pair.length==2,"Require flag=value");
			final String key=pair[0].toLowerCase(Locale.ROOT);
			require(Arrays.asList(OPTIONS).contains(key) && options.put(key,pair[1])==null,"Unknown/duplicate option: "+key);
		}
		final int limit=options.containsKey("maxfamilies") ? Parse.parseIntKMG(options.get("maxfamilies")) : 1000000;
		require(limit>0 && limit<=1000000,"maxfamilies must be 1..1000000");
		if(options.containsKey("mode")){
			if(options.get("mode").equals("compact")){compact(options);return;}
			if(options.get("mode").equals("monomer")){monomer(options,limit);return;}
			if(options.get("mode").equals("subset")){subset(options,limit);return;}
			// Experimental November 16 2026 exports stay behind explicit modes.
			// Member-start M1 added +39/+350 ALL/40-99 raw objective points
			// at20/40 and +7/+333 at40/80 in one seed; fuller +1 padding added
			// positives in all six paired cells but raised4/100 negatives. Neither was adopted.
			// See results/member_start_fit_20261116 and
			// results/identity_padding_fuller_replays_20261116.
			if(options.get("mode").equals("identitypadding")){identityPadding(options,limit);return;}
			if(options.get("mode").equals("memberstart")){memberStart(options,limit);return;}
			require(options.get("mode").equals("raise"),"Existing references require mode=compact, mode=raise, mode=memberstart, mode=monomer, mode=subset or mode=identitypadding");
			raise(options,limit);return;
		}
		for(String key:new String[]{"parent","parentsha80","overrides","overridessha80","overrideslist","overrideslistsha80","identityfloor","threads",
			"rawcalibration","rawcalibrationsha80",
			"survivors","survivorssha80","ddlprofiles","ddlprofilessha80"}){
			require(!options.containsKey(key),"Parent-reference option requires an explicit mode: "+key);
		}
		final Pins pins=new Pins();
		final Path refIndex=input(options,"referenceindex",pins),upperIndex=input(options,"upperindex",pins),consensus=input(options,"consensus",pins);
		require(options.containsKey("supplementalupper")==options.containsKey("supplementaluppersha80"),"Supplemental upper table and pin must be supplied together");
		final Path supplemental=options.containsKey("supplementalupper") ? input(options,"supplementalupper",pins) : null;
		final boolean additional=options.containsKey("additionalupper");
		require(additional==options.containsKey("additionaluppersha80") && additional==options.containsKey("additionalqueries"),
			"Additional upper table, pin and original query population must be supplied together");
		final Path additionalUpper=additional ? input(options,"additionalupper",pins) : null;
		final int additionalQueries=additional ? Parse.parseIntKMG(required(options,"additionalqueries")) : 0;
		require(!additional || additionalQueries>0,"Additional negative sample must have a positive query population");
		final boolean listed=options.containsKey("additionalupperlist");
		require(listed==options.containsKey("additionalupperlistsha80") && !(listed && additional),
			"Use either a pinned additional-upper list or one pinned table/population, not both");
		final Path additionalList=listed ? input(options,"additionalupperlist",pins) : null;
		final double blosumFloor=options.containsKey("blosumfloor") ? Double.parseDouble(options.get("blosumfloor")) : 0;
		require(Double.isFinite(blosumFloor) && blosumFloor>=0 && blosumFloor<=Float.MAX_VALUE,"BLOSUM floor must be finite, nonnegative and representable");
		final double identityUpperMultiplier=identityPadding(options,"identityuppermultiplier",
			FamilyGateCodec.DEFAULT_IDENTITY_UPPER_MULTIPLIER,1);
		final double identityUpperAdditive=identityPadding(options,"identityupperadditive",
			FamilyGateCodec.DEFAULT_IDENTITY_UPPER_ADDITIVE,0);
		final boolean customIdentityPadding=identityUpperMultiplier!=FamilyGateCodec.DEFAULT_IDENTITY_UPPER_MULTIPLIER ||
			identityUpperAdditive!=FamilyGateCodec.DEFAULT_IDENTITY_UPPER_ADDITIVE;
		final boolean scalar=options.containsKey("scalars");
		require(scalar==options.containsKey("scalarssha80") && scalar==options.containsKey("scalarz"),
			"Scalar table, pin and z must be supplied together");
		final Path scalars=scalar ? input(options,"scalars",pins) : null;
		final double scalarZ=scalar ? Double.parseDouble(required(options,"scalarz")) : 0;
		require(!scalar || Double.isFinite(scalarZ) && scalarZ>0,"Scalar z must be finite and positive");
		require(options.containsKey("floorreference")==options.containsKey("floorreferencesha80"),
			"Historical floor reference and pin must be supplied together");
		final Path floorReference=options.containsKey("floorreference") ? input(options,"floorreference",pins) : null;
		final boolean replacement=options.containsKey("alignmentupper");
		require(replacement==options.containsKey("alignmentuppersha80") && replacement==options.containsKey("alignmentqueries"),
			"Complete alignment-upper replacement table, pin and query population must be supplied together");
		final Path alignmentUpper=replacement ? input(options,"alignmentupper",pins) : null;
		final int alignmentQueries=replacement ? Parse.parseIntKMG(required(options,"alignmentqueries")) : 0;
		require(!replacement || alignmentQueries>0,"Complete alignment upper replacement must declare a positive query population");
		final Path out=absolute(required(options,"out"));
		require(!Files.exists(out) && !Files.isSymbolicLink(out) && out.getParent()!=null &&
			Files.isDirectory(out.getParent()) && out.getParent().toRealPath().equals(out.getParent()),"Require fresh output directory with real parent");
		final Bucket[] refs=index(refIndex,false),uppers=index(upperIndex,true);
		final Entry[] entries=new Entry[limit]; final HashMap<String,Entry> byId=new HashMap<String,Entry>();
		int families=0,live=0,empty=0;
		for(int b=0;b<256;b++){
			final Bucket r=refs[b],u=uppers[b];
			if(r==null){require(u==null,"Upper bucket without reference");continue;}
			require(u!=null && r.same(u),"Reference/upper bucket census mismatch: "+b);
			families=Math.addExact(families,r.families);live=Math.addExact(live,r.live);empty=Math.addExact(empty,r.empty);
			require(families<=limit,"Family census exceeds allocation bound");
			pins.add(r.manifest,r.pin);pins.add(r.lowerManifest,r.lowerPin);pins.add(u.manifest,u.pin);
			final HashMap<String,Path> files=manifest(r.manifest,REFERENCE_FILES,"FAMILY_FULL_REFERENCES_PASS",pins);
			final HashMap<String,Path> lower=manifest(r.lowerManifest,LOWER_FILES,"FAMILY_ALL_MEMBER_CALIBRATION_PASS",pins);
			final HashMap<String,Path> upper=manifest(u.manifest,UPPER_FILES,"FAMILY_NEGATIVE_UPPER_CALIBRATION_PASS",pins);
			validateSummary(files.get("summary.tsv"),r,b);
			final int upperQueries=upperQueries(upper.get("summary.tsv"),u);
			census(files.get("members.tsv"),r,b,entries,byId);
			profiles(files.get("profiles.tsv"),b,byId);
			final HashMap<String,ProteinSpectrum.Spectrum> spectra=ProteinSpectrum.load(files.get("spectra.tsv").toString());
			final ArrayList<FamilyGateExport.Row> rows=new ArrayList<FamilyGateExport.Row>();
			final HashMap<String,FamilyGateExport.Row> gates=new HashMap<String,FamilyGateExport.Row>();
			FamilyGateExport.readLower(lower.get("cutoffs.tsv"),rows,gates,limit);
			FamilyGateExport.readUpper(upper.get("uppers.tsv"),gates,upperQueries,500);
			require(rows.size()==r.families,"Lower/reference family count mismatch");
			int seenLive=0;
			for(FamilyGateExport.Row gate:rows){
				final Entry e=byId.get(gate.id);
				require(e!=null && e.bucket==b && e.rank==gate.rank && e.members==gate.members && e.gates==null,"Unmatched/duplicate lower family");
				e.gates=gate;
				if(e.members==0){continue;}
				require(gate.hasUpper && e.ddl!=null && e.profileModes==3,"Incomplete live reference/calibration");
				final ProteinSpectrum.Spectrum s=spectra.remove("tail:"+e.id);
				require(s!=null && s.k==2,"Missing dimer spectrum: "+e.id);
				e.spectrum=frequencies(s);seenLive++;
			}
			require(seenLive==r.live && spectra.isEmpty(),"Spectrum/live census mismatch");
		}
		require(families>0 && families==live+empty && byId.size()==families,"Incomplete collection census");
		for(int rank=0;rank<entries.length;rank++){
			require(rank<families ? entries[rank]!=null && entries[rank].gates!=null : entries[rank]==null,"Noncontiguous roster rank");
		}
		final long[] alignment=alignmentUpper==null ? new long[2] : replaceAlignmentUpper(alignmentUpper,entries,families,live,alignmentQueries);
		final long[] supplement=supplemental==null ? new long[2] : mergeSupplemental(supplemental,entries,families,0);
		final long[] extra=listed ? mergeAdditionalList(additionalList,entries,families,pins) :
			additionalUpper==null ? new long[2] : mergeSupplemental(additionalUpper,entries,families,additionalQueries);
		readConsensus(consensus,byId,families);
		final FamilyReference[] records=new FamilyReference[families];
		for(int rank=0;rank<families;rank++){
			final Entry e=entries[rank];final FamilyGateExport.Row g=e.gates;
			if(e.members==0){records[rank]=FamilyReference.empty(e.id,rank);}
			else{
				records[rank]=new FamilyReference(e.id,rank,e.members,e.residues,e.words,new ProteinSequence(e.id,e.consensus),e.ddl,e.spectrum,
					FamilyGateCodec.lower(g.lower[0]),FamilyGateCodec.lower(g.lower[1]),FamilyGateCodec.lower(g.lower[2]),
					FamilyGateCodec.cosineAccept(g.lower[0],g.upper[0]),
					FamilyGateCodec.upper(g.lower[1],g.upper[1],identityUpperMultiplier,identityUpperAdditive),
					FamilyGateCodec.finalBlosum(Math.max(g.lower[2],blosumFloor),g.upper[2]),e.consensus);
			}
		}
		if(scalar){applyScalarBounds(scalars,records,scalarZ);}
		final int[] floors=floorReference==null ? null : applyHistoricalFloors(floorReference,records,1,limit);
		pins.verify();Files.createDirectory(out);
		final FamilyReferenceCodec.Header header=new FamilyReferenceCodec.Header(FamilyReferenceCodec.SCHEMA,5,1024,0,families,live,empty,
			FamilyReferenceCodec.FORMULA_VERSION,FamilyReferenceCodec.SPECTRUM_ENCODING,scalarZ);
		final ByteStreamWriter writer=new ByteStreamWriter(out.resolve("references.tsv").toString(),false,false,false);writer.start();
		final ByteBuilder buffer=new ByteBuilder();
		try{
			writer.print(FamilyReferenceCodec.appendHeader(header,buffer));
			for(FamilyReference record:records){writer.print(FamilyReferenceCodec.appendFiveLines(record,buffer.clear()));}
		}finally{require(!writer.poisonAndWait(),"Reference writer failed; no PASS");}
		write(out.resolve("recipe.txt"),"schema="+FamilyReferenceCodec.SCHEMA+"\nformula_version="+FamilyReferenceCodec.FORMULA_VERSION+
			"\ncosine_formula=round_up_float32(max(raw_lower,0.72,1.1*raw_upper,raw_upper+0.05));all_live_finite;no_cap;strict_greater\n"+
			"blosum_formula=round_up_float32(max(lower,blosum_floor,1.1*upper,upper+0.005))\nblosum_floor="+blosumFloor+"\n"+
			"population=all_winner_assigned_members\nprofile=k5_width1024_legacy16_exponent6_seed0_nonempty\n"+
			"spectrum=(float)((double)count/total);zero_total_allzero;A48_raw_float_bits\n"+
			"norm=float_sqrt(ascending_float_sum_of_squared_frequencies)\norder=roster_rank\nempty=four_NA_payload_lines\n"+
			(customIdentityPadding ? "identity_formula=round_up_float32(max(lower,"+identityUpperMultiplier+"*upper,upper+"+
				identityUpperAdditive+",30));strict_greater\nidentity_upper_multiplier="+identityUpperMultiplier+
				"\nidentity_upper_additive="+identityUpperAdditive+"\n" :
				"identity_formula=round_up_float32(max(lower,1.1*upper,upper+2,30));strict_greater\n")+
			"empty_numeric_fields=finite_zero;not_measured;excluded_from_live_roster\n"+
			"computed_query_metrics=float\nrecalibration=false\n"+
			"alignment_upper_replacement="+(alignmentUpper==null ? "none" : "complete_identity_blosum_raw_maxima")+"\n"+
			"alignment_queries="+alignmentQueries+"\nalignment_families="+alignment[0]+"\nalignment_pair_observations="+alignment[1]+"\n"+
			"legacy_cosine_upper="+(alignmentUpper==null ? "active" : "retained_from_upperindex;not_newly_calibrated;requires_cosine_ranking_downstream")+"\n"+
			"raw_upper_base="+(alignmentUpper==null ? "original_top500" : "complete_alignment_replacement")+"\n"+
			"raw_upper_union="+(supplemental==null ? "base_only" : "max(base,deployed_shortlist_supplement)")+"\n"+
			"supplemental_families="+supplement[0]+"\nsupplemental_pair_observations="+supplement[1]+"\n"+
			"additional_upper_union=max(existing_raw_upper,observed_additional_upper);absent_family_retains_existing_evidence\n"+
			"additional_queries="+(listed ? "per_table_in_pinned_manifest" : Integer.toString(additionalQueries))+"\n"+
			"additional_tables="+(listed ? extra[2] : additional ? 1 : 0)+"\nadditional_families="+extra[0]+"\nadditional_pair_observations="+extra[1]+"\n"+
			"additional_observation_count=sum_across_sources;overlapping_query_family_pairs_not_deduplicated\n"+
			"unique_union_pair_count=not_computed;old_top500_membership_not_available\n"+
			(floorReference==null ? "" : "historical_floor_reference="+floorReference+"\n"+
				"historical_floor_families="+floors[0]+"\nhistorical_identity_raised="+floors[1]+"\nhistorical_blosum_raised="+floors[2]+"\n"+
				"historical_gate_rule=max(new_encoded_gate,historical_final_gate);no_gate_fell\n"));
		write(out.resolve("scalar_recipe.txt"),"scalar_z="+scalarZ+"\nscalar_formula=mean_plus_minus_z_times_population_sd;inclusive_bounds\n"+
			"scalar_order=hydrophobicity,net_charge,aromaticity,fraction_small,pI_tail\n"+
			"scalar_encoding=decimal_double;interleaved_min_max;five_bit_enabled_mask\n"+
			"scalar_disabled=N_lt_30_or_sd_zero;mask_bit_clear_and_both_bounds_zero\n"+
			"scalar_action=reject_shortlisted_pair_only;never_accept;undefined_query_pI_skips_pI\n"+
			"retention=per_scalar_gaussian_projection;not_empirical_joint_retention\n");
		writeGateQuantiles(out,entries,records,live);
		final ByteBuilder inputs=new ByteBuilder("path\tsha80\n");
		for(java.util.Map.Entry<Path,String> p:pins.values.entrySet()){inputs.append(p.getKey().toString()).tab().append(p.getValue()).nl();}
		write(out.resolve("inputs.tsv"),inputs.toString());pins.verify();
		write(out.resolve("PASS"),PASS+"\n");
		final ByteBuilder result=new ByteBuilder("file\tsha80\n");
		for(String name:new String[]{"references.tsv","recipe.txt","scalar_recipe.txt","gate_quantiles.tsv","inputs.tsv","PASS"}){result.append(name).tab().append(DigestSuffix.file(out.resolve(name).toString())).nl();}
		write(out.resolve("result.sha80.tsv"),result.toString());
		System.out.println(PASS+" families="+families+" live="+live+" empty="+empty);
	}

	/** Joins the accepted compact side tables and renders ten-line family payloads. */
	private static void compact(HashMap<String,String> options) throws Exception{
		compact(options,compactCensus(required(options,"view")));
	}

	private static void compact(HashMap<String,String> options,CompactCensus census) throws Exception{
		for(String key:options.keySet()){
			require(Arrays.asList(COMPACT_OPTIONS).contains(key),"Compact binding cannot rederive or render payloads from option: "+key);
		}
		final StructuralInputs inputs=new StructuralInputs();
		final Path references=structuralInput(options,"references",inputs),
			ddlProfiles=structuralInput(options,"ddlprofiles",inputs),
			profileVectors=structuralInput(options,"profilevectors",inputs),
			profileSummary=structuralInput(options,"profilesummary",inputs),
			lengthMedians=structuralInput(options,"lengthmedians",inputs),
			familyHydro=structuralInput(options,"familyhydro",inputs),
			familyGc=structuralInput(options,"familygc",inputs),
			familyGcMap=structuralInput(options,"familygcmap",inputs),
			scalarBounds=structuralInput(options,"scalarbounds",inputs),
			alphabet=structuralInput(options,"alphabet",inputs),
			blacklist=structuralInput(options,"blacklist",inputs),
			neutralBlacklist=structuralInput(options,"neutralblacklist",inputs),
			matrixResource=structuralInput(options,"matrixresource",inputs),
			out=absolute(required(options,"out"));
		final boolean explicitFamilies=options.containsKey("outfamilies");
		final Path familiesOut=explicitFamilies ? absolute(required(options,"outfamilies")) :
			out.resolve("families.tsv");
		final int zipThreads=options.containsKey("threads") ?
			Parse.parseIntKMG(options.get("threads")) : Shared.threads();
		require(zipThreads>0,"Compact BGZF writer requires a positive thread count");
		require(!Files.exists(out) && !Files.isSymbolicLink(out) && out.getParent()!=null &&
			Files.isDirectory(out.getParent()) && out.getParent().toRealPath().equals(out.getParent()),
			"Require fresh compact output directory with real parent");
		if(explicitFamilies){
			require(!Files.exists(familiesOut) && !Files.isSymbolicLink(familiesOut) && familiesOut.getParent()!=null &&
				Files.isDirectory(familiesOut.getParent()) &&
				familiesOut.getParent().toRealPath().equals(familiesOut.getParent()),
				"Require fresh compact families output with real parent");
			require(familiesOut.toString().endsWith(".gz"),"Explicit compact family output must be BGZF .gz");
		}
		final ProteinKmerGeometry geometry=compactGeometry(compactAabits(alphabet));
		final String matrixBinding=ProteinSubstitutionMatrix.load(ProteinSubstitutionMatrix.POSSUM55_NAME,
			matrixResource.toString()).binding();
		require(required(options,"matrixbinding").equals(matrixBinding),
			"Compact matrix resource must match its exact POSSUM55 matrixbinding");
		final LinkedHashMap<String,String> fields=new LinkedHashMap<String,String>();
		fields.put("families",Integer.toString(census.families));
		fields.put("view",census.view);
		fields.put("alphabet",compactLocator(options,"alphabetlocator","alphabet.config"));
		fields.put("blacklist",compactLocator(options,"blacklistlocator","blacklist.faa.gz"));
		fields.put("matrixresource",compactLocator(options,"matrixlocator","POSSUM55.tsv"));
		fields.put("live",Integer.toString(census.live));
		fields.put("empty",Integer.toString(census.empty));
		fields.put("k","6");
		fields.put("width",Integer.toString(FamilyCompactCodec.DDL_WIDTH));
		fields.put("seed","0");
		fields.put("ddlreject","6");
		fields.put("ddlaccept","1025");
		fields.put("scalargates","true");
		fields.put("indexearlyreject","false");
		fields.put("namespace","tail");
		fields.put("mode","nonempty");
		fields.put("sketch","positional-v1");
		fields.put("exponent","6");
		fields.put("geometry",geometry.name());
		fields.put("partitions","4");
		fields.put("overlap","2");
		fields.put("tag_bits","2");
		fields.put("segment_words","49");
		fields.put("formulaVersion",FamilyReferenceCodec.FORMULA_VERSION);
		fields.put("rankingfeatures","combined25");
		fields.put("cosine","ranking");
		fields.put("matrix",ProteinSubstitutionMatrix.POSSUM55_NAME);
		fields.put("matrixbinding",matrixBinding);
		fields.put("matrixdisplay","POSSUM55");
		fields.put("matrixdescription",ProteinSubstitutionMatrix.POSSUM55_DESCRIPTION);
		fields.put("matrixprovenance","bbgap_possum55_20261116_v1");
		fields.put("top","50");
		fields.put("scalarscale","0.1");
		fields.put("scalarz","3.8906");
		fields.put("identityacceptfloor","30");
		fields.put("blosumacceptfloor","0.13");
		fields.put("blacklistk","6");
		fields.put("blacklistmode","reduced_literal");
		fields.put("blacklistwords","2048");
		fields.put("blacklistprovenance","bbgap_blacklist_k6_2048_20261116_v1");
		fields.put("alignmentmode","pin_right_trim_left_aligned_query");
		fields.put("alignmentminspan",Integer.toString(census.alignmentMinSpan));
		fields.put("alignmentminpaired",Integer.toString(census.alignmentMinPaired));
		fields.put("heuristic","hitfourth_h020_c150");
		fields.put("minqueryaa",Integer.toString(census.minAa));
		fields.put("maxqueryaa",Integer.toString(census.maxAa));
		fields.put("minmemberaa",Integer.toString(census.minAa));
		fields.put("maxmemberaa",Integer.toString(census.maxAa));
		fields.put("taillow",Integer.toString(census.minAa));
		fields.put("tailhigh",Integer.toString(census.maxAa));
		fields.put("marginscale",census.marginScale);
		fields.put("familygclambda","24");
		fields.put("gcfilter","true");
		fields.put("gc_z","8");
		fields.put("gc_minimum_members","30");
		fields.put("familygcfamilies",Integer.toString(census.live));
		fields.put("familygcmembers",Long.toString(census.members));
		fields.put("packageprovenance","bbgap_release_20261116_v1");
		addScalarBounds(fields,scalarBounds,census);
		validateCompactProfileList(ddlProfiles,profileVectors);
		final FamilyReferenceLoader.Loaded loaded=FamilyReferenceLoader.loadMetricBound(references.toString(),1,census.families);
		require(loaded.header.totalFamilies==census.families && loaded.header.liveFamilies==census.live &&
			loaded.header.emptyFamilies==census.empty && loaded.members==census.members,
			"Compact references disagree with the accepted view census");
		require(loaded.header.width==FamilyCompactCodec.DDL_WIDTH &&
			loaded.header.spectrumWidth==FamilyCompactCodec.MONOMER_WIDTH,
			"Compact references must already contain 512-register DDLs and 20 monomer frequencies");
		final FamilyReference[] records=loaded.all;
		final HashMap<String,FamilyReference> byId=new HashMap<String,FamilyReference>(records.length*2);
		for(FamilyReference f:records){
			require(f.rank<records.length && records[f.rank]==f &&
				byId.put(f.id,f)==null,"Compact reference ranks and family IDs must be unique");
		}
		final char[][] compactProfiles=compactProfiles(profileVectors,records,byId);
		validateCompactProfileSummary(profileSummary,fields,census,loaded.validWords,
			blacklist,neutralBlacklist,geometry);
		final double[] lengthMedian=new double[records.length],kdMean=new double[records.length],
			gcMean=new double[records.length],gcSd=new double[records.length];
		final long[] gcMeasured=new long[records.length];
		final int[] stable=new int[records.length];
		Arrays.fill(stable,-1);
		readCompactLengthMedians(lengthMedians,records,byId,lengthMedian);
		readCompactHydrophobicities(familyHydro,records,kdMean);
		readCompactGc(familyGc,familyGcMap,records,stable,gcMeasured,gcMean,gcSd,census);
		final FamilyCompactCodec.Header header=new FamilyCompactCodec.Header(fields);
		final ByteBuilder buffer=FamilyCompactCodec.appendHeader(header,new ByteBuilder());
		inputs.verify();
		if(explicitFamilies){ZipState.requireExternalBgzip();}
		Files.createDirectory(out);
		writeCompactFamilies(familiesOut,explicitFamilies,zipThreads,
			buffer,records,stable,lengthMedian,kdMean,gcMeasured,gcMean,gcSd,compactProfiles);
		write(out.resolve("recipe.txt"),"mode=compact\nview="+census.view+
			"\nreferences="+references+"\nddlprofiles="+ddlProfiles+"\nprofilevectors="+profileVectors+
			"\nprofilesummary="+profileSummary+"\nlengthmedians="+lengthMedians+
			"\nfamilyhydro="+familyHydro+"\nfamilygc="+familyGc+"\nfamilygcmap="+familyGcMap+
			"\nscalarbounds="+scalarBounds+"\nalphabet="+alphabet+"\nblacklist="+blacklist+
			"\nneutralblacklist="+neutralBlacklist+
			"\nmatrixresource="+matrixResource+
			"\noutput="+(explicitFamilies ? familiesOut.toString() : "families.tsv")+
			"\npayload_rendering=complete\nbiological_rescoring=false\n");
		write(out.resolve("inputs.tsv"),structuralInputRecords(new Path[]{references,ddlProfiles,profileVectors,
			profileSummary,lengthMedians,familyHydro,familyGc,familyGcMap,scalarBounds,alphabet,blacklist,
			neutralBlacklist,matrixResource}));
		inputs.verify();
		write(out.resolve("PASS"),COMPACT_PASS+"\n");
		System.out.println(COMPACT_PASS+" view="+census.view+" families="+census.families+" live="+census.live);
	}

	private static void writeCompactFamilies(Path out,boolean bgzip,int zipThreads,ByteBuilder buffer,
			FamilyReference[] records,int[] stable,double[] lengthMedian,double[] kdMean,
			long[] gcMeasured,double[] gcMean,double[] gcSd,char[][] compactProfiles){
		final ZipState zip=bgzip ? new ZipState() : null;
		ByteStreamWriter writer=null;
		try{
			if(zip!=null){zip.forceBgzip(zipThreads);}
			final FileFormat ff=FileFormat.testOutput(out.toString(),FileFormat.TXT,null,true,false,false,false);
			writer=new ByteStreamWriter(ff);
			writer.start();
			writer.print(buffer);
			for(FamilyReference f:records){
				final FamilyCompactCodec.Record record=f.members==0 ?
					new FamilyCompactCodec.Record(f.id,f.rank,f.rank,0,0,0,0,0,0,0,0,0,0,
						FamilyCompactCodec.AlignmentCalibration.UNAVAILABLE,new double[10],0,0,0,0,0,
						new byte[0],new char[FamilyCompactCodec.DDL_WIDTH],
						new float[FamilyCompactCodec.MONOMER_WIDTH]) :
					new FamilyCompactCodec.Record(f.id,f.rank,stable[f.rank],f.members,f.residues,
						f.validWords,f.cosineLower,f.identityLower,f.blosumLower,f.cosineAccept,
						f.identityUpper,f.blosumFinal,f.scalarMask,compactAlignment(f.alignmentCalibration),
						f.scalarBounds,lengthMedian[f.rank],kdMean[f.rank],gcMeasured[f.rank],
						gcMean[f.rank],gcSd[f.rank],f.consensusRaw,compactProfiles[f.rank],f.spectrum);
				writer.print(FamilyCompactCodec.appendRecord(record,buffer.clear()));
			}
		}finally{
			try{
				if(writer!=null){require(!writer.poisonAndWait(),"Compact family writer failed; no PASS");}
			}finally{
				if(zip!=null){zip.restore();}
			}
		}
	}

	private static void validateCompactProfileList(Path list,Path profileVectors){
		final ByteFile in=ByteFile.makeByteFile(list.toString(),false);
		final LineParser1 lp=new LineParser1('\t');
		int rows=0;
		try{
			for(byte[] line=in.nextLine();line!=null;line=in.nextLine()){
				lp.set(line);require(lp.terms()==3,"Compact DDL profile list requires path, sha80 and summary identity");
				Path path=Paths.get(lp.parseString(0));
				if(!path.isAbsolute()){path=list.getParent().resolve(path);}
				require(path.normalize().equals(profileVectors),"Compact DDL profile list must point at profilevectors");
				rows++;
			}
		}finally{require(!in.close(),"Compact DDL profile list read failed");}
		require(rows==1,"Compact package must bind exactly one DDL profile vector table");
	}

	private static char[][] compactProfiles(Path path,FamilyReference[] records,
			HashMap<String,FamilyReference> byId){
		final char[][] profiles=new char[records.length][];
		final ByteFile in=ByteFile.makeByteFile(path.toString(),false);
		final LineParser1 lp=new LineParser1('\t');
		int count=0;long words=0;
		try{
			header(in,DDLFamilyProfiles.HEADER16);
			for(byte[] line=in.nextLine();line!=null;line=in.nextLine()){
				lp.set(line);
				require(lp.terms()==10 && lp.termEquals("tail",0) &&
					integer(lp,2)==FamilyCompactCodec.DDL_WIDTH,"Compact profiles require 512-wide tail rows");
				require(lp.termEquals("literal",6) || lp.termEquals("nonempty",6),
					"Compact profile row has unknown mode");
				if(lp.termEquals("literal",6)){continue;}
				final FamilyReference f=byId.get(lp.parseString(1));
				require(f!=null && f.members>0 && profiles[f.rank]==null &&
					number(lp,3)==f.members && number(lp,4)==f.residues &&
					number(lp,5)==f.validWords,"Compact profile identity/count mismatch");
				final char[] values=new char[FamilyCompactCodec.DDL_WIDTH];
				decimalCharArray(lp,7,values);
				profiles[f.rank]=values;count++;words=Math.addExact(words,f.validWords);
			}
		}finally{require(!in.close(),"Compact profile vector read failed");}
		int expected=0;long expectedWords=0;
		for(FamilyReference f:records){
			if(f.members>0){expected++;expectedWords=Math.addExact(expectedWords,f.validWords);}
		}
		require(count==expected && words==expectedWords,"Compact profiles must cover every live family exactly once");
		return profiles;
	}

	private static void decimalCharArray(LineParser1 lp,int field,char[] target){
		lp.setBounds(field);
		final byte[] bytes=lp.line();int slot=0;long value=0;boolean digitSeen=false;
		for(int i=lp.a();i<=lp.b();i++){
			if(i==lp.b() || bytes[i]==','){
				require(digitSeen && slot<target.length && value<=65535,
					"Compact profile values must be 512 uint16 decimal terms");
				target[slot++]=(char)value;value=0;digitSeen=false;
			}else{
				final int digit=bytes[i]-'0';
				require(digit>=0 && digit<=9 && value<=(65535-digit)/10,
					"Invalid compact profile uint16 value");
				value=value*10+digit;digitSeen=true;
			}
		}
		require(slot==target.length,"Compact profile array must have 512 terms");
	}

	private static void validateCompactProfileSummary(Path path,LinkedHashMap<String,String> fields,
			CompactCensus census,long validWords,Path sourceBlacklist,Path neutralBlacklist,
			ProteinKmerGeometry geometry) throws Exception{
		final java.util.List<String> lines=Files.readAllLines(path,StandardCharsets.US_ASCII);
		require(lines.size()==2,"Compact profile summary requires one header and one data row");
		final String[] header=lines.get(0).split("\t",-1),values=lines.get(1).split("\t",-1);
		require(header.length==values.length,"Compact profile summary field count mismatch");
		final HashMap<String,String> map=new HashMap<String,String>();
		for(int i=0;i<header.length;i++){require(map.put(header[i],values[i])==null,"Duplicate profile summary key");}
		for(String key:new String[]{"namespace","k","seed","exponent","built_families","members",
				"valid_words","sketch","partitions","overlap","tag_bits","segment_words",
				"geometry","blacklist_words","blacklist_k","blacklist_mode","blacklist_provenance"}){
			require(map.containsKey(key),"Compact profile summary missing "+key);
		}
		require(map.get("namespace").equals(fields.get("namespace")) &&
			map.get("k").equals(fields.get("k")) && map.get("seed").equals(fields.get("seed")) &&
			map.get("exponent").equals(fields.get("exponent")) &&
			map.get("sketch").equals(fields.get("sketch")) &&
			map.get("partitions").equals(fields.get("partitions")) &&
			map.get("overlap").equals(fields.get("overlap")) &&
			map.get("tag_bits").equals(fields.get("tag_bits")) &&
			map.get("segment_words").equals(fields.get("segment_words")) &&
			map.get("geometry").equals(fields.get("geometry")) &&
			map.get("blacklist_words").equals(fields.get("blacklistwords")) &&
			map.get("blacklist_k").equals(fields.get("blacklistk")) &&
			map.get("blacklist_mode").equals(fields.get("blacklistmode")),
			"Compact profile summary metadata disagrees with the package header");
		require(Integer.toString(census.live).equals(map.get("built_families")) &&
			Long.toString(census.members).equals(map.get("members")) &&
			Long.toString(validWords).equals(map.get("valid_words")),
			"Compact profile summary census disagrees with joined references");
		final long[] sourceKeys=compactBlacklistKeys(sourceBlacklist,geometry,map.get("blacklist_provenance"));
		final long[] neutralKeys=compactBlacklistKeys(neutralBlacklist,geometry,fields.get("blacklistprovenance"));
		require(Arrays.equals(sourceKeys,neutralKeys),"Compact source and neutral blacklist loaded keys differ");
	}

	/** Validates compact blacklist policy/provenance, then returns the decoded reduced keys. */
	private static long[] compactBlacklistKeys(Path path,ProteinKmerGeometry geometry,String provenance) throws Exception{
		DDLPentamerBlacklist.loadReducedStructural(path,"unmeasured",2048,geometry,provenance);
		final LongHashSet set=new LongHashSet(4096);
		final ByteFile in=ByteFile.makeByteFile(path.toString(),true);
		try{
			for(byte[] line=in.nextLine();line!=null;line=in.nextLine()){
				if(line.length>0 && line[0]!='>'){
					final long key=geometry.word(DDLPentamerBlacklist.pack(line,6));
					require(set.add(key),"Duplicate compact blacklist loaded key: "+path);
				}
			}
		}finally{require(!in.close(),"Compact blacklist reader failed");}
		require(set.size()==2048,"Compact blacklist loaded key count mismatch: "+path);
		final long[] keys=set.toArray();
		Arrays.sort(keys);
		return keys;
	}

	private static void readCompactLengthMedians(Path path,FamilyReference[] records,
			HashMap<String,FamilyReference> byId,double[] medians){
		final ByteFile in=ByteFile.makeByteFile(path.toString(),false);
		final LineParser1 lp=new LineParser1('\t');final boolean[] seen=new boolean[records.length];
		int count=0;
		try{
			header(in,"family\tmembers\tlength_p10\tlength_median\tlength_p90\tcount_under_200");
			for(byte[] line=in.nextLine();line!=null;line=in.nextLine()){
				lp.set(line);require(lp.terms()==6,"Length medians require six fields");
				final FamilyReference f=byId.get(lp.parseString(0));
				require(f!=null && number(lp,1)==f.members,"Length median family/count mismatch");
				if(f.members==0){
					require(lp.termEquals("NA",2) && lp.termEquals("NA",3) &&
						lp.termEquals("NA",4),"Empty family length quantiles must be NA");
					continue;
				}
				final double median=Double.parseDouble(lp.parseString(3));
				require(Double.isFinite(median) && median>=40 && median<=200 &&
					median*2==Math.rint(median*2) && !seen[f.rank],
					"Invalid or duplicate compact median length");
				medians[f.rank]=median;seen[f.rank]=true;count++;
			}
		}finally{require(!in.close(),"Compact length median read failed");}
		require(count==live(records,seen),"Compact length medians must cover every live family");
	}

	private static void readCompactHydrophobicities(Path path,FamilyReference[] records,double[] kdMeans){
		final ByteFile in=ByteFile.makeByteFile(path.toString(),false);
		final LineParser1 lp=new LineParser1('\t');final boolean[] seen=new boolean[records.length];
		int rows=0,count=0;
		try{
			header(in,"family_rank\tfamily_id\tN\thydrophobicity_mean");
			for(byte[] line=in.nextLine();line!=null;line=in.nextLine()){
				lp.set(line);require(lp.terms()==4,"Hydrophobicity rows require four fields");
				final int rank=integer(lp,0);require(rank==rows++ && rank<records.length,
					"Hydrophobicity rows must be dense physical ranks");
				final FamilyReference f=records[rank];
				require(lp.termEquals(f.id,1) && number(lp,2)==f.members,
					"Hydrophobicity row identity/count mismatch");
				if(f.members==0){require(lp.termEquals("NA",3),"Empty hydrophobicity must be NA");continue;}
				final double value=Double.parseDouble(lp.parseString(3));
				require(Double.isFinite(value) && !seen[rank],"Invalid or duplicate hydrophobicity mean");
				kdMeans[rank]=value;seen[rank]=true;count++;
			}
		}finally{require(!in.close(),"Compact hydrophobicity read failed");}
		require(rows==records.length && count==live(records,seen),
			"Hydrophobicities must cover the full roster and every live family");
	}

	private static void readCompactGc(Path report,Path map,FamilyReference[] records,int[] stable,
			long[] measured,double[] means,double[] sd,CompactCensus census){
		final ByteFile rf=ByteFile.makeByteFile(report.toString(),false),mf=ByteFile.makeByteFile(map.toString(),false);
		final LineParser1 r=new LineParser1('\t'),m=new LineParser1('\t');
		int rows=0;long members=0;
		try{
			header(rf,FamilyGcInputs.TAIL_HEADER_TEXT);
			header(mf,"stable_family_id\troster_rank\tfamily_id\tassigned\tresidues");
			for(byte[] row=rf.nextLine();row!=null;row=rf.nextLine()){
				final byte[] mapped=mf.nextLine();require(mapped!=null,"Compact GC map ended early");
				r.set(row);m.set(mapped);require(r.terms()==10 && m.terms()==5,"Compact GC rows require ten/five fields");
				final int rank=integer(r,1);require(rank>=0 && rank<records.length,"Compact GC rank out of range");
				final FamilyReference f=records[rank];
				require(f.members>0 && stable[rank]<0,"Compact GC rows must name each live family once");
				final int stableId=integer(r,0);
				require(integer(m,0)==stableId && integer(m,1)==rank && r.termEquals(f.id,2) &&
					m.termEquals(f.id,2) && number(r,3)==f.members && number(m,3)==f.members &&
					number(r,4)==f.residues && number(m,4)==f.residues && number(r,5)==f.members &&
					integer(r,8)==census.minAa && integer(r,9)==census.maxAa,
					"Compact GC report/map identity, counts or interval mismatch");
				stable[rank]=stableId;measured[rank]=f.members;
				means[rank]=fraction(r,6);sd[rank]=fraction(r,7);
				require(f.members>1 || sd[rank]==0,"Singleton GC sample SD must be zero");
				members=Math.addExact(members,f.members);rows++;
			}
			require(mf.nextLine()==null,"Compact GC map contains rows absent from report");
		}finally{
			final boolean reportFailed=rf.close(),mapFailed=mf.close();
			require(!reportFailed && !mapFailed,"Compact GC report/map read failed");
		}
		require(rows==census.live && members==census.members,"Compact GC report must cover every live family");
	}

	private static double fraction(LineParser1 lp,int field){
		final double value=Double.parseDouble(lp.parseString(field));
		require(Double.isFinite(value) && value>=0 && value<=1,"Compact fraction out of range");
		return value;
	}

	private static int live(FamilyReference[] records,boolean[] seen){
		int count=0;
		for(FamilyReference f:records){
			if(f.members>0){require(seen[f.rank],"Missing compact live row");count++;}
		}
		return count;
	}

	private static FamilyCompactCodec.AlignmentCalibration compactAlignment(FamilyReference.AlignmentCalibration value){
		if(value==FamilyReference.AlignmentCalibration.ENABLED){return FamilyCompactCodec.AlignmentCalibration.ENABLED;}
		if(value==FamilyReference.AlignmentCalibration.UNAVAILABLE){return FamilyCompactCodec.AlignmentCalibration.UNAVAILABLE;}
		throw new IllegalArgumentException("Unknown reference alignment calibration");
	}

	/** Family-level nearest-rank distributions distinguish raw calibrated bounds from encoded final acceptance gates. */
	private static void writeGateQuantiles(Path out,Entry[] entries,FamilyReference[] records,int live) throws Exception{
		final String[] names={"identity_lower_raw","blosum_lower_raw","identity_upper_union_raw","blosum_upper_union_raw","identity_accept_f32","blosum_final_f32"};
		final double[][] values=new double[names.length][live];int slot=0;
		for(FamilyReference f:records){
			if(f.members==0){continue;}
			final FamilyGateExport.Row g=entries[f.rank].gates;
			values[0][slot]=g.lower[1];values[1][slot]=g.lower[2];values[2][slot]=g.upper[1];values[3][slot]=g.upper[2];
			values[4][slot]=f.identityUpper;values[5][slot++]=f.blosumFinal;
		}
		require(slot==live,"Quantile population must include every nonempty family exactly once");
		final StringBuilder table=new StringBuilder("metric\tn\tmin\tp10\tp25\tmedian\tp75\tp90\tmax\n");
		for(int i=0;i<names.length;i++){
			Arrays.sort(values[i]);table.append(names[i]).append('\t').append(live);
			for(int p:new int[]{0,10,25,50,75,90,100}){table.append('\t').append(live==0 ? "NA" : Double.toString(FamilyGateStats.quantile(values[i],p)));}
			table.append('\n');
		}
		write(out.resolve("gate_quantiles.tsv"),table.toString());
	}

	/** Stores the exact ranking projection once; no gate, scalar or biological score is recomputed. */
	private static void monomer(HashMap<String,String> options,int limit) throws Exception{
		for(String key:options.keySet()){
			require(Arrays.asList("mode","parent","parentsha80","out","maxfamilies","threads").contains(key),
				"Monomer conversion cannot change gates, scalars or calibration: "+key);
		}
		final int threads=options.containsKey("threads") ? Parse.parseIntKMG(options.get("threads")) : 1;
		require(threads>0 && threads<=32,"Reference conversion loading requires 1..32 threads");
		final Pins pins=new Pins();
		final Path parent=input(options,"parent",pins),out=absolute(required(options,"out"));
		require(!Files.exists(out) && !Files.isSymbolicLink(out) && out.getParent()!=null &&
			Files.isDirectory(out.getParent()) && out.getParent().toRealPath().equals(out.getParent()),
			"Require fresh monomer output directory with real parent");
		final FamilyReferenceLoader.Loaded loaded=FamilyReferenceLoader.load(parent.toString(),threads,limit);
		require(loaded.header.spectrumWidth==FamilyReferenceCodec.SPECTRUM_WIDTH,
			"Monomer conversion requires a 400-coordinate dimer parent, not an already projected package");
		final FamilyReferenceCodec.Header h=loaded.header;
		final FamilyReferenceCodec.Header compactHeader=new FamilyReferenceCodec.Header(FamilyReferenceCodec.MONOMER_SCHEMA,
			h.k,h.width,h.seed,h.totalFamilies,h.liveFamilies,h.emptyFamilies,h.formulaVersion,h.spectrumEncoding,h.scalarZ);
		pins.verify();Files.createDirectory(out);
		final ByteStreamWriter writer=new ByteStreamWriter(out.resolve("references.tsv").toString(),false,false,false);writer.start();
		final ByteBuilder buffer=new ByteBuilder();
		try{
			writer.print(FamilyReferenceCodec.appendHeader(compactHeader,buffer));
			for(FamilyReference f:loaded.all){
				final FamilyReference compact=projectMonomers(f);
				writer.print(FamilyReferenceCodec.appendFiveLines(compact,buffer.clear()));
			}
		}finally{require(!writer.poisonAndWait(),"Monomer reference writer failed; no PASS");}
		write(out.resolve("summary.tsv"),"families\tlive\tempty\tspectrum_values_before\tspectrum_values_after\n"+
			h.totalFamilies+"\t"+h.liveFamilies+"\t"+h.emptyFamilies+"\t"+((long)h.liveFamilies*400)+"\t"+((long)h.liveFamilies*20)+"\n");
		write(out.resolve("recipe.txt"),"mode=monomer\nparent="+parent+"\nparentsha80="+required(options,"parentsha80")+
			"\nschema="+FamilyReferenceCodec.MONOMER_SCHEMA+"\nprojection=FamilyTailCascade.dimerRowSums;ascending_float32_sum_per_row\n"+
			"normalization=none\nunchanged=roster,gates,scalar_bounds,consensus,DDL,counts\n"+
			"supported_search=combined25_ranking_without_dimer_diagnostics\nbiological_rescoring=false\n");
		write(out.resolve("inputs.tsv"),"path\tsha80\n"+parent+"\t"+required(options,"parentsha80")+"\n");
		pins.verify();write(out.resolve("PASS"),MONOMER_PASS+"\n");
		final ByteBuilder manifest=new ByteBuilder("file\tsha80\n");
		for(String name:new String[]{"references.tsv","summary.tsv","recipe.txt","inputs.tsv","PASS"}){
			manifest.append(name).tab().append(DigestSuffix.file(out.resolve(name).toString())).nl();
		}
		write(out.resolve("result.sha80.tsv"),manifest.toString());
		System.out.println(MONOMER_PASS+" families="+h.totalFamilies+" live="+h.liveFamilies);
	}

	/** Densifies a pinned survivor-only package and filters the matching DDL rows without recomputing registers. */
	private static void subset(HashMap<String,String> options,int limit) throws Exception{
		for(String key:options.keySet()){
			require(Arrays.asList("mode","parent","parentsha80","survivors","survivorssha80","ddlprofiles","ddlprofilessha80","out","maxfamilies","threads").contains(key),
				"Subset export cannot change gates, scalars, calibration or profile payloads: "+key);
		}
		final int threads=options.containsKey("threads") ? Parse.parseIntKMG(options.get("threads")) : 1;
		require(threads>0 && threads<=32,"Reference subset loading requires 1..32 threads");
		final Pins pins=new Pins();
		final Path parent=input(options,"parent",pins),survivors=input(options,"survivors",pins),
			ddlProfiles=input(options,"ddlprofiles",pins),out=absolute(required(options,"out"));
		require(!Files.exists(out) && !Files.isSymbolicLink(out) && out.getParent()!=null &&
			Files.isDirectory(out.getParent()) && out.getParent().toRealPath().equals(out.getParent()),
			"Require fresh subset output directory with real parent");
		final FamilyReferenceLoader.Loaded loaded=FamilyReferenceLoader.load(parent.toString(),threads,limit);
		final HashSet<String> keep=survivors(survivors);
		final ArrayList<FamilyReference> survivorsOut=new ArrayList<FamilyReference>();
		final ByteBuilder roster=new ByteBuilder("old_rank\tfamily_id\tnew_rank\tmembers\n");
		for(FamilyReference f:loaded.all){
			if(f.members==0 || !keep.remove(f.id)){continue;}
			roster.append(f.rank).tab().append(f.id).tab().append(survivorsOut.size()).tab().append(f.members).nl();
			survivorsOut.add(f);
		}
		require(survivorsOut.size()>0 && keep.isEmpty(),"Survivor list contains families that are absent or empty in the parent");
		final FamilyReference[] records=survivorsOut.toArray(new FamilyReference[survivorsOut.size()]);
		pins.verify();Files.createDirectory(out);
		final Path profilesOut=out.resolve("profiles.tsv.gz"),summaryOut=out.resolve("summary.tsv");
		final int profiled=filterProfiles(ddlProfiles,records,profilesOut,summaryOut,pins);
		require(profiled==records.length,"Subset profiles must cover every survivor exactly once");
		writeSubsetReferences(parent,out.resolve("references.tsv"),records);
		final FamilyReferenceLoader.Loaded subsetLoaded=FamilyReferenceLoader.load(out.resolve("references.tsv").toString(),threads,records.length);
		require(subsetLoaded.live.length==records.length && subsetLoaded.header.emptyFamilies==0,
			"Subset reference failed dense live-only validation");
		final FamilyReferenceCodec.Header h=loaded.header;
		write(out.resolve("ddlprofiles.tsv"),profilesOut+"\t"+DigestSuffix.file(profilesOut.toString())+"\t"+DigestSuffix.file(summaryOut.toString())+"\n");
		write(out.resolve("survivor_roster.tsv"),roster.toString());
		write(out.resolve("summary.txt"),"parent_families="+h.totalFamilies+"\nparent_live="+h.liveFamilies+
			"\nsurvivors="+records.length+"\nprofiles="+profiled+"\n");
		write(out.resolve("recipe.txt"),"mode=subset\nparent="+parent+"\nsurvivors="+survivors+
			"\nddlprofiles="+ddlProfiles+"\nrank_policy=parent_live_order_with_empty_and_retired_families_removed\n"+
			"unchanged=gates,scalar_bounds,consensus,DDL,spectrum,profile_rows\nbiological_rescoring=false\n");
		final ByteBuilder inputs=new ByteBuilder("path\tsha80\n");
		for(java.util.Map.Entry<Path,String> p:pins.values.entrySet()){inputs.append(p.getKey().toString()).tab().append(p.getValue()).nl();}
		write(out.resolve("inputs.tsv"),inputs.toString());pins.verify();write(out.resolve("PASS"),SUBSET_PASS+"\n");
		final ByteBuilder manifest=new ByteBuilder("file\tsha80\n");
		for(String name:new String[]{"references.tsv","profiles.tsv.gz","summary.tsv","ddlprofiles.tsv","survivor_roster.tsv","summary.txt","recipe.txt","inputs.tsv","PASS"}){
			manifest.append(name).tab().append(DigestSuffix.file(out.resolve(name).toString())).nl();
		}
		write(out.resolve("result.sha80.tsv"),manifest.toString());
		System.out.println(SUBSET_PASS+" families="+records.length);
	}

	/** Raises final acceptance gates on a pinned parent; all other reference fields are borrowed unchanged. */
	private static void raise(HashMap<String,String> options,int limit) throws Exception{
		for(String key:options.keySet()){
			require(Arrays.asList("mode","parent","parentsha80","overrides","overridessha80","overrideslist","overrideslistsha80","identityfloor","blosumfloor","out","maxfamilies","threads").contains(key),
				"Raw calibration/scalar options cannot be mixed with parent-reference overrides: "+key);
		}
		final int threads=options.containsKey("threads") ? Parse.parseIntKMG(options.get("threads")) : 1;
		require(threads>0 && threads<=32,"Reference override loading requires 1..32 threads");
		final Pins pins=new Pins();
		final boolean single=options.containsKey("overrides"),listed=options.containsKey("overrideslist"),
			floor=options.containsKey("identityfloor") || options.containsKey("blosumfloor");
		require(single==options.containsKey("overridessha80") && listed==options.containsKey("overrideslistsha80") &&
			!(single && listed) && (single || listed || floor),"Supply one pinned override table/list or explicit acceptance floors");
		final Path parent=input(options,"parent",pins),overrides=single ? input(options,"overrides",pins) : null,
			list=listed ? input(options,"overrideslist",pins) : null,out=absolute(required(options,"out"));
		final double identityFloor=options.containsKey("identityfloor") ? Double.parseDouble(options.get("identityfloor")) : 0,
			blosumFloor=options.containsKey("blosumfloor") ? Double.parseDouble(options.get("blosumfloor")) : 0;
		require(Double.isFinite(identityFloor) && identityFloor>=0 && identityFloor<=100 && Double.isFinite(blosumFloor) &&
			blosumFloor>=0 && blosumFloor<=Float.MAX_VALUE,"Acceptance floors must be finite, identity in0..100 and BLOSUM nonnegative");
		require(!Files.exists(out) && !Files.isSymbolicLink(out) && out.getParent()!=null &&
			Files.isDirectory(out.getParent()) && out.getParent().toRealPath().equals(out.getParent()),"Require fresh override output directory with real parent");
		final FamilyReferenceLoader.Loaded loaded=FamilyReferenceLoader.load(parent.toString(),threads,limit);
		final FamilyReference[] records=loaded.all.clone();
		final OverrideSet updates=new OverrideSet(records);
		if(single){updates.merge(overrides);}
		if(listed){
			final ByteFile in=ByteFile.makeByteFile(list.toString(),false);final LineParser1 lp=new LineParser1('\t');
			final HashSet<Path> paths=new HashSet<Path>();
			try{
				header(in,"path\tsha80");
				for(byte[] line=in.nextLine();line!=null;line=in.nextLine()){
					lp.set(line);require(lp.terms()==2,"Override list requires path and sha80");
					final Path table=absolute(lp.parseString(0));require(paths.add(table),"Duplicate override table path");
					pins.add(table,lp.parseString(1));updates.merge(table);
				}
			}finally{require(!in.close(),"Override list read failed");}
			require(!paths.isEmpty(),"Override list must contain at least one table");
		}
		if(floor){updates.floors(identityFloor,blosumFloor);}
		final int[] counts=updates.apply();
		final double[] maxima=new double[4];
		for(int i=0;i<records.length;i++){
			final FamilyReference old=loaded.all[i],now=records[i];
			if(old.members==0){continue;}
			maxima[0]=Math.max(maxima[0],old.identityUpper);maxima[1]=Math.max(maxima[1],now.identityUpper);
			maxima[2]=Math.max(maxima[2],old.blosumFinal);maxima[3]=Math.max(maxima[3],now.blosumFinal);
		}
		pins.verify();Files.createDirectory(out);
		final ByteStreamWriter writer=new ByteStreamWriter(out.resolve("references.tsv").toString(),false,false,false);writer.start();
		final ByteBuilder buffer=new ByteBuilder();
		try{
			writer.print(FamilyReferenceCodec.appendHeader(loaded.header,buffer));
			for(FamilyReference record:records){writer.print(FamilyReferenceCodec.appendFiveLines(record,buffer.clear()));}
		}finally{require(!writer.poisonAndWait(),"Raised reference writer failed; no PASS");}
		write(out.resolve("summary.tsv"),"rows_read\tidentity_changed\tblosum_changed\teither_changed\tunchanged_rows\tidentity_max_before\tidentity_max_after\tblosum_max_before\tblosum_max_after\n"+
			counts[0]+"\t"+counts[1]+"\t"+counts[2]+"\t"+counts[3]+"\t"+(counts[0]-counts[3])+"\t"+
			maxima[0]+"\t"+maxima[1]+"\t"+maxima[2]+"\t"+maxima[3]+"\n");
		write(out.resolve("recipe.txt"),"mode=raise\nparent="+parent+"\noverrides="+overrides+
			"\noverrideslist="+list+"\noverride_tables="+updates.tables+"\nsource_rows="+updates.sourceRows+
			"\nrows_read=unique_family_rows_in_tables_or_floor_application\n"+
			"acceptance=explicit_final_thresholds_not_raw_calibration_maxima\nrounding=max_union_original_doubles_then_upward_float32_once;reject_source_decrease_vs_same_parent\n"+
			"unchanged=header,roster,lowers,cosine,consensus,DDL,spectrum,scalar_bounds\n"+
			"header_formula=original_base_formula;explicit_final_overrides_recorded_in_pinned_table\n"+
			"missing_row=retain_parent;unchanged_stage=exact_parent_float_as_double\n"+
			"identityfloor="+identityFloor+"\nblosumfloor="+blosumFloor+"\nglobal_floors=max_with_parent_and_all_overrides;absent_is_inert_zero\nbiological_rescoring=false\n");
		final StringBuilder merged=new StringBuilder(ACCEPT_OVERRIDE_HEADER).append('\n');
		final StringBuilder sources=new StringBuilder("roster_rank\tfamily_id\tidentity_source\tblosum_source\n");
		for(int i=0;i<records.length;i++){
			if(records[i].identityUpper!=loaded.all[i].identityUpper || records[i].blosumFinal!=loaded.all[i].blosumFinal){
				merged.append(i).append('\t').append(records[i].id).append('\t').append(updates.identity[i]).append('\t').append(updates.blosum[i]).append('\n');
				sources.append(i).append('\t').append(records[i].id).append('\t').append(updates.identitySource[i]).append('\t').append(updates.blosumSource[i]).append('\n');
			}
		}
		write(out.resolve("merged_overrides.tsv"),merged.toString());
		write(out.resolve("override_sources.tsv"),sources.toString());
		final ByteBuilder inputs=new ByteBuilder("path\tsha80\n");
		for(java.util.Map.Entry<Path,String> p:pins.values.entrySet()){inputs.append(p.getKey().toString()).tab().append(p.getValue()).nl();}
		write(out.resolve("inputs.tsv"),inputs.toString());pins.verify();write(out.resolve("PASS"),RAISE_PASS+"\n");
		final ByteBuilder manifest=new ByteBuilder("file\tsha80\n");
		for(String name:new String[]{"references.tsv","summary.tsv","recipe.txt","merged_overrides.tsv","override_sources.tsv","inputs.tsv","PASS"}){
			manifest.append(name).tab().append(DigestSuffix.file(out.resolve(name).toString())).nl();
		}
		write(out.resolve("result.sha80.tsv"),manifest.toString());System.out.println(RAISE_PASS+" rows="+counts[0]+" changed="+counts[3]);
	}

	/** Exports legal M1 ordinary-fitter starts from each stored member lower and the runtime floors.
	 * Member-start20/40 improved one caller seed by +39 overall/+350 40-99aa raw objective,
	 * within seed noise, and was not adopted; see results/member_start_fit_20261116.
	 */
	private static void memberStart(HashMap<String,String> options,int limit) throws Exception{
		for(String key:options.keySet()){
			require(Arrays.asList("mode","parent","parentsha80","identityfloor","blosumfloor","out","maxfamilies","threads").contains(key),
				"Member-start export accepts only a parent reference and acceptance floors: "+key);
		}
		final int threads=options.containsKey("threads") ? Parse.parseIntKMG(options.get("threads")) : 1;
		require(threads>0 && threads<=32,"Member-start loading requires 1..32 threads");
		final double identityFloor=options.containsKey("identityfloor") ? Double.parseDouble(options.get("identityfloor")) : 30,
			blosumFloor=options.containsKey("blosumfloor") ? Double.parseDouble(options.get("blosumfloor")) : 0;
		require(Double.isFinite(identityFloor) && identityFloor>=30 && identityFloor<=100 && Double.isFinite(blosumFloor) &&
			blosumFloor>=0 && blosumFloor<=Float.MAX_VALUE,"Member-start floors must be finite, identity in30..100 and BLOSUM nonnegative");
		final String parentPin=options.containsKey("parentsha80") ? options.get("parentsha80") : "unmeasured";
		final boolean structuralParent=parentPin.equals("unmeasured");
		final Pins pins=structuralParent ? null : new Pins();
		final StructuralInputs inputs=structuralParent ? new StructuralInputs() : null;
		final Path parent=structuralParent ? structuralInput(options,"parent",inputs) : input(options,"parent",pins),
			out=absolute(required(options,"out"));
		require(!Files.exists(out) && !Files.isSymbolicLink(out) && out.getParent()!=null &&
			Files.isDirectory(out.getParent()) && out.getParent().toRealPath().equals(out.getParent()),
			"Require fresh member-start output directory with real parent");
		final FamilyReferenceLoader.Loaded loaded=FamilyReferenceLoader.loadMetricBound(parent.toString(),threads,limit);
		final FamilyReference[] records=loaded.all.clone();
		final int[] counts=applyMemberStarts(records,identityFloor,blosumFloor);
		final double[] maxima=new double[4];
		for(int i=0;i<records.length;i++){
			final FamilyReference old=loaded.all[i],now=records[i];
			if(old.members==0 || old.alignmentCalibration!=FamilyReference.AlignmentCalibration.ENABLED){continue;}
			maxima[0]=Math.max(maxima[0],old.identityUpper);maxima[1]=Math.max(maxima[1],now.identityUpper);
			maxima[2]=Math.max(maxima[2],old.blosumFinal);maxima[3]=Math.max(maxima[3],now.blosumFinal);
		}
		if(structuralParent){inputs.verify();}else{pins.verify();}
		Files.createDirectory(out);
		final ByteStreamWriter writer=new ByteStreamWriter(out.resolve("references.tsv").toString(),false,false,false);writer.start();
		final ByteBuilder buffer=new ByteBuilder();
		try{
			writer.print(FamilyReferenceCodec.appendHeader(loaded.header,buffer));
			for(FamilyReference record:records){writer.print(FamilyReferenceCodec.appendFiveLines(record,buffer.clear(),loaded.header));}
		}finally{require(!writer.poisonAndWait(),"Member-start reference writer failed; no PASS");}
		validateMemberStartOutput(out.resolve("references.tsv"),loaded.header,records,threads,records.length);
		write(out.resolve("summary.tsv"),"rows_read\tidentity_changed\tblosum_changed\teither_changed\tunchanged_rows\tidentity_floor_bound\tblosum_floor_bound\talignment_unavailable\tidentity_max_before\tidentity_max_after\tblosum_max_before\tblosum_max_after\n"+
			counts[0]+"\t"+counts[1]+"\t"+counts[2]+"\t"+counts[3]+"\t"+(counts[0]-counts[3])+"\t"+
			counts[4]+"\t"+counts[5]+"\t"+counts[6]+"\t"+
			maxima[0]+"\t"+maxima[1]+"\t"+maxima[2]+"\t"+maxima[3]+"\n");
		write(out.resolve("recipe.txt"),"mode=memberstart\nparent="+parent+"\nparentsha80="+parentPin+
			"\nrows_read=live_alignment_enabled_families\nacceptance=max(runtime_floor,stored_member_lower)\n"+
			"member_multiplier=1.0\nidentityfloor="+identityFloor+"\nblosumfloor="+blosumFloor+"\n"+
			"cosine=retain_parent_final\nalignment_unavailable=retain_parent_inert_gates\n"+
			"unchanged=header,roster,lowers,cosine,consensus,DDL,spectrum,scalar_bounds,alignment_state\n"+
			"biological_rescoring=false\n");
		write(out.resolve("inputs.tsv"),structuralParent ? structuralInputRecord(parent) :
			"path\tsha80\n"+parent+"\t"+parentPin+"\n");
		if(structuralParent){inputs.verify();}else{pins.verify();}
		write(out.resolve("PASS"),MEMBER_START_PASS+"\n");
		System.out.println(MEMBER_START_PASS+" rows="+counts[0]+" changed="+counts[3]);
	}

	private static int[] applyMemberStarts(FamilyReference[] records,double identityFloor,double blosumFloor){
		assert(records!=null && Double.isFinite(identityFloor) && Double.isFinite(blosumFloor)) : "Member-start CLI validates records and floors";
		final int[] counts=new int[7];
		for(int i=0;i<records.length;i++){
			final FamilyReference f=records[i];
			if(f.members==0){continue;}
			if(f.alignmentCalibration!=FamilyReference.AlignmentCalibration.ENABLED){counts[6]++;continue;}
			counts[0]++;
			final float id=FamilyGateCodec.roundUpper(Math.max(identityFloor,f.identityLower));
			final float bl=FamilyGateCodec.roundUpper(Math.max(blosumFloor,f.blosumLower));
			final boolean ic=id!=f.identityUpper,bc=bl!=f.blosumFinal;
			if(id==FamilyGateCodec.roundUpper(identityFloor)){counts[4]++;}
			if(bl==FamilyGateCodec.roundUpper(blosumFloor)){counts[5]++;}
			if(ic){counts[1]++;}if(bc){counts[2]++;}if(!ic && !bc){continue;}counts[3]++;
			records[i]=copyAcceptGates(f,id,bl);
		}
		assert(counts[3]<=counts[0] && counts[1]<=counts[3] && counts[2]<=counts[3]) : "Every member-start change belongs to one enabled family";
		return counts;
	}

	private static FamilyReference copyAcceptGates(FamilyReference f,float identityUpper,float blosumFinal){
		if(f.explicitAlignmentCalibration){
			return new FamilyReference(f.id,f.rank,f.members,f.residues,f.validWords,f.consensus,f.ddl,f.spectrum,
				f.cosineLower,f.identityLower,f.blosumLower,f.cosineAccept,identityUpper,blosumFinal,
				f.consensusRaw,f.scalarMask,f.scalarBounds,f.alignmentCalibration);
		}
		return new FamilyReference(f.id,f.rank,f.members,f.residues,f.validWords,f.consensus,f.ddl,f.spectrum,
			f.cosineLower,f.identityLower,f.blosumLower,f.cosineAccept,identityUpper,blosumFinal,
			f.consensusRaw,f.scalarMask,f.scalarBounds);
	}

	private static void validateMemberStartOutput(Path path,FamilyReferenceCodec.Header header,
			FamilyReference[] records,int threads,int limit){
		final FamilyReferenceLoader.Loaded loaded=FamilyReferenceLoader.loadMetricBound(path.toString(),threads,limit);
		require(sameHeader(header,loaded.header) && loaded.all.length==records.length,
			"Member-start output header/census mismatch");
		for(int i=0;i<records.length;i++){
			require(sameReference(records[i],loaded.all[i]),"Member-start output row mismatch: "+i);
		}
	}

	/** Replays one saved raw_calibration.tsv family roster through a new identity-only margin.
	 * results/identity_padding_fuller_replays_20261116 found fuller +1 gained
	 * 7266..28218 positives in the six paired cells, but4/100 negatives rose
	 * SHORT3->10 and LONG0->1, so no default changed.
	 */
	private static void identityPadding(HashMap<String,String> options,int limit) throws Exception{
		for(String key:options.keySet()){
			require(Arrays.asList("mode","parent","rawcalibration","identityuppermultiplier","identityupperadditive","out","maxfamilies","threads").contains(key),
				"Identity padding replay accepts only a parent reference, saved raw calibration and identity margin: "+key);
		}
		final int threads=options.containsKey("threads") ? Parse.parseIntKMG(options.get("threads")) : 1;
		require(threads>0 && threads<=32,"Identity padding loading requires 1..32 threads");
		final double identityUpperMultiplier=identityPadding(options,"identityuppermultiplier",
			FamilyGateCodec.DEFAULT_IDENTITY_UPPER_MULTIPLIER,1);
		final double identityUpperAdditive=identityPadding(options,"identityupperadditive",
			FamilyGateCodec.DEFAULT_IDENTITY_UPPER_ADDITIVE,0);
		final StructuralInputs inputs=new StructuralInputs();
		final Path parent=structuralInput(options,"parent",inputs),raw=structuralInput(options,"rawcalibration",inputs),
			out=absolute(required(options,"out"));
		require(!Files.exists(out) && !Files.isSymbolicLink(out) && out.getParent()!=null &&
			Files.isDirectory(out.getParent()) && out.getParent().toRealPath().equals(out.getParent()),
			"Require fresh identity-padding output directory with real parent");
		final FamilyReferenceLoader.Loaded loaded=FamilyReferenceLoader.loadMetricBound(parent.toString(),threads,limit);
		final FamilyReference[] records=loaded.all.clone();
		final int[] counts=applyRawIdentityPadding(raw,loaded.header,records,identityUpperMultiplier,identityUpperAdditive);
		Files.createDirectory(out);
		final ByteStreamWriter writer=new ByteStreamWriter(out.resolve("references.tsv").toString(),false,false,false);writer.start();
		final ByteBuilder buffer=new ByteBuilder();
		try{
			writer.print(FamilyReferenceCodec.appendHeader(loaded.header,buffer));
			for(FamilyReference record:records){writer.print(FamilyReferenceCodec.appendFiveLines(record,buffer.clear(),loaded.header));}
		}finally{require(!writer.poisonAndWait(),"Identity-padded reference writer failed; no PASS");}
		write(out.resolve("summary.tsv"),"rows_read\tlive_rows\talignment_enabled\tunavailable_retained\tidentity_changed\tidentity_unchanged\n"+
			counts[0]+"\t"+counts[1]+"\t"+counts[2]+"\t"+counts[3]+"\t"+counts[4]+"\t"+(counts[2]-counts[4])+"\n");
		write(out.resolve("recipe.txt"),"mode=identitypadding\nparent="+parent+"\nrawcalibration="+raw+
			"\nacceptance=saved_raw_identity_lower_and_upper_maxima\nidentity_formula=round_up_float32(max(lower,"+
			identityUpperMultiplier+"*upper,upper+"+identityUpperAdditive+",30));strict_greater\n"+
			"absent_upper=max(raw_lower,30);no_upper_margin_terms\nunavailable_alignment=preserve_inert_parent_alignment_fields\n"+
			"unchanged=header,roster,lowers,cosine,blosum,consensus,DDL,spectrum,scalar_bounds,alignment_state\nbiological_rescoring=false\n");
		validateIdentityPaddingOutput(out.resolve("references.tsv"),loaded.header,records,threads,limit);
		inputs.verify();
		write(out.resolve("PASS"),IDENTITY_PADDING_PASS+"\n");
		System.out.println(IDENTITY_PADDING_PASS+" rows="+counts[0]+" changed="+counts[4]);
	}

	/** Applies final encoded gates from a prior package as per-family floors for the same calibrated targets. */
	private static int[] applyHistoricalFloors(Path path,FamilyReference[] records,int threads,int limit){
		final FamilyReferenceLoader.Loaded floors=FamilyReferenceLoader.load(path.toString(),threads,limit);
		require(floors.all.length==records.length && floors.header.totalFamilies==records.length &&
			floors.header.liveFamilies==floors.live.length,"Historical floor reference census mismatch");
		final int[] counts=new int[3];
		for(int i=0;i<records.length;i++){
			final FamilyReference current=records[i],floor=floors.all[i];
			require(current!=null && floor!=null && current.rank==i && floor.rank==i && current.id.equals(floor.id),
				"Historical floor roster mismatch at rank "+i);
			require((current.members==0)==(floor.members==0),"Historical floor live/empty mismatch: "+current.id);
			if(current.members==0){continue;}
			require(floor.members>0 && Arrays.equals(current.consensusRaw,floor.consensusRaw),
				"Historical floor consensus mismatch: "+current.id);
			final float id=(float)Math.max(current.identityUpper,floor.identityUpper);
			final float bl=(float)Math.max(current.blosumFinal,floor.blosumFinal);
			if(id>current.identityUpper || bl>current.blosumFinal){
				records[i]=new FamilyReference(current.id,current.rank,current.members,current.residues,current.validWords,
					current.consensus,current.ddl,current.spectrum,current.cosineLower,current.identityLower,current.blosumLower,
					current.cosineAccept,id,bl,current.consensusRaw,current.scalarMask,current.scalarBounds);
				counts[0]++; if(id>current.identityUpper){counts[1]++;} if(bl>current.blosumFinal){counts[2]++;}
			}
		}
		return counts;
	}

	/** Returns rows, identity changes, BLOSUM changes and union changes. Only touched live rows are replaced. */
	public static int[] applyAcceptOverrides(Path path,FamilyReference[] records){
		require(path!=null && records!=null,"Acceptance overrides require a table and parent roster");
		final OverrideSet updates=new OverrideSet(records);updates.merge(path);return updates.apply();
	}

	/** Every source is checked against the same immutable parent before any union is installed. */
	private static final class OverrideSet{
		OverrideSet(FamilyReference[] records_){
			require(records_!=null,"Overrides require a parent roster");records=records_;seen=new boolean[records.length];
			identity=new double[records.length];blosum=new double[records.length];
			identitySource=new String[records.length];blosumSource=new String[records.length];
			Arrays.fill(identitySource,"parent");Arrays.fill(blosumSource,"parent");
			for(int i=0;i<records.length;i++){
				final FamilyReference f=records[i];require(f!=null && f.rank==i,"Override parent must have a complete ordered roster");
				require(!f.explicitAlignmentCalibration,"Legacy override tables do not bind the new alignment metric: "+f.id);
				identity[i]=f.identityUpper;blosum[i]=f.blosumFinal;
			}
		}
		void merge(Path path){
			final ByteFile in=ByteFile.makeByteFile(path.toString(),false);final LineParser1 lp=new LineParser1('\t');
			final boolean[] local=new boolean[records.length];tables++;
			final String source=path.toString();
			try{
				header(in,ACCEPT_OVERRIDE_HEADER);
				for(byte[] line=in.nextLine();line!=null;line=in.nextLine()){
					lp.set(line);require(lp.terms()==4,"Acceptance override requires four columns");
					final int rank=integer(lp,0);require(rank<records.length && !local[rank],"Unknown or duplicate per-table override rank: "+rank);
					final FamilyReference f=records[rank];
					require(f.members>0 && lp.termEquals(f.id,1),"Override must match exact live parent rank/family: "+rank);
					final double id=Double.parseDouble(lp.parseString(2)),bl=Double.parseDouble(lp.parseString(3));
					require(FamilyGateCodec.roundUpper(id)>=f.identityUpper && FamilyGateCodec.roundUpper(bl)>=f.blosumFinal,
						"Every override table must preserve its parent's encoded gates: "+f.id);
					if(id>identity[rank] || id==identity[rank] && id>f.identityUpper && source.compareTo(identitySource[rank])<0){identity[rank]=id;identitySource[rank]=source;}
					if(bl>blosum[rank] || bl==blosum[rank] && bl>f.blosumFinal && source.compareTo(blosumSource[rank])<0){blosum[rank]=bl;blosumSource[rank]=source;}
					local[rank]=seen[rank]=true;sourceRows=Math.addExact(sourceRows,1);
				}
			}finally{require(!in.close(),"Acceptance override reader failed");}
		}
		void floors(double id,double bl){
			assert(Double.isFinite(id) && Double.isFinite(bl) && id>=0 && bl>=0) : "CLI floor bounds were validated before parent mutation";
			for(int i=0;i<records.length;i++){
				if(records[i].members==0){continue;}seen[i]=true;
				if(id>identity[i]){identity[i]=id;identitySource[i]="identityfloor";}
				if(bl>blosum[i]){blosum[i]=bl;blosumSource[i]="blosumfloor";}
			}
		}
		int[] apply(){
			final int[] counts=new int[4];
			for(int rank=0;rank<records.length;rank++){
				if(!seen[rank]){continue;}counts[0]++;final FamilyReference f=records[rank];
				final float id=FamilyGateCodec.roundUpper(identity[rank]),bl=FamilyGateCodec.roundUpper(blosum[rank]);
				final boolean ic=id!=f.identityUpper,bc=bl!=f.blosumFinal;
				if(ic){counts[1]++;}if(bc){counts[2]++;}if(!ic && !bc){continue;}counts[3]++;
				records[rank]=new FamilyReference(f.id,f.rank,f.members,f.residues,f.validWords,
					f.consensus,f.ddl,f.spectrum,f.cosineLower,f.identityLower,f.blosumLower,f.cosineAccept,
					id,bl,f.consensusRaw,f.scalarMask,f.scalarBounds);
			}
			assert(counts[3]<=counts[0] && counts[1]<=counts[3] && counts[2]<=counts[3]) : "Every changed gate belongs to one considered family";
			return counts;
		}
		final FamilyReference[] records;final boolean[] seen;final double[] identity,blosum;
		final String[] identitySource,blosumSource;
		int tables;long sourceRows;
	}

	/** Joins every all-member scalar row to the exact roster; bounds never alter alignment gates or ranking payloads. */
	public static void applyScalarBounds(Path path,FamilyReference[] records,double z){
		require(path!=null && records!=null && Double.isFinite(z) && z>0,"Scalar export requires a table, roster and positive finite z");
		final StringBuilder expected=new StringBuilder("family_rank\tfamily_id\tN");
		for(String name:new String[]{"hydrophobicity","net_charge","aromaticity","fraction_small","pI_tail"}){
			for(String metric:new String[]{"mean","sd","min","max"}){expected.append('\t').append(name).append('_').append(metric);}
		}
		final ByteFile in=ByteFile.makeByteFile(path.toString(),false); final LineParser1 lp=new LineParser1('\t');
		final boolean[] seen=new boolean[records.length]; final double[] statistics=new double[20]; int count=0;
		try{
			header(in,expected.toString());
			for(byte[] line=in.nextLine();line!=null;line=in.nextLine()){
				lp.set(line); require(lp.terms()==23,"All-member scalar table requires 23 columns");
				final int rank=integer(lp,0); final long n=number(lp,2);
				require(rank<records.length && !seen[rank],"Unknown or duplicate scalar family rank: "+rank);
				final FamilyReference f=records[rank];
				require(f!=null && f.rank==rank && lp.termEquals(f.id,1) && f.members==n,"Scalar roster/member mismatch at rank "+rank);
				seen[rank]=true; count++;
				for(int i=0;i<5;i++){
					final int offset=3+4*i;
					if(n==0){
						for(int j=0;j<4;j++){require(lp.termEquals("NA",offset+j),"Zero-member scalar statistics must be NA");}
						continue;
					}
					for(int j=0;j<4;j++){statistics[4*i+j]=finite(lp,offset+j);}
				}
				records[rank]=withScalarBounds(f,n==0 ? null : statistics,z);
			}
		}finally{require(!in.close(),"Scalar table read failed");}
		require(count==records.length,"Scalar table must cover the complete family roster");
	}

	/** Stores the existing ascending float row sums without renormalizing or changing other fields. */
	public static FamilyReference projectMonomers(FamilyReference f){
		require(f!=null, "Monomer projection requires an explicit reference, including empty records");
		if(f.members==0){return f;}
		require(f.spectrum.length==FamilyReferenceCodec.SPECTRUM_WIDTH,
			"Monomer projection requires400 dimer coordinates, not an already projected record");
		final float[] projected=new float[FamilyReferenceCodec.MONOMER_WIDTH];
		FamilyTailCascade.dimerRowSums(f.spectrum, projected);
		return copyPayload(f, projected, f.scalarMask, f.scalarBounds);
	}

	/** Builds one live reference from measured original-double maxima, after the caller has united
	 * all required negative panels. Cosine must be the separate raw-dimer maximum, not combined25.
	 * Scalar z=0 explicitly omits scalar gates; positive z uses the prepared rounded statistics.
	 * Empty families use FamilyReference.empty instead. Consensus and plain modal arrays are borrowed
	 * immutable; frequencies and scalar bounds are newly owned. Positional rows remain a separate payload.
	 */
	public static FamilyReference fromPrepared(FamilyMemberPreparation.Prepared f, double cosineUpper,
			double identityUpper, double blosumUpper, double blosumFloor, double scalarZ){
		return fromPrepared(f,cosineUpper,identityUpper,blosumUpper,blosumFloor,scalarZ,
			FamilyGateCodec.DEFAULT_IDENTITY_UPPER_MULTIPLIER,FamilyGateCodec.DEFAULT_IDENTITY_UPPER_ADDITIVE);
	}

	/** Builds one legacy reference with explicit identity padding; padding is validated
	 * before any evidence-specific branch so direct callers fail closed on bad margins.
	 */
	public static FamilyReference fromPrepared(FamilyMemberPreparation.Prepared f, double cosineUpper,
			double identityUpper, double blosumUpper, double blosumFloor, double scalarZ,
			double identityUpperMultiplier, double identityUpperAdditive){
		requireIdentityPadding(identityUpperMultiplier,identityUpperAdditive);
		require(f!=null && f.members>0, "Final calibrated references require a live prepared family; empty families have no measured gates");
		require(f.alignmentMetric.same(FamilyAlignmentMetric.LEGACY),
			"Legacy export cannot combine explicit-metric member cutoffs with unbound negative uppers: "+f.id);
		require(Double.isFinite(blosumFloor) && blosumFloor>=0 && blosumFloor<=Float.MAX_VALUE,
			"BLOSUM floor must be finite, nonnegative and representable");
		require(Double.isFinite(scalarZ) && scalarZ>=0, "Scalar z must be finite and nonnegative; zero explicitly disables scalar gates");
		final double cosineLower=f.lower[0].cutoff, identityLower=f.lower[1].cutoff, blosumLower=f.lower[2].cutoff;
		final FamilyReference result=new FamilyReference(f.id, f.rank, f.members, f.residues, f.validWords,
			f.consensus, f.plain, frequencies(f.spectrum),
			FamilyGateCodec.lower(cosineLower), FamilyGateCodec.lower(identityLower), FamilyGateCodec.lower(blosumLower),
			FamilyGateCodec.cosineAccept(cosineLower, cosineUpper),
			FamilyGateCodec.upper(identityLower, identityUpper, identityUpperMultiplier, identityUpperAdditive),
			FamilyGateCodec.finalBlosum(Math.max(blosumLower, blosumFloor), blosumUpper), f.consensusRaw);
		return scalarZ==0 ? result : withScalarBounds(result, f.statistics, scalarZ);
	}

	/** Merges actual per-panel evidence for one explicitly bound family before applying margins once.
	 * An unavailable negative upper contributes no margin terms. Missing member alignment cutoffs
	 * disable both alignment channels, independently of cosine/scalar state. The enclosing driver
	 * verifies that this array contains every required panel and retains each panel's outcome counts.
	 */
	public static FamilyReference fromPrepared(FamilyMemberPreparation.Prepared f, double cosineUpper,
			FamilyNegativeUpperCalibration.MetricResult[] evidence, double blosumFloor, double scalarZ){
		return fromPrepared(f,cosineUpper,evidence,blosumFloor,scalarZ,
			FamilyGateCodec.DEFAULT_IDENTITY_UPPER_MULTIPLIER,FamilyGateCodec.DEFAULT_IDENTITY_UPPER_ADDITIVE);
	}

	/** Builds one metric-bound reference with explicit identity padding. Invalid padding
	 * is rejected even when all alignment channels are unavailable or no upper was measured.
	 */
	public static FamilyReference fromPrepared(FamilyMemberPreparation.Prepared f, double cosineUpper,
			FamilyNegativeUpperCalibration.MetricResult[] evidence, double blosumFloor, double scalarZ,
			double identityUpperMultiplier, double identityUpperAdditive){
		requireIdentityPadding(identityUpperMultiplier,identityUpperAdditive);
		require(f!=null && f.members>0 && !f.alignmentMetric.same(FamilyAlignmentMetric.LEGACY) &&
			f.measuredLower!=null && f.measuredLower.length==3,"Metric export requires measured-state member preparation");
		require(evidence!=null && evidence.length>0,"Every required negative panel must supply explicit evidence, including unavailable results");
		require(Double.isFinite(blosumFloor) && blosumFloor>=0 && blosumFloor<=Float.MAX_VALUE,
			"BLOSUM floor must be finite, nonnegative and representable");
		require(Double.isFinite(scalarZ) && scalarZ>=0,"Scalar z must be finite and nonnegative");
		boolean upperMeasured=false; double identityUpper=0, blosumUpper=0;
		for(FamilyNegativeUpperCalibration.MetricResult result:evidence){
			require(result!=null && result.familyRank==f.rank && result.familyId.equals(f.id) && result.metric.same(f.alignmentMetric),
				"Negative evidence must bind this physical family and its exact alignment metric: "+f.id);
			if(!result.measured()){continue;}
			if(!upperMeasured || result.identityUpper()>identityUpper){identityUpper=result.identityUpper();}
			if(!upperMeasured || result.blosumUpper()>blosumUpper){blosumUpper=result.blosumUpper();}
			upperMeasured=true;
		}
		final boolean enabled=f.measuredLower[1].measured();
		require(enabled==f.measuredLower[2].measured(),"Identity and BLOSUM member availability must describe the same qualified alignments");
		final double cosineLower=f.measuredLower[0].cutoff();
		float identityLower=0, blosumLower=0, identityAccept=0, blosumAccept=0;
		if(enabled){
			final double id=f.measuredLower[1].cutoff(), bl=f.measuredLower[2].cutoff();
			identityLower=FamilyGateCodec.lower(id); blosumLower=FamilyGateCodec.lower(bl);
			identityAccept=upperMeasured ? FamilyGateCodec.upper(id, identityUpper, identityUpperMultiplier, identityUpperAdditive) :
				FamilyGateCodec.roundUpper(Math.max(id,30));
			blosumAccept=upperMeasured ? FamilyGateCodec.finalBlosum(Math.max(bl,blosumFloor),blosumUpper) :
				FamilyGateCodec.roundUpper(Math.max(0,Math.max(bl,blosumFloor)));
		}
		final FamilyReference result=new FamilyReference(f.id,f.rank,f.members,f.residues,f.validWords,
			f.consensus,f.plain,frequencies(f.spectrum),FamilyGateCodec.lower(cosineLower),identityLower,blosumLower,
			FamilyGateCodec.cosineAccept(cosineLower,cosineUpper),identityAccept,blosumAccept,f.consensusRaw,0,new double[10],
			enabled ? FamilyReference.AlignmentCalibration.ENABLED : FamilyReference.AlignmentCalibration.UNAVAILABLE);
		return scalarZ==0 ? result : withScalarBounds(result,f.statistics,scalarZ);
	}

	/** Copies only scalar gates, borrowing (never retaining) five rounded mean/SD/min/max groups.
	 * Use FamilyMemberScalars.Moments.exportStatistics for preparation to preserve the existing decimal boundary.
	 * Empty references require null statistics and are returned unchanged. Live references require twenty finite values.
	 */
	public static FamilyReference withScalarBounds(FamilyReference f,double[] statistics,double z){
		require(f!=null && Double.isFinite(z) && z>0,"Scalar bounds require a reference and positive finite z");
		if(f.members==0){require(statistics==null,"Empty families have no measured scalar statistics");return f;}
		require(statistics!=null && statistics.length==20,"Live scalar bounds require five mean/SD/min/max groups");
		final double[] bounds=new double[10]; int mask=0;
		for(int i=0;i<5;i++){
			final int offset=4*i;
			final double mean=statistics[offset],sd=statistics[offset+1],min=statistics[offset+2],max=statistics[offset+3];
			require(Double.isFinite(mean) && Double.isFinite(sd) && Double.isFinite(min) && Double.isFinite(max),
				"All-member scalar statistics must be finite for N>0");
			require(sd>=0 && min<=mean && mean<=max && (f.members!=1 || sd==0),"Invalid scalar population moments: "+f.id);
			if(f.members<30 || sd==0){continue;}
			final double lo=mean-z*sd,hi=mean+z*sd;
			require(Double.isFinite(lo) && Double.isFinite(hi) && lo<hi,"Scalar bounds overflow or collapse: "+f.id);
			mask|=1<<i; bounds[2*i]=lo; bounds[2*i+1]=hi;
		}
		return copyPayload(f,f.spectrum,mask,bounds);
	}
	/** Shape/scalar copies retain explicit framing, including unavailable alignment channels. */
	private static FamilyReference copyPayload(FamilyReference f,float[] spectrum,int scalarMask,double[] scalarBounds){
		if(f.explicitAlignmentCalibration){
			return new FamilyReference(f.id,f.rank,f.members,f.residues,f.validWords,f.consensus,f.ddl,spectrum,
				f.cosineLower,f.identityLower,f.blosumLower,f.cosineAccept,f.identityUpper,f.blosumFinal,
				f.consensusRaw,scalarMask,scalarBounds,f.alignmentCalibration);
		}
		return new FamilyReference(f.id,f.rank,f.members,f.residues,f.validWords,f.consensus,f.ddl,spectrum,
			f.cosineLower,f.identityLower,f.blosumLower,f.cosineAccept,f.identityUpper,f.blosumFinal,
			f.consensusRaw,scalarMask,scalarBounds);
	}

	private static FamilyReference copyIdentityUpper(FamilyReference f,float identityUpper){
		if(f.explicitAlignmentCalibration){
			return new FamilyReference(f.id,f.rank,f.members,f.residues,f.validWords,f.consensus,f.ddl,f.spectrum,
				f.cosineLower,f.identityLower,f.blosumLower,f.cosineAccept,identityUpper,f.blosumFinal,
				f.consensusRaw,f.scalarMask,f.scalarBounds,f.alignmentCalibration);
		}
		return new FamilyReference(f.id,f.rank,f.members,f.residues,f.validWords,f.consensus,f.ddl,f.spectrum,
			f.cosineLower,f.identityLower,f.blosumLower,f.cosineAccept,identityUpper,f.blosumFinal,
			f.consensusRaw,f.scalarMask,f.scalarBounds);
	}

	private static int[] applyRawIdentityPadding(Path raw,FamilyReferenceCodec.Header header,FamilyReference[] records,
			double multiplier,double additive){
		final ByteFile in=ByteFile.makeByteFile(raw.toString(),false);
		final LineParser1 lp=new LineParser1('\t');final boolean[] seen=new boolean[records.length];
		final int[] counts=new int[5];
		final boolean metric=!header.alignmentMetric.same(FamilyAlignmentMetric.LEGACY);
		try{
			final byte[] rawHeader=in.nextLine();require(rawHeader!=null,"Raw calibration is empty");
			lp.set(rawHeader);require(lp.terms()>=(metric ? 18 : 12),"Raw calibration is missing saved identity or metric columns");
			final String[] prefix={"roster_rank","family_id","members","cosine_lower","cosine_retained","identity_lower",
				"identity_retained","blosum_lower","blosum_retained","cosine_upper","identity_upper","blosum_upper"};
			for(int i=0;i<prefix.length;i++){require(lp.termEquals(prefix[i],i),"Unexpected raw calibration column "+(i+1));}
			if(metric){
				final String[] bound={"alignment_mode","alignment_min_span","member_alignment_measured","member_alignment_failed",
					"identity_quota_achievable","blosum_quota_achievable"};
				for(int i=0;i<bound.length;i++){require(lp.termEquals(bound[i],12+i),"Unexpected metric raw calibration column "+(13+i));}
			}
			for(byte[] line=in.nextLine();line!=null;line=in.nextLine()){
				lp.set(line);require(lp.terms()>=(metric ? 18 : 12),"Raw calibration row is missing identity or metric columns");
				final int rank=integer(lp,0);require(rank<records.length && !seen[rank],"Unknown/duplicate raw calibration rank: "+rank);
				final FamilyReference f=records[rank];final long members=number(lp,2);
				require(f.rank==rank && lp.termEquals(f.id,1) && f.members==members,"Raw calibration roster mismatch at rank "+rank);
				seen[rank]=true;counts[0]++;
				if(f.members==0){
					require(rawNa(lp,5) && rawNa(lp,10),"Empty raw calibration rows must retain NA alignment gates: "+f.id);
					continue;
				}
				counts[1]++;
				if(metric){
					require(lp.termEquals(header.alignmentMetric.name(),12) && integer(lp,13)==header.alignmentMetric.minimumSpan,
						"Raw calibration alignment metric disagrees with parent header: "+f.id);
				}
				if(f.alignmentCalibration==FamilyReference.AlignmentCalibration.UNAVAILABLE){
					require(rawNa(lp,5),"Unavailable member alignments must retain NA raw identity lower: "+f.id);
					if(!rawNa(lp,10)){rawIdentity(lp,10);}
					counts[3]++;continue;
				}
				counts[2]++;
				final double lower=rawIdentity(lp,5);final boolean hasUpper=!rawNa(lp,10);
				final double upper=hasUpper ? rawIdentity(lp,10) : 0;
				require(FamilyGateCodec.lower(lower)==f.identityLower,"Raw identity lower no longer matches parent row: "+f.id);
				final float identity=!hasUpper ? FamilyGateCodec.roundUpper(Math.max(lower,30)) :
					FamilyGateCodec.upper(lower,upper,multiplier,additive);
				if(identity!=f.identityUpper){counts[4]++;}
				records[rank]=copyIdentityUpper(f,identity);
			}
		}finally{require(!in.close(),"Raw calibration reader failed");}
		for(int i=0;i<seen.length;i++){require(seen[i],"Raw calibration is missing rank "+i);}
		return counts;
	}

	private static void validateIdentityPaddingOutput(Path path,FamilyReferenceCodec.Header header,
			FamilyReference[] records,int threads,int limit){
		final FamilyReferenceLoader.Loaded loaded=FamilyReferenceLoader.loadMetricBound(path.toString(),threads,limit);
		require(sameHeader(header,loaded.header) && loaded.all.length==records.length,
			"Identity padding output header/census mismatch");
		for(int i=0;i<records.length;i++){
			require(sameReference(records[i],loaded.all[i]),"Identity padding output row mismatch: "+i);
		}
	}

	private static boolean sameHeader(FamilyReferenceCodec.Header a,FamilyReferenceCodec.Header b){
		return a.schema.equals(b.schema) && a.k==b.k && a.width==b.width && a.seed==b.seed &&
			a.totalFamilies==b.totalFamilies && a.liveFamilies==b.liveFamilies && a.emptyFamilies==b.emptyFamilies &&
			a.formulaVersion.equals(b.formulaVersion) && a.spectrumEncoding.equals(b.spectrumEncoding) &&
			a.scalarZ==b.scalarZ && a.alignmentMetric.same(b.alignmentMetric);
	}

	private static boolean sameReference(FamilyReference a,FamilyReference b){
		if(a==b){return true;}
		return a.id.equals(b.id) && a.rank==b.rank && a.members==b.members && a.residues==b.residues &&
			a.validWords==b.validWords && a.scalarMask==b.scalarMask && a.cosineLower==b.cosineLower &&
			a.identityLower==b.identityLower && a.blosumLower==b.blosumLower && a.cosineAccept==b.cosineAccept &&
			a.identityUpper==b.identityUpper && a.blosumFinal==b.blosumFinal &&
			a.alignmentCalibration==b.alignmentCalibration && Arrays.equals(a.consensusRaw,b.consensusRaw) &&
			Arrays.equals(a.ddl,b.ddl) && Arrays.equals(a.spectrum,b.spectrum) && Arrays.equals(a.scalarBounds,b.scalarBounds);
	}

	private static double rawIdentity(LineParser1 lp,int field){
		final double value=Double.parseDouble(lp.parseString(field));
		require(Double.isFinite(value) && value>=0 && value<=100,"Raw identity gate must be a finite percent");
		return value;
	}
	private static boolean rawNa(LineParser1 lp,int field){
		if(lp.termEquals("NA",field)){return true;}
		return false;
	}
	/** Mirrors FamilyGateCodec.upper margin validity for evidence paths that may skip the upper formula. */
	private static void requireIdentityPadding(double multiplier,double additive){
		require(Double.isFinite(multiplier) && multiplier>=1 && Double.isFinite(additive) && additive>=0,
			"Identity margins require a finite multiplier >= 1 and a finite additive >= 0");
	}
	private static double finite(LineParser1 lp,int term){
		// These retained statistics determine bounds; use the correctly rounded parser rather than the approximate fast decimal path.
		final double value=Double.parseDouble(lp.parseString(term));
		require(Double.isFinite(value),"All-member scalar statistics must be finite for N>0"); return value;
	}

	/** Adds raw maxima before any margin/rounding; absent families retain their original calibration. */
	private static long[] mergeSupplemental(Path path,Entry[] entries,int families,int additionalQueries){
		return mergeSupplemental(path,entries,families,additionalQueries,null);
	}

	/** Replaces the legacy bootstrap alignment maxima with a complete packaged upper panel. */
	private static long[] replaceAlignmentUpper(Path path,Entry[] entries,int families,int live,int population){
		final ByteFile in=ByteFile.makeByteFile(path.toString(),false);
		final LineParser1 lp=new LineParser1('\t');final boolean[] seen=new boolean[families];
		long rows=0,observations=0;
		try{
			header(in,SUPPLEMENTAL_HEADER);
			for(byte[] line=in.nextLine();line!=null;line=in.nextLine()){
				lp.set(line);require(lp.terms()==5,"Complete alignment upper row must have five fields");
				final int rank=integer(lp,0);require(rank<families && !seen[rank],"Unknown/duplicate complete upper rank: "+rank);
				final Entry e=entries[rank];
				require(e.members>0 && e.gates.hasUpper && lp.termEquals(e.id,1),"Complete upper rank/family must match a live calibrated reference: "+rank);
				final int aligned=integer(lp,2);
				require(aligned>0 && aligned<=population,"Complete upper observations must be within the declared population: "+e.id);
				final double identity=Double.parseDouble(lp.parseString(3)),blosum=Double.parseDouble(lp.parseString(4));
				require(Double.isFinite(identity) && identity>=0 && identity<=100 && Double.isFinite(blosum),"Complete upper raw maxima must be finite, identity in percent: "+e.id);
				e.gates.negativeCount=population;e.gates.alignedCount=aligned;
				e.gates.upper[1]=identity;e.gates.upper[2]=blosum;
				seen[rank]=true;rows++;observations=Math.addExact(observations,aligned);
			}
		}finally{require(!in.close(),"Complete alignment upper reader failed");}
		require(rows==live,"Complete alignment upper must cover every live family exactly once");
		return new long[]{rows,observations};
	}
	private static long[] mergeSupplemental(Path path,Entry[] entries,int families,int additionalQueries,boolean[] unionSeen){
		final ByteFile in=ByteFile.makeByteFile(path.toString(),false);
		final LineParser1 lp=new LineParser1('\t');final boolean[] seen=new boolean[families];
		long rows=0,observations=0;
		try{
			header(in,SUPPLEMENTAL_HEADER);
			for(byte[] line=in.nextLine();line!=null;line=in.nextLine()){
				lp.set(line);require(lp.terms()==5,"Supplemental upper row must have five fields");
				final int rank=integer(lp,0);require(rank<families && !seen[rank],"Unknown/duplicate supplemental rank: "+rank);
				final Entry e=entries[rank];
				require(e.members>0 && e.gates.hasUpper && lp.termEquals(e.id,1),"Supplemental rank/family must match a live calibrated reference: "+rank);
				final int aligned=integer(lp,2);
				final int population=additionalQueries>0 ? additionalQueries : e.gates.negativeCount;
				require(aligned>0 && aligned<=population,"Supplemental observations must be distinct queries within their declared negative panel: "+e.id);
				final double identity=Double.parseDouble(lp.parseString(3)),blosum=Double.parseDouble(lp.parseString(4));
				require(Double.isFinite(identity) && identity>=0 && identity<=100 && Double.isFinite(blosum),"Supplemental raw maxima must be finite, identity in percent: "+e.id);
				e.gates.upper[1]=Math.max(e.gates.upper[1],identity);
				e.gates.upper[2]=Math.max(e.gates.upper[2],blosum);
				seen[rank]=true;rows++;observations=Math.addExact(observations,aligned);
				if(unionSeen!=null){unionSeen[rank]=true;}
			}
		}finally{require(!in.close(),"Supplemental upper reader failed");}
		assert(rows<=families && observations>=rows) : "Each supplemental family contributes at least one query observation";
		return new long[]{rows,observations};
	}

	/** Keeps evidence sources distinct while max-unioning their raw bounds; counts observations, not unique pairs. */
	private static long[] mergeAdditionalList(Path path,Entry[] entries,int families,Pins pins){
		final ByteFile in=ByteFile.makeByteFile(path.toString(),false);final LineParser1 lp=new LineParser1('\t');
		final HashSet<Path> paths=new HashSet<Path>();final boolean[] seen=new boolean[families];long observations=0;
		try{
			header(in,"path\tsha80\tnegative_count");
			for(byte[] line=in.nextLine();line!=null;line=in.nextLine()){
				lp.set(line);require(lp.terms()==3,"Additional upper list requires path, sha80 and negative_count");
				final Path table=absolute(lp.parseString(0));final int count=integer(lp,2);
				require(count>0 && paths.add(table),"Additional evidence needs a positive population and unique table path");
				pins.add(table,lp.parseString(1));
				final long[] result=mergeSupplemental(table,entries,families,count,seen);
				observations=Math.addExact(observations,result[1]);
			}
		}finally{require(!in.close(),"Additional upper list read failed");}
		require(!paths.isEmpty(),"Additional upper list must name at least one evidence table");
		int observed=0;for(boolean value:seen){if(value){observed++;}}
		return new long[]{observed,observations,paths.size()};
	}

	/** L1 frequencies; count arithmetic remains double until the final float cast. */
	public static float[] frequencies(ProteinSpectrum.Spectrum s){
		require(s!=null && s.k==2,"Dimer spectrum required");s.validate();
		final float[] f=new float[400];if(s.total>0){for(int i=0;i<f.length;i++){f[i]=(float)((double)s.counts[i]/s.total);}}return f;
	}
	private static Bucket[] index(Path path,boolean upper){
		final Bucket[] buckets=new Bucket[256];final ByteFile in=ByteFile.makeByteFile(path.toString(),false);final LineParser1 lp=new LineParser1('\t');
		try{
			header(in,FamilyGateStats.INDEX_HEADERS[upper ? 1 : 0]);
			for(byte[] line=in.nextLine();line!=null;line=in.nextLine()){
				lp.set(line);require(lp.terms()==(upper ? 10 : 11),"Malformed collection index");
				final int b=integer(lp,0);require(b<256 && buckets[b]==null,"Duplicate/out-of-range bucket");
				final int f=integer(lp,1),n=integer(lp,2),z=integer(lp,3);require((long)n+z==f,"Invalid bucket census");
				final long m=number(lp,4),aa=number(lp,5);if(upper){require(integer(lp,6)==n,"Upper row count differs from census");}
				final int field=upper ? 7 : 6;
				buckets[b]=new Bucket(f,n,z,m,aa,absolute(lp.parseString(field)),lp.parseString(field+1),
					upper ? null : absolute(lp.parseString(8)),upper ? null : lp.parseString(9));
				require(lp.termEquals("NEW",lp.terms()-1) || lp.termEquals("REUSED_ACCEPTED",lp.terms()-1) ||
					!upper && lp.termEquals("NEW_RETAINED_SCORES",lp.terms()-1),"Unknown input lineage");
			}
		}finally{require(!in.close(),"Index reader failed");}return buckets;
	}
	private static HashMap<String,Path> manifest(Path path,String[] required,String pass,Pins pins){
		final HashMap<String,Path> paths=new HashMap<String,Path>();final ByteFile in=ByteFile.makeByteFile(path.toString(),false);final LineParser1 lp=new LineParser1('\t');
		try{
			header(in,"file\tsha80");
			for(byte[] line=in.nextLine();line!=null;line=in.nextLine()){
				lp.set(line);require(lp.terms()==2,"Malformed result manifest");final String file=lp.parseString(0);
				final String logical=file.endsWith(".gz") ? file.substring(0,file.length()-3) : file;
				require(Arrays.asList(required).contains(logical) && (!file.endsWith(".gz") || logical.equals("profiles.tsv") || logical.equals("spectra.tsv") || logical.equals("cutoffs.tsv")),"Unexpected manifest file");
				final Path payload=path.getParent().resolve(file);require(paths.put(logical,payload)==null,"Duplicate logical artifact");pins.add(payload,lp.parseString(1));
			}
		}finally{require(!in.close(),"Manifest reader failed");}
		require(paths.size()==required.length,"Incomplete manifest");
		final ByteFile marker=ByteFile.makeByteFile(paths.get("PASS").toString(),false);
		try{header(marker,pass);require(marker.nextLine()==null,"Extra PASS content");}finally{require(!marker.close(),"PASS reader failed");}return paths;
	}
	private static void validateSummary(Path path,Bucket bucket,int shard){
		final ByteFile in=ByteFile.makeByteFile(path.toString(),false);final LineParser1 keys=new LineParser1('\t'),values=new LineParser1('\t');
		try{
			final byte[] a=in.nextLine(),b=in.nextLine();require(a!=null && b!=null && in.nextLine()==null,"Summary requires exactly two lines");keys.set(a);values.set(b);
			require(keys.terms()==values.terms(),"Summary columns differ");final HashMap<String,Integer> columns=new HashMap<String,Integer>();
			for(int i=0;i<keys.terms();i++){require(columns.put(keys.parseString(i),i)==null,"Duplicate summary column");}
			for(String[] pair:new String[][]{{"namespace","tail"},{"k","5"},{"seed","0"},{"exponent","6"},{"registerbits","16"},{"buckets","1024"},
				{"shard",Integer.toString(shard)},{"built_families",Integer.toString(bucket.live)},{"zero_families",Integer.toString(bucket.empty)},
				{"members",Long.toString(bucket.members)},{"residues",Long.toString(bucket.residues)}}){
				final Integer slot=columns.get(pair[0]);require(slot!=null && values.termEquals(pair[1],slot),"Reference summary mismatch: "+pair[0]);
			}
		}finally{require(!in.close(),"Summary reader failed");}
	}
	private static int upperQueries(Path path,Bucket bucket){
		final ByteFile in=ByteFile.makeByteFile(path.toString(),false);final LineParser1 lp=new LineParser1('\t');
		try{
			header(in,"families\tprocessed_families\tnegatives\tcosine_calls\talignments\ttop\tthreads\tshard\tshards\tload_seconds\tddl_score_seconds\tcosine_seconds\talign_seconds\ttotal_seconds");
			final byte[] row=in.nextLine();require(row!=null && in.nextLine()==null,"Upper summary requires exactly one data row");
			lp.set(row);require(lp.terms()==14,"Upper summary field count");
			final int families=integer(lp,0),processed=integer(lp,1),queries=integer(lp,2);
			final int alignments=integer(lp,4),top=integer(lp,5);
			final int observedShard=integer(lp,7),shards=integer(lp,8);
			require(families==bucket.live && processed==bucket.live &&
				observedShard==0 && shards==1,"Upper summary census/shard mismatch");
			require(number(lp,3)==(long)bucket.live*queries && top==500 && alignments==(long)bucket.live*top,
				"Upper summary population/count mismatch");
			return queries;
		}finally{require(!in.close(),"Upper summary reader failed");}
	}
	private static void census(Path path,Bucket bucket,int shard,Entry[] entries,HashMap<String,Entry> byId){
		final ByteFile in=ByteFile.makeByteFile(path.toString(),false);final LineParser1 lp=new LineParser1('\t');int count=0,live=0;long members=0,residues=0;
		try{
			header(in,"roster_rank\tfamily_id\tmembers\tresidues\tvalid_words");
			for(byte[] line=in.nextLine();line!=null;line=in.nextLine()){
				lp.set(line);require(lp.terms()==5,"Malformed census row");final int rank=integer(lp,0);final String id=lp.parseString(1);
				final long m=number(lp,2),aa=number(lp,3),words=number(lp,4);
				require(rank<entries.length && entries[rank]==null && !id.isEmpty() && words<=aa && (m>0 || aa==0 && words==0),"Invalid census row");
				final Entry e=new Entry(id,rank,shard,m,aa,words);require(byId.put(id,e)==null,"Duplicate census ID");entries[rank]=e;
				count++;if(m>0){live++;}members=Math.addExact(members,m);residues=Math.addExact(residues,aa);
			}
		}finally{require(!in.close(),"Census reader failed");}
		require(count==bucket.families && live==bucket.live && members==bucket.members && residues==bucket.residues,"Bucket census totals differ");
	}
	private static void profiles(Path path,int shard,HashMap<String,Entry> byId){
		final ByteFile in=ByteFile.makeByteFile(path.toString(),false);final LineParser1 lp=new LineParser1('\t');final int[] values=new int[1024],votes=new int[1024],occupied=new int[1024];
		try{
			header(in,DDLFamilyProfiles.HEADER16);
			for(byte[] line=in.nextLine();line!=null;line=in.nextLine()){
				lp.set(line);require(lp.terms()==10 && lp.termEquals("tail",0) && integer(lp,2)==1024,"Wrong profile layout");
				final Entry e=byId.get(lp.parseString(1));require(e!=null && e.bucket==shard && e.members>0 && number(lp,3)==e.members && number(lp,4)==e.residues && number(lp,5)==e.words,"Profile/census mismatch");
				final boolean literal=lp.termEquals("literal",6);require(literal || lp.termEquals("nonempty",6),"Unknown profile mode");
				final int mask=literal ? 1 : 2;require((e.profileModes&mask)==0,"Duplicate profile mode");e.profileModes|=mask;
				decimalArray(lp,7,values);decimalArray(lp,8,votes);decimalArray(lp,9,occupied);
				for(int i=0;i<1024;i++){
					require(values[i]<=65535 && votes[i]<=e.members && occupied[i]<=e.members,"Profile value/count range");
					require(values[i]==0 ? votes[i]==(literal ? e.members-occupied[i] : 0) : votes[i]>0 && votes[i]<=occupied[i],"Profile support disagreement");
					require(literal ? values[i]==0 || votes[i]>e.members-occupied[i] : (values[i]==0)==(occupied[i]==0),"Profile mode/occupancy disagreement");
				}
				if(!literal){e.ddl=new char[1024];for(int i=0;i<1024;i++){e.ddl[i]=(char)values[i];}}
			}
		}finally{require(!in.close(),"Profile reader failed");}
	}
	private static void decimalArray(LineParser1 lp,int field,int[] target){
		lp.setBounds(field);final byte[] bytes=lp.line();int start=lp.a(),slot=0;long value=0;boolean digitSeen=false;
		for(int i=start;i<=lp.b();i++){
			if(i==lp.b() || bytes[i]==','){
				require(digitSeen && slot<target.length,"Empty/excess profile array term");target[slot++]=(int)value;value=0;digitSeen=false;
			}else{final int digit=bytes[i]-'0';require(digit>=0 && digit<=9 && value<=(Integer.MAX_VALUE-digit)/10,"Invalid profile array integer");value=value*10+digit;digitSeen=true;}
		}
		require(slot==target.length,"Profile array must have 1024 terms");
	}
	private static void readConsensus(Path path,HashMap<String,Entry> byId,int total){
		//TODO: Possible bug - large early-rejected FASTA inputs can strand FastaStreamerST's producer.
		//FastaStreamerST.close() does not cancel/join; its terminal outputQueue.put at FastaStreamerST.java:232
		//can block after validation exits this loop early; this was not identified as the leaked object.
		Shared.AMINO_IN=true;final Streamer in=StreamerFactory.makeStreamer(FileFormat.testInput(path.toString(),FileFormat.FASTA,null,true,true),null,true,-1,false,true,1);
		int count=0;in.start();
		try{
			for(ListNum<Read> batch=in.nextList();batch!=null;batch=in.nextList()){
				for(Read read:batch.list){
					require(read!=null && read.id!=null && read.bases!=null && read.mate==null,"Invalid consensus record");
					int end=0;while(end<read.id.length() && !Character.isWhitespace(read.id.charAt(end))){end++;}
					final Entry e=byId.get(read.id.substring(0,end));require(e!=null && e.consensus==null,"Unknown/duplicate consensus ID");
					require(read.bases.length>0 && read.bases.length<=MAX_CONSENSUS_AA,
						"Consensus outside 1.."+MAX_CONSENSUS_AA+"AA reference contract for "+e.id+": "+read.bases.length);
					new ProteinSequence(e.id,read.bases);e.consensus=read.bases;count++;
				}
			}
		}finally{require(!ReadWrite.closeStream(in),"Consensus reader failed");}
		require(count==total,"Consensus does not cover full roster including empty families");
	}
	private static HashSet<String> survivors(Path path){
		final ByteFile in=ByteFile.makeByteFile(path.toString(),false);
		final LineParser1 lp=new LineParser1('\t');
		final HashSet<String> keep=new HashSet<String>();
		try{
			final byte[] header=in.nextLine();require(header!=null,"Survivor list requires a header");
			lp.set(header);final int column=column(lp,"family_id");
			for(byte[] line=in.nextLine();line!=null;line=in.nextLine()){
				lp.set(line);require(lp.terms()>column,"Malformed survivor row");
				final String id=lp.parseString(column);require(!id.isEmpty() && keep.add(id),"Empty/duplicate survivor: "+id);
			}
		}finally{require(!in.close(),"Survivor list read failed");}
		require(!keep.isEmpty(),"Survivor list must not be empty");
		return keep;
	}

	private static void writeSubsetReferences(Path parent,Path out,FamilyReference[] records){
		final HashMap<String,Integer> denseRanks=new HashMap<String,Integer>(records.length*2);
		for(int i=0;i<records.length;i++){require(denseRanks.put(records[i].id,i)==null,"Duplicate survivor reference ID");}
		final ByteFile in=ByteFile.makeByteFile(parent.toString(),false);
		final ByteStreamWriter writer=new ByteStreamWriter(out.toString(),false,false,false);writer.start();
		final LineParser1 lp=new LineParser1('\t');final ByteBuilder buffer=new ByteBuilder(4096);
		int written=0;
		boolean complete=false;
		try{
			final byte[] header=in.nextLine();require(header!=null,"Missing reference header");
			writer.print(rewriteCensus(header,records.length,buffer));
			for(int i=1;i<FamilyReferenceCodec.LINES_PER_FAMILY;i++){
				final byte[] line=in.nextLine();require(line!=null,"Incomplete reference header labels");
				writer.print(buffer.clear().append(line).nl());
			}
			for(byte[] meta=in.nextLine();meta!=null;meta=in.nextLine()){
				final byte[][] payload=new byte[FamilyReferenceCodec.LINES_PER_FAMILY-1][];
				for(int i=0;i<payload.length;i++){
					payload[i]=in.nextLine();require(payload[i]!=null,"Incomplete reference record");
				}
				lp.set(meta);require(lp.terms()==6 && lp.termEquals("#family",0),"Malformed #family metadata");
				final Integer rank=denseRanks.remove(lp.parseString(1));
				if(rank==null){continue;}
				writer.print(rewriteRank(lp,rank.intValue(),buffer));
				for(byte[] line:payload){writer.print(buffer.clear().append(line).nl());}
				written++;
			}
			complete=true;
		}finally{
			final boolean readError=in.close(),writeError=writer.poisonAndWait();
			if(complete && (readError || writeError)){throw new IllegalStateException("Reference subset I/O failed; no PASS");}
		}
		require(written==records.length && denseRanks.isEmpty(),"Parent stream missed a survivor reference");
	}

	private static ByteBuilder rewriteCensus(byte[] line,int families,ByteBuilder out){
		final String[] terms=new String(line,StandardCharsets.US_ASCII).split("\t",-1);
		require((terms.length==11 || terms.length==14) && terms[0].equals("#family_reference"),"Malformed reference header");
		terms[5]=terms[6]=Integer.toString(families);terms[7]="0";
		out.clear().append(terms[0]);
		for(int i=1;i<terms.length;i++){out.tab().append(terms[i]);}
		return out.nl();
	}

	private static ByteBuilder rewriteRank(LineParser1 lp,int rank,ByteBuilder out){
		final byte[] line=lp.line();
		out.clear();
		for(int i=0;i<lp.terms();i++){
			if(i>0){out.tab();}
			lp.setBounds(i);
			if(i==2){out.append(rank);}
			else{out.append(line,lp.a(),lp.b()-lp.a());}
		}
		return out.nl();
	}

	/** Copies only tail/width1024 rows for the accepted survivor set; register bytes stay unchanged. */
	private static int filterProfiles(Path list,FamilyReference[] records,Path profilesOut,Path summaryOut,Pins pins) throws Exception{
		final HashMap<String,Integer> ordinals=new HashMap<String,Integer>(records.length*2);
		for(int i=0;i<records.length;i++){require(ordinals.put(records[i].id,i)==null,"Duplicate survivor reference ID");}
		final boolean[] seen=new boolean[records.length];
		final ByteStreamWriter writer=new ByteStreamWriter(profilesOut.toString(),false,false,false);writer.start();
		final ByteBuilder buffer=new ByteBuilder(8192);
		ProfileSummary summary=null;
		String profileHeader=null;int profileTerms=0,profiled=0;
		long profileMembers=0,profileResidues=0,profileWords=0;
		final ByteFile in=ByteFile.makeByteFile(list.toString(),false);final LineParser1 lp=new LineParser1('\t');
		final HashSet<Path> paths=new HashSet<Path>();
		boolean complete=false;
		try{
			for(byte[] entry=in.nextLine();entry!=null;entry=in.nextLine()){
				lp.set(entry);require(lp.terms()==3,"Profile list requires path, sha80 and summary_sha80");
				Path path=Paths.get(lp.parseString(0));
				if(!path.isAbsolute()){path=list.getParent().resolve(path);}
				path=path.normalize();require(paths.add(path),"Duplicate profile path");
				pins.add(path,lp.parseString(1));
				final Path sourceSummary=path.getParent().resolve("summary.tsv");pins.add(sourceSummary,lp.parseString(2));
				summary=summary==null ? new ProfileSummary(sourceSummary) : summary.compatible(sourceSummary);
				final ByteFile profiles=ByteFile.makeByteFile(path.toString(),false);final LineParser1 row=new LineParser1('\t');
				try{
					final byte[] header=profiles.nextLine();require(header!=null,"Missing profile header");
					final String text=new String(header,StandardCharsets.US_ASCII);
					if(profileHeader==null){
						require(text.equals(DDLFamilyProfiles.HEADER16) || text.equals(DDLFamilyProfiles.HEADER16+"\tpositions"),
							"Subset profile input must be legacy16 DDL rows");
						profileHeader=text;writer.print(text+"\n");
						profileTerms=text.equals(DDLFamilyProfiles.HEADER16) ? 10 : 11;
					}else{require(profileHeader.equals(text),"Mixed profile headers");}
					for(byte[] line=profiles.nextLine();line!=null;line=profiles.nextLine()){
						row.set(line);require(row.terms()==profileTerms,"Malformed profile row");
						if(!row.termEquals("tail",0) || row.parseInt(2)!=1024){continue;}
						final Integer ordinal=ordinals.get(row.parseString(1));
						if(ordinal==null){continue;}
						buffer.clear().append(line).nl();writer.print(buffer);
						if(row.termEquals("nonempty",6)){
							require(!seen[ordinal],"Duplicate survivor nonempty profile");
							profileMembers=Math.addExact(profileMembers,number(row,3));
							profileResidues=Math.addExact(profileResidues,number(row,4));
							profileWords=Math.addExact(profileWords,number(row,5));
							seen[ordinal]=true;profiled++;
						}else{require(row.termEquals("literal",6),"Unknown profile mode");}
					}
				}finally{require(!profiles.close(),"Profile reader failed: "+path);}
			}
			complete=true;
		}finally{
			final boolean readError=in.close(),writeError=writer.poisonAndWait();
			if(complete && (readError || writeError)){throw new IllegalStateException("Profile subset I/O failed; no PASS");}
		}
		for(int i=0;i<seen.length;i++){require(seen[i],"Missing survivor profile: "+records[i].id);}
		require(summary!=null,"Profile list must not be empty");
		summary.write(summaryOut,profiled,profileMembers,profileResidues,profileWords);
		return profiled;
	}
	private static final class Entry{
		Entry(String id,int rank,int bucket,long members,long residues,long words){this.id=id;this.rank=rank;this.bucket=bucket;this.members=members;this.residues=residues;this.words=words;}
		final String id;final int rank,bucket;final long members,residues,words;
		int profileModes;char[] ddl;float[] spectrum;byte[] consensus;FamilyGateExport.Row gates;
	}
	private static final class Bucket{
		Bucket(int f,int n,int z,long m,long aa,Path manifest,String pin,Path lower,String lowPin){families=f;live=n;empty=z;members=m;residues=aa;this.manifest=manifest;this.pin=pin;lowerManifest=lower;lowerPin=lowPin;}
		boolean same(Bucket b){return b.families==families && b.live==live && b.empty==empty && b.members==members && b.residues==residues;}
		final int families,live,empty;final long members,residues;final Path manifest,lowerManifest;final String pin,lowerPin;
	}
	private static final class ProfileSummary{
		ProfileSummary(Path path) throws Exception{
			final java.util.List<String> lines=Files.readAllLines(path,StandardCharsets.US_ASCII);
			require(lines.size()==2,"Profile summary requires header and one data row");
			header=lines.get(0).split("\t",-1);values=lines.get(1).split("\t",-1);
			require(header.length==values.length,"Profile summary field count");
			for(int i=0;i<header.length;i++){require(columns.put(header[i],i)==null,"Duplicate profile summary column");}
			for(String key:new String[]{"namespace","k","seed","exponent","sketch"}){
				require(columns.containsKey(key),"Profile summary missing "+key);
			}
		}
		ProfileSummary compatible(Path path) throws Exception{
			final ProfileSummary next=new ProfileSummary(path);
			for(String key:new String[]{"namespace","k","seed","exponent","registerbits","buckets","sketch","partitions","overlap","tag_bits","segment_words","position_basis","blacklist_sha80","blacklist_words"}){
				final Integer a=columns.get(key),b=next.columns.get(key);
				require((a==null)==(b==null) && (a==null || values[a].equals(next.values[b])),"Profile summary mismatch: "+key);
			}
			return this;
		}
		void write(Path path,int families,long members,long residues,long words) throws Exception{
			set("shard","0");set("shards","1");set("manifest_families",families);set("built_families",families);
			set("zero_families","0");set("members",members);set("residues",residues);set("valid_words",words);
			set("load_seconds","0");set("build_seconds","0");set("elapsed_seconds","0");
			FamilyReferenceExport.write(path,String.join("\t",header)+"\n"+String.join("\t",values)+"\n");
		}
		void set(String key,long value){set(key,Long.toString(value));}
		void set(String key,String value){final Integer column=columns.get(key);if(column!=null){values[column]=value;}}
		final String[] header,values;final HashMap<String,Integer> columns=new HashMap<String,Integer>();
	}
	private static final class Pins{
		void add(Path path,String expected){
			expected=DigestSuffix.requireSuffix(expected,"input");require(Files.isRegularFile(path) && !Files.isSymbolicLink(path),"Require regular input: "+path);
			final String old=values.put(path,expected);require(old==null || old.equals(expected),"Conflicting input pins");
			require(DigestSuffix.file(path.toString()).equals(expected),"Input sha80 mismatch: "+path);
		}
		void verify(){for(java.util.Map.Entry<Path,String> p:values.entrySet()){require(DigestSuffix.file(p.getKey().toString()).equals(p.getValue()),"Input changed: "+p.getKey());}}
		final LinkedHashMap<Path,String> values=new LinkedHashMap<Path,String>();
	}
	private static long number(LineParser1 lp,int field){
		lp.setBounds(field);require(lp.a()<lp.b(),"Empty decimal field");long x=0;
		for(int i=lp.a();i<lp.b();i++){final int d=lp.line()[i]-'0';require(d>=0 && d<=9 && x<=(Long.MAX_VALUE-d)/10,"Invalid/overflowed count");x=x*10+d;}return x;
	}
	private static int integer(LineParser1 lp,int field){long x=number(lp,field);require(x<=Integer.MAX_VALUE,"Integer overflow");return (int)x;}
	private static Path input(HashMap<String,String> options,String key,Pins pins){final Path path=absolute(required(options,key));pins.add(path,required(options,key+"sha80"));return path;}
	private static Path structuralInput(HashMap<String,String> options,String key,StructuralInputs inputs) throws Exception{
		return inputs.track(absolute(required(options,key)));
	}
	private static String structuralInputRecords(Path[] paths) throws Exception{
		final ByteBuilder out=new ByteBuilder("logical\tpath\tsize\tmodified_time\tfile_key\n");
		for(Path path:paths){
			final BasicFileAttributes attrs=Files.readAttributes(path,BasicFileAttributes.class,LinkOption.NOFOLLOW_LINKS);
			require(attrs.isRegularFile(),"Require regular structural input: "+path);
			out.append(path.getFileName().toString()).tab().append(path.toString()).tab().append(attrs.size()).tab()
				.append(attrs.lastModifiedTime().toString()).tab().append(String.valueOf(attrs.fileKey())).nl();
		}
		return out.toString();
	}
	private static String structuralInputRecord(Path path) throws Exception{
		final BasicFileAttributes attrs=Files.readAttributes(path,BasicFileAttributes.class,LinkOption.NOFOLLOW_LINKS);
		require(attrs.isRegularFile(),"Require regular structural input: "+path);
		return "path\tsize\tmodified_time\tfile_key\n"+path+"\t"+attrs.size()+"\t"+
			attrs.lastModifiedTime().toString()+"\t"+String.valueOf(attrs.fileKey())+"\n";
	}
	private static CompactCensus compactCensus(String view){
		if("short".equals(view)){return new CompactCensus(view,159994,159994,0,11149289,40,99,28,20,"0.985");}
		if("tails".equals(view)){return new CompactCensus(view,120000,120000,0,117316888,100,100,60,43,"1");}
		throw new IllegalArgumentException("Compact view must be short or tails");
	}
	private static ProteinKmerGeometry compactGeometry(String aabits){
		final ProteinKmerGeometry geometry=ProteinKmerGeometry.make(6,aabits);
		final String name=geometry.name();
		require(name.startsWith("dimer6+1,3,5+"),"Accepted compact alphabet must be the Q70 dimer6+1,3,5 geometry");
		final String suffix=name.substring("dimer6+1,3,5+".length());
		int groups=1;
		for(int i=0;i<suffix.length();i++){if(suffix.charAt(i)=='/'){groups++;}}
		require(groups==120,"Accepted compact alphabet must preserve exactly 120 ordered Q70 groups");
		return geometry;
	}
	private static String compactAabits(Path path) throws Exception{
		final java.util.List<String> lines=Files.readAllLines(path,StandardCharsets.US_ASCII);
		require(lines.size()==2 && lines.get(0).equals("#bbgap_alphabet_config\t1"),
			"Compact alphabet resource requires its schema and one aabits line");
		final String[] pair=Parse.splitOnFirst(lines.get(1),'=');
		require(pair.length==2 && pair[0].equals("aabits") && !pair[1].isEmpty(),
			"Compact alphabet resource must contain one aabits mapping");
		return pair[1];
	}
	private static void addScalarBounds(LinkedHashMap<String,String> fields,Path path,CompactCensus census){
		final ByteFile in=ByteFile.makeByteFile(path.toString(),false);
		final LineParser1 lp=new LineParser1('\t');
		final HashSet<String> seen=new HashSet<String>();
		int rows=0;
		try{
			header(in,"scalar\trows\tmin\tmax\tmin_family\tmax_family");
			for(byte[] line=in.nextLine();line!=null;line=in.nextLine()){
				lp.set(line);require(lp.terms()==6,"Scalar bounds require six fields");
				final String name=lp.parseString(0);
				require(Arrays.asList(SCALAR_NAMES).contains(name) && seen.add(name),"Unknown or duplicate scalar bound: "+name);
				require(number(lp,1)==COMBINED25_SCALAR_FAMILIES,
					"Compact scalar bounds must come from the shared 79565-family combined25 source: "+name);
				final double min=Double.parseDouble(lp.parseString(2)),max=Double.parseDouble(lp.parseString(3));
				require(Double.isFinite(min) && Double.isFinite(max) && min<max,"Scalar bounds must be finite and increasing: "+name);
				fields.put(name+"_min",lp.parseString(2));
				fields.put(name+"_max",lp.parseString(3));
				rows++;
			}
		}finally{require(!in.close(),"Scalar bounds reader failed");}
		require(rows==SCALAR_NAMES.length,"Compact scalar bounds must cover all five combined25 scalars");
	}
	private static String compactLocator(HashMap<String,String> options,String key,String fallback){
		final String value=options.containsKey(key) ? required(options,key) : fallback;
		final Path path=Paths.get(value).normalize();
		require(value.indexOf('\\')<0 && value.charAt(0)!='?' && !path.isAbsolute() && path.getNameCount()>0 &&
			path.equals(Paths.get(value)),"Compact shared resource locator must be normalized and package-root relative: "+value);
		for(Path term:path){
			final String text=term.toString();
			require(!text.equals("..") && !text.equals("."),"Compact shared resource locator cannot escape or contain dots: "+value);
		}
		return value;
	}
	private static Path absolute(String text){final Path path=Paths.get(text);require(path.isAbsolute() && path.normalize().equals(path),"Require normalized absolute path");return path;}
	private static int column(LineParser1 lp,String name){
		for(int i=0;i<lp.terms();i++){if(lp.termEquals(name,i)){return i;}}
		throw new IllegalArgumentException("Missing column: "+name);
	}
	private static String required(HashMap<String,String> options,String key){final String value=options.get(key);require(value!=null && !value.isEmpty(),"Missing option: "+key);return value;}
	private static double identityPadding(HashMap<String,String> options,String key,double fallback,double minimum){
		final double value=options.containsKey(key) ? Double.parseDouble(options.get(key)) : fallback;
		require(Double.isFinite(value) && value>=minimum,"Identity padding is below its finite minimum: "+key);
		return value;
	}
	private static void header(ByteFile input,String text){require(Arrays.equals(input.nextLine(),text.getBytes(StandardCharsets.US_ASCII)),"Unexpected header/PASS: "+text);}
	private static void write(Path path,String text) throws Exception{Files.write(path,text.getBytes(StandardCharsets.US_ASCII),StandardOpenOption.CREATE_NEW);}
	private static void require(boolean condition,String message){if(!condition){throw new IllegalArgumentException(message);}}
	public static final String PASS="FAMILY_REFERENCE_EXPORT_PASS";
	public static final String RAISE_PASS="FAMILY_REFERENCE_RAISE_PASS";
	public static final String MEMBER_START_PASS="FAMILY_REFERENCE_MEMBER_START_PASS";
	public static final String IDENTITY_PADDING_PASS="FAMILY_REFERENCE_IDENTITY_PADDING_PASS";
	public static final String MONOMER_PASS="FAMILY_REFERENCE_MONOMER_PASS";
	public static final String SUBSET_PASS="FAMILY_REFERENCE_SUBSET_PASS";
	public static final String COMPACT_PASS="FAMILY_REFERENCE_COMPACT_BINDING_PASS";
	public static final String ACCEPT_OVERRIDE_HEADER="roster_rank\tfamily_id\tidentity_accept\tblosum_final";
	public static final String SUPPLEMENTAL_HEADER="roster_rank\tfamily_id\taligned_count\tidentity_pct_upper\tblosum_a_sq_upper";
	/** Existing collection lower/bootstrap recipes pass maxreferenceaa=1000. */
	private static final int MAX_CONSENSUS_AA=1000;
	private static final int COMBINED25_SCALAR_FAMILIES=79565;
	private static final String[] SCALAR_NAMES={"hydrophobicity","net_charge","aromaticity","fraction_small","pI_tail"};
	private static final String[] OPTIONS={"referenceindex","referenceindexsha80","upperindex","upperindexsha80","consensus","consensussha80","out","maxfamilies","supplementalupper","supplementaluppersha80","scalars","scalarssha80","scalarz","mode","parent","parentsha80","overrides","overridessha80","overrideslist","overrideslistsha80","identityfloor","threads","additionalupper","additionaluppersha80","additionalqueries","additionalupperlist","additionalupperlistsha80","blosumfloor","identityuppermultiplier","identityupperadditive","floorreference","floorreferencesha80","survivors","survivorssha80","ddlprofiles","ddlprofilessha80","alignmentupper","alignmentuppersha80","alignmentqueries","rawcalibration","rawcalibrationsha80","view","references","profilevectors","profilesummary","lengthmedians","familyhydro","familygc","familygcmap","scalarbounds","alphabet","blacklist","neutralblacklist","matrixresource","matrixbinding","alphabetlocator","blacklistlocator","matrixlocator","outfamilies"};
	private static final String[] COMPACT_OPTIONS={"mode","view","out","references","ddlprofiles","profilevectors","profilesummary",
		"lengthmedians","familyhydro","familygc","familygcmap","scalarbounds","alphabet","blacklist","neutralblacklist","matrixresource",
		"matrixbinding","alphabetlocator","blacklistlocator","matrixlocator","outfamilies","threads"};
	private static final String[] REFERENCE_FILES={"profiles.tsv","spectra.tsv","members.tsv","summary.tsv","recipe.txt","timing.tsv","PASS"};
	private static final String[] LOWER_FILES={"cutoffs.tsv","spectra.tsv","summary.tsv","recipe.txt","timing.tsv","PASS"};
	private static final String[] UPPER_FILES={"uppers.tsv","summary.tsv","timing.tsv","recipe.txt","PASS"};
	private static final class CompactCensus{
		CompactCensus(String view_,int families_,int live_,int empty_,long members_,int minAa_,int maxAa_,
				int alignmentMinSpan_,int alignmentMinPaired_,String marginScale_){
			view=view_;families=families_;live=live_;empty=empty_;members=members_;minAa=minAa_;maxAa=maxAa_;
			alignmentMinSpan=alignmentMinSpan_;alignmentMinPaired=alignmentMinPaired_;marginScale=marginScale_;
		}
		final String view,marginScale;
		final int families,live,empty,minAa,maxAa,alignmentMinSpan,alignmentMinPaired;
		final long members;
	}

	private static final class ZipState{
		ZipState(){
			threads=Shared.threads(); buffers=Shared.READ_BUFFER_NUM_BUFFERS();
			usePigz=ReadWrite.USE_PIGZ; useUnpigz=ReadWrite.USE_UNPIGZ; useBgzf=ReadWrite.USE_BGZF;
			useBgzip=ReadWrite.USE_BGZIP; forceBgzip=ReadWrite.FORCE_BGZIP;
			allowNativeBgzf=ReadWrite.ALLOW_NATIVE_BGZF;
			zipLevel=ReadWrite.ZIPLEVEL; zipThreads=ReadWrite.MAX_ZIP_THREADS();
			zipThreadMult=ReadWrite.ZIP_THREAD_MULT(); setZipThreads=ReadWrite.SET_ZIP_THREADS;
			setZipThreadMult=ReadWrite.SET_ZIP_THREAD_MULT;
		}
		void forceBgzip(int zipThreads_){
			requireExternalBgzip();
			Shared.setThreads(zipThreads_);
			ReadWrite.USE_PIGZ=false; ReadWrite.USE_UNPIGZ=false;
			ReadWrite.USE_BGZIP=true; ReadWrite.USE_BGZF=true; ReadWrite.FORCE_BGZIP=true;
			ReadWrite.ALLOW_NATIVE_BGZF=false;
			ReadWrite.ZIPLEVEL=9; ReadWrite.setZipThreadMult(1);
			ReadWrite.setZipThreads(zipThreads_);
		}
		void restore(){
			Shared.setThreads(threads);
			Shared.setBuffers(buffers);
			ReadWrite.USE_PIGZ=usePigz; ReadWrite.USE_UNPIGZ=useUnpigz; ReadWrite.USE_BGZF=useBgzf;
			ReadWrite.USE_BGZIP=useBgzip; ReadWrite.FORCE_BGZIP=forceBgzip;
			ReadWrite.ALLOW_NATIVE_BGZF=allowNativeBgzf;
			ReadWrite.ZIPLEVEL=zipLevel; ReadWrite.setZipThreadMult(zipThreadMult);
			ReadWrite.setZipThreads(zipThreads);
			ReadWrite.SET_ZIP_THREAD_MULT=setZipThreadMult;
			ReadWrite.SET_ZIP_THREADS=setZipThreads;
		}
		static void requireExternalBgzip(){
			require(Data.BGZIP() && Data.BGZIP_VERSION_levelFlag && Data.BGZIP_VERSION_threadsFlag,
				"Explicit compact BGZF output requires external bgzip with -l and -@ support");
		}
		final int threads,buffers,zipLevel,zipThreads;
		final float zipThreadMult;
		final boolean usePigz,useUnpigz,useBgzf,useBgzip,forceBgzip,allowNativeBgzf,setZipThreads,setZipThreadMult;
	}
}
