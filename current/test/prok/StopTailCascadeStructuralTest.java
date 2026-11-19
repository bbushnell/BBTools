package test.prok;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;

import ddl.DDLFamilyProfiles;
import ddl.ProteinKmerGeometry;
import dna.GeneticCode;
import prok.StopTailCascadeResources;
import prok.StopTailCascadeSupport;
import prok.StopTailCandidateScan;
import prot.FamilyAlignmentMetric;
import prot.FamilyCalibration.VerificationMode;
import prot.FamilyReference;
import prot.FamilyReferenceCodec;
import prot.FamilyScoringPackage;
import prot.FamilyTailCascade;
import prot.FamilyTailCascade.Assignment;
import prot.ProteinSequence;
import shared.Shared;
import structures.ByteBuilder;

/** No-hash caller resource binding and actual two-set override propagation.
 * Literal three-hit profiles and 61 families distinguish floor and shortlist policies.
 * @author Keqing
 */
public final class StopTailCascadeStructuralTest{
	/** Reuse the exact pinned structural fixture for information preflight integration.
	 * The three-hit profiles deliberately reject under floor5; this builds no fitted model. */
	public static Path writeInformationResources(Path root)throws Exception{
		final byte[] raw=new byte[100];
		raw[0]='M';
		for(int i=1; i<raw.length; i++){raw[i]=(byte)ACIDS.charAt((i*7+i/20)%20);}
		return writeInformationResources(root, raw, false);
	}
	/** Accepted identical-reference variant for an independent actual-driver fitting oracle. */
	public static Path writeInformationResources(Path root, byte[] raw, boolean accept)throws Exception{
		return writeInformationResources(root, raw, accept, 28);
	}
	/** Length-collector fixture: LONG may use span70; existing callers retain28. */
	public static Path writeInformationResources(Path root, byte[] raw, boolean accept, int longSpan)throws Exception{
		check(raw.length==100, "The composite fixture requires a literal100-residue reference");
		Files.createDirectory(root);
		final Fixture full=fixture(root.resolve("long"), raw, 100, LONG_PIN, 1, accept, longSpan);
		final Fixture tiny=fixture(root.resolve("short"), Arrays.copyOfRange(raw, 1, 100), 40, SHORT_PIN, .985, accept);
		final Path child=write(tiny.root.resolve("manifest.tsv"), "#bbgap_tail_cascade\t2\n"+tiny.settings);
		return write(full.root.resolve("manifest.tsv"), "#bbgap_tail_cascade\t7\n"+full.settings+
			"short40manifest\t"+child+"\nshort40sha80\t"+SHORT_PIN+"\n");
	}

	/** Builds fresh LONG/SHORT fixtures and checks routing, overrides, alias protection and input stability. */
	public static void main(String[] args)throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Supply a fresh absolute Dori fixture directory");}
		Shared.setThreads(1);
		final Path root=Paths.get(args[0]);
		check(root.isAbsolute() && !Files.exists(root), "Fixture root must be new and absolute");
		Files.createDirectory(root);
		final byte[] raw=new byte[100];
		raw[0]='M';
		for(int i=1; i<raw.length; i++){raw[i]=(byte)ACIDS.charAt((i*7+i/20)%20);}
		final Fixture full=fixture(root.resolve("long"), raw, 100, LONG_PIN, 1);
		final Fixture tiny=fixture(root.resolve("short"), Arrays.copyOfRange(raw, 1, 100), 40, SHORT_PIN, .985);
		final Path child=write(tiny.root.resolve("manifest.tsv"), "#bbgap_tail_cascade\t2\n"+tiny.settings);
		final Path manifest=write(full.root.resolve("manifest.tsv"), "#bbgap_tail_cascade\t7\n"+full.settings+
			"short40manifest\t"+child+"\nshort40sha80\t"+SHORT_PIN+"\n");
		final byte[] manifestBefore=Files.readAllBytes(manifest), childBefore=Files.readAllBytes(child);
		final StopTailCascadeResources resources=StopTailCascadeResources.readStructural(manifest.toString(), LONG_PIN);
		check(resources.sha80.equals(LONG_PIN) && resources.queryLengths().equals("40-99,100") && resources.allFamilySets,
			"Unmeasured identities and complete version-seven routing are retained without digest verification");
		check(StopTailCascadeResources.readStructural(manifest.toString(), METADATA_PIN).sha80.equals(METADATA_PIN),
			"Structural loading also preserves an explicitly inherited root identity");
		final Path aliasList=Files.createLink(root.resolve("long-list-alias.tsv"), full.root.resolve("list.tsv"));
		int aliasCase=0;
		for(Path sameList:new Path[]{full.root.resolve("list.tsv"), aliasList}){
			final Path wrongChild=write(tiny.root.resolve("same-index-"+aliasCase+".tsv"), "#bbgap_tail_cascade\t2\n"+
				tiny.settings.replace("ddlprofiles\t"+tiny.root.resolve("list.tsv"), "ddlprofiles\t"+sameList));
			final Path wrongRoot=write(full.root.resolve("same-index-"+aliasCase+".tsv"), "#bbgap_tail_cascade\t7\n"+
				full.settings+"short40manifest\t"+wrongChild+"\nshort40sha80\tunmeasured\n");
			reject(()->StopTailCascadeResources.readStructural(wrongRoot.toString(), LONG_PIN), "distinct short-family index");
			aliasCase++;
		}
		final Assignment baseline=resources.load(1).newWorker().assign("baseline", raw);
		check(!baseline.assigned() && baseline.rawHitFloor==5 && baseline.admittedFamilies==0 && baseline.dpAttempts==0,
			"Sixty-one literal three-hit profiles must remain below the original floor five");
		for(int cap:new int[]{50, 100}){
			final Counts counts=new Counts();
			final Assignment actual=resources.loadWithOverrides(1, cap, 3, 0).newWorker().assignTrace("override", raw, counts);
			final int selected=Math.min(cap, 61);
			check(actual.rawHitFloor==3 && actual.admittedFamilies==61 && counts.selected==selected && actual.dpAttempts==selected,
				"Floor three admits all61; the strict requested cap controls actual attempts when no upper gate passes");
			check(counts.lower==96 && counts.identityUpper==110 && counts.blosumUpper>1,
				"Floor/cap overrides leave the observed lower and direct upper gates unchanged");
			same(actual, full.pack.loadWithOverrides(1, cap, 3).newWorker().assign("prior_api", raw));
		}
		final Counts reduced=new Counts();
		resources.loadWithOverrides(1, 100, 3, 4).newWorker().assignTrace("lower", raw, reduced);
		check(reduced.lower==92 && reduced.identityUpper==110 && reduced.blosumUpper>1,
			"Lower relaxation is four percentage points while direct upper gates stay fixed");
		same(baseline, resources.load(1).newWorker().assign("still_original", raw));
		check(Arrays.equals(manifestBefore, Files.readAllBytes(manifest)) && Arrays.equals(childBefore, Files.readAllBytes(child)),
			"No override rewrites a source manifest");
		final byte[] shortQuery=Arrays.copyOfRange(raw, 1, 100);
		final FamilyTailCascade.LoadedClassifier shortControl=tiny.pack.loadWithOverrides(1, 50, 3).withEvaluationMarginScale(.985);
		final Counts beforeScores=new Counts(), afterScores=new Counts();
		final Assignment beforeLower=shortControl.newWorker().assignTrace("short_before", shortQuery, beforeScores),
			afterLower=shortControl.withIdentityLowerReduction(4).newWorker().assignTrace("short_after", shortQuery, afterScores);
		check(beforeScores.identity>92 && beforeScores.identity<93 && beforeLower.terminal.equals("identity_reject") &&
			afterLower.terminal.equals("blosum_reject") && beforeScores.identity==afterScores.identity && beforeScores.blosum==afterScores.blosum &&
			beforeScores.lower==96 && afterScores.lower==92 && beforeScores.identityUpper==afterScores.identityUpper && beforeScores.blosumUpper==afterScores.blosumUpper,
			"Seven interior substitutions in99 residues witness SHORT lower96-to92 relaxation without changing alignment or upper gates");
		final byte[] dna=dna(raw);
		try(StopTailCascadeSupport support=resources.loadSupportWithOverrides(1, 50, 3, 0)){
			final StopTailCandidateScan.Scan scan=support.scan("unrelaxed", dna, GeneticCode.forTable(11));
			final StopTailCascadeSupport.ContigEvidence evidence=support.classify(scan);
			boolean witness=false;
			for(int i=0; i<scan.candidates.size(); i++){
				if(Arrays.equals(scan.candidates.get(i).copyTail(), raw)){
					witness=true;
					same(beforeLower, evidence.familySetAssignment(i, true));
				}
			}
			check(witness, "The resource control must reach the unrelaxed SHORT rejection");
		}
		for(int cap:new int[]{50, 100}){
			final FamilyTailCascade.LoadedClassifier expectedLong=full.pack.loadWithOverrides(1, cap, 3).withIdentityLowerReduction(4);
			final FamilyTailCascade.LoadedClassifier expectedShort=tiny.pack.loadWithOverrides(1, cap, 3)
				.withEvaluationMarginScale(.985).withIdentityLowerReduction(4);
			try(StopTailCascadeSupport support=resources.loadSupportWithOverrides(1, cap, 3, 4)){
				final StopTailCandidateScan.Scan scan=support.scan("synthetic", dna, GeneticCode.forTable(11));
				final StopTailCascadeSupport.ContigEvidence evidence=support.classify(scan);
				boolean witness=false;
				for(int i=0; i<scan.candidates.size(); i++){
					final byte[] tail=scan.candidates.get(i).copyTail();
					final byte[] suffix=Arrays.copyOfRange(tail, Math.max(0, tail.length-99), tail.length);
					same(expectedShort.newWorker().assign("short_direct", suffix), evidence.familySetAssignment(i, true));
					if(tail.length==100){same(expectedLong.newWorker().assign("long_direct", tail), evidence.familySetAssignment(i, false));}
					if(Arrays.equals(tail, raw)){
						witness=true;
						check(evidence.familySetAssignment(i, true).terminal.equals("blosum_reject"),
							"SHORT must receive the lower reduction too; leaving its old lower would produce identity_reject");
						check(evidence.familySetAssignment(i, true).dpAttempts==Math.min(cap, 61) &&
							evidence.familySetAssignment(i, false).dpAttempts==Math.min(cap, 61), "Both sets receive the same cap and floor");
					}
				}
				check(witness, "The actual scanner must recover the literal100-residue test query");
			}
		}
		resources.verify();
		for(Path file:new Path[]{manifest, child, full.root.resolve("list.tsv"), full.root.resolve("profiles.tsv"),
			full.root.resolve("summary.tsv"), full.root.resolve("references.tsv"), tiny.root.resolve("references.tsv"),
			tiny.root.resolve("medians.tsv"), tiny.root.resolve("hydro.tsv"), tiny.root.resolve("bounds.tsv")}){
			final StopTailCascadeResources frozen=StopTailCascadeResources.readStructural(manifest.toString(), LONG_PIN);
			final byte[] before=Files.readAllBytes(file), changed=Arrays.copyOf(before, before.length+1);
			changed[before.length]='\n';
			try{
				Files.write(file, changed);
				reject(frozen::verify, "Structural input changed");
			}finally{Files.write(file, before);}
		}
		final StopTailCascadeResources fresh=StopTailCascadeResources.readStructural(manifest.toString(), LONG_PIN);
		final Path alias=Files.createLink(root.resolve("payload-alias.tsv"), tiny.root.resolve("profiles.tsv"));
		reject(()->fresh.protectOutputs(new String[]{alias.toString()}), "aliases stop-tail resource");
		reject(()->fresh.loadWithOverrides(1, -1, 3, 0), "Caller overrides");
		reject(()->fresh.loadWithOverrides(1, 50, -2, 0), "Caller overrides");
		reject(()->fresh.loadWithOverrides(1, 50, 3, Float.NaN), "Caller overrides");
		reject(()->StopTailCascadeResources.readStructural(manifest.toString(), "invalid"), "declared identity");
		final Path filtered=write(root.resolve("filtered.tsv"), new String(manifestBefore, StandardCharsets.US_ASCII)+
			"blacklist\tmissing.faa\nblacklistsha80\t"+LONG_PIN+"\n");
		// An enabled blacklist without policy fields must fail before opening its payload.
		reject(()->StopTailCascadeResources.readStructural(filtered.toString(), LONG_PIN), "requires the complete explicit blacklist policy");
		plainShort(root, full, tiny, raw, baseline);
		System.out.println("STOP_TAIL_CASCADE_STRUCTURAL_TEST_PASS checks="+checks+" both_sets=true frozen_manifests=true hashes_computed=0");
	}

	/** A real plain SHORT index must reproduce direct-package decisions while the
	 * composite keeps its original LONG policy. Layout ambiguity fails before loading. */
	private static void plainShort(Path root, Fixture full, Fixture positional, byte[] raw, Assignment longBaseline)throws Exception{
		final byte[] query=Arrays.copyOfRange(raw, 1, 100);
		final Fixture plain=fixture(root.resolve("plain_short"), query, 40, SHORT_PIN, .985, true, 28, true);
		final Path child=write(plain.root.resolve("manifest.tsv"), "#bbgap_tail_cascade\t2\n"+plain.settings);
		final String prefix="#bbgap_tail_cascade\t7\n"+full.settings;
		final Path manifest=write(full.root.resolve("plain_child.tsv"), prefix+
			"short40manifest\t"+child+"\nshort40sha80\t"+SHORT_PIN+"\n");
		final StopTailCascadeResources resources=StopTailCascadeResources.readStructural(manifest.toString(), LONG_PIN);
		final Assignment direct=plain.pack.load(1).withEvaluationMarginScale(.985).newWorker().assign("plain_direct", query);
		check(direct.assigned(), "The real plain profile must admit and accept its identical SHORT query");
		same(direct, resources.loadShort40(1).newWorker().assign("plain_child", query));
		same(longBaseline, resources.load(1).newWorker().assign("unchanged_long", raw));
		boolean witness=false;
		try(StopTailCascadeSupport support=resources.loadSupport(1)){
			final StopTailCandidateScan.Scan scan=support.scan("plain_composite", dna(raw), GeneticCode.forTable(11));
			final StopTailCascadeSupport.ContigEvidence evidence=support.classify(scan);
			for(int i=0; i<scan.candidates.size(); i++){
				if(Arrays.equals(scan.candidates.get(i).copyTail(), raw)){
					witness=true;
					same(direct, evidence.familySetAssignment(i, true));
					same(longBaseline, evidence.familySetAssignment(i, false));
				}
			}
		}
		check(witness, "The mixed positional-LONG/plain-SHORT resource must reach the literal scanner query");
		for(String layout:new String[]{"partitions\t4\n", "overlap\t2\n"}){
			final Path bad=write(plain.root.resolve(layout.startsWith("partitions") ? "bad_parts.tsv" : "bad_overlap.tsv"),
				"#bbgap_tail_cascade\t2\n"+plain.settings+layout);
			final Path badRoot=write(full.root.resolve(bad.getFileName()), prefix+"short40manifest\t"+bad+"\nshort40sha80\t"+SHORT_PIN+"\n");
			reject(()->StopTailCascadeResources.readStructural(badRoot.toString(), LONG_PIN), "Plain SHORT caller child must omit");
			final Path missing=write(positional.root.resolve(bad.getFileName()), "#bbgap_tail_cascade\t2\n"+positional.settings.replace(layout, ""));
			final Path missingRoot=write(full.root.resolve("missing_"+bad.getFileName()), prefix+"short40manifest\t"+missing+"\nshort40sha80\t"+SHORT_PIN+"\n");
			reject(()->StopTailCascadeResources.readStructural(missingRoot.toString(), LONG_PIN), "Missing explicit classifier setting");
		}
	}

	/** All61 profiles contain exactly three query register values; no acceptance gate can pass. */
	private static Fixture fixture(Path root, byte[] raw, int min, String pin, double margin)throws Exception{
		return fixture(root, raw, min, pin, margin, false);
	}
	/** Chooses accepted identical references or deliberately rejecting references with span28. */
	private static Fixture fixture(Path root, byte[] raw, int min, String pin, double margin, boolean accept)throws Exception{
		return fixture(root, raw, min, pin, margin, accept, 28);
	}
	/** Builds positional profiles with the requested minimum aligned-query span. */
	private static Fixture fixture(Path root, byte[] raw, int min, String pin, double margin, boolean accept, int minimumSpan)throws Exception{
		return fixture(root, raw, min, pin, margin, accept, minimumSpan, false);
	}
	/** Emit either positional or plain modal rows; existing fixture overloads stay positional. */
	private static Fixture fixture(Path root, byte[] raw, int min, String pin, double margin, boolean accept, int minimumSpan, boolean plain)throws Exception{
		Files.createDirectory(root);
		final ProteinKmerGeometry geometry=ProteinKmerGeometry.reduced(5, GROUPS);
		final char[] row=new char[512];
		if(plain){DDLFamilyProfiles.sketch(raw, geometry, 0, row, null);}
		else{DDLFamilyProfiles.sketchPositional(raw, geometry, 0, 4, 2, row, null);}
		int nonzero=0;
		for(int i=0; i<row.length; i++){if(row[i]!=0 && ++nonzero>3 && !accept){row[i]=0;}}
		check(nonzero>=(accept ? 5 : 3), accept ?
			"Accepted fixtures need at least five distinct registers for their unchanged floor5" :
			"Three distinct registers are required for the independent admission oracle");
		final FamilyAlignmentMetric metric=new FamilyAlignmentMetric(FamilyAlignmentMetric.Mode.PIN_RIGHT_TRIM_LEFT_ALIGNED_QUERY, minimumSpan);
		final byte[] reference=raw.clone();
		if(min==40 && !accept){for(int position=12; position<=84; position+=12){reference[position]=(byte)(reference[position]=='W' ? 'D' : 'W');}}
		final FamilyReferenceCodec.Header header=new FamilyReferenceCodec.Header(FamilyReferenceCodec.ALIGNMENT_SCHEMA, 5, 512, 0, 61, 61, 0,
			FamilyReferenceCodec.FORMULA_VERSION, FamilyReferenceCodec.SPECTRUM_ENCODING, 3.8906, metric);
		final ByteBuilder refs=FamilyReferenceCodec.appendHeader(header, new ByteBuilder());
		final ByteBuilder profiles=new ByteBuilder().append(DDLFamilyProfiles.HEADER16).nl();
		final StringBuilder medians=new StringBuilder("family\tmembers\tlength_p10\tlength_median\tlength_p90\tcount_under_200\n"),
			hydro=new StringBuilder("family_rank\tfamily_id\tN\thydrophobicity_mean\n");
		for(int f=0; f<61; f++){
			final String id="f"+f;
			final float[] spectrum=new float[400];
			spectrum[0]=1;
			final int residues=30*raw.length, words=30*(raw.length-4);
			FamilyReferenceCodec.appendFiveLines(new FamilyReference(id, f, 30, residues, words, new ProteinSequence(id, reference), row, spectrum,
				0, accept ? 30 : 96, 0, 1.1f, accept ? 30 : 110, accept ? 0 : 1.1f, reference, 0, new double[10], FamilyReference.AlignmentCalibration.ENABLED), refs, header);
			profiles.append("tail\t").append(id).append("\t512\t30\t").append(residues).tab().append(words).append("\tnonempty\t");
			for(int field=0; field<3; field++){
				if(field>0){profiles.tab();}
				for(int i=0; i<row.length; i++){
					if(i>0){profiles.append(',');}
					profiles.append(field==0 ? (int)row[i] : row[i]==0 ? 0 : 1);
				}
			}
			profiles.nl();
			medians.append(id).append("\t30\t").append(raw.length).append('\t').append(raw.length).append('\t').append(raw.length).append("\t30\n");
			hydro.append(f).append('\t').append(id).append("\t30\t1.8\n");
		}
		Files.write(root.resolve("references.tsv"), refs.toBytes());
		Files.write(root.resolve("profiles.tsv"), profiles.toBytes());
		write(root.resolve("summary.tsv"), plain ? "namespace\tk\tseed\texponent\tregisterbits\tsketch\tgeometry\n"+
			"tail\t5\t0\t6\t16\tplain\t"+geometry.name()+"\n" :
			"namespace\tk\tseed\texponent\tregisterbits\tsketch\tpartitions\toverlap\ttag_bits\tsegment_words\tgeometry\n"+
			"tail\t5\t0\t6\t16\tpositional-v1\t4\t2\t2\t"+DDLFamilyProfiles.positionalSegmentWords(4, geometry)+"\t"+geometry.name()+"\n");
		write(root.resolve("list.tsv"), root.resolve("profiles.tsv")+"\t"+pin+"\t"+pin+"\n");
		write(root.resolve("medians.tsv"), medians.toString());
		write(root.resolve("hydro.tsv"), hydro.toString());
		final StringBuilder bounds=new StringBuilder("scalar\trows\tmin\tmax\tmin_family\tmax_family\n");
		for(String name:new String[]{"hydrophobicity", "net_charge", "aromaticity", "fraction_small", "pI_tail"}){bounds.append(name).append("\t79565\t-20\t20\tlo\thi\n");}
		write(root.resolve("bounds.tsv"), bounds.toString());
		final StringBuilder settings=new StringBuilder("top\t50\nk\t5\nseed\t0\nwidth\t512\nddlreject\t5\nddlaccept\t1025\n"+
			(plain ? "positional\tf\n" : "positional\tt\npartitions\t4\noverlap\t2\n")+
			"indexearlyreject\tf\ncosine\tranking\nrankingfeatures\tcombined25\n"+
			"scalargates\tt\nscalarz\t3.8906\nidentityacceptfloor\t30\nblosumacceptfloor\t0\nscalarscale\t0.1\n"+
			"alignmentmode\tpin_right_trim_left_aligned_query\nalignmentminspan\t"+minimumSpan+"\naabits\tall:"+GROUPS+"\n");
		for(String[] binding:new String[][]{{"familyreferences", "references.tsv"}, {"ddlprofiles", "list.tsv"}, {"scalarbounds", "bounds.tsv"}}){
			settings.append(binding[0]).append('\t').append(root.resolve(binding[1])).append('\n').append(binding[0]).append("sha80\t")
				.append(binding[0].equals("scalarbounds") ? METADATA_PIN : pin).append('\n');
		}
		final String common=settings.toString(), query="minqueryaa="+min+"\nmaxqueryaa="+raw.length+"\n";
		final Path config=write(root.resolve("classifier.config"), common.replace('\t', '=')+query+
			"heuristic=hitfourth_h020_c150\nlengthmedians="+root.resolve("medians.tsv")+"\nlengthmedianssha80="+METADATA_PIN+
			"\nfamilyhydro="+root.resolve("hydro.tsv")+"\nfamilyhydrosha80="+METADATA_PIN+"\n");
		settings.append("marginscale\t").append(margin).append('\n').append(query.replace('=', '\t'))
			.append("heuristicarm\thitfourth_h020_c150\nfamilymedians\t").append(root.resolve("medians.tsv"))
			.append("\nfamilymedianssha80\t").append(METADATA_PIN).append("\nfamilyhydro\t").append(root.resolve("hydro.tsv"))
			.append("\nfamilyhydrosha80\t").append(METADATA_PIN).append('\n');
		return new Fixture(root, settings.toString(), new FamilyScoringPackage(config, VerificationMode.STRUCTURAL));
	}

	/** Encodes each fixture residue with its literal codon and appends a terminal TAA stop. */
	private static byte[] dna(byte[] protein){
		final StringBuilder result=new StringBuilder();
		for(byte aa:protein){
			final int index=ACIDS.indexOf(aa);
			check(index>=0, "Fixture codon dictionary must cover every residue");
			result.append(CODONS[index]);
		}
		return result.append("TAA").toString().getBytes(StandardCharsets.US_ASCII);
	}
	/** Compares the decision fields and counters asserted by these routing cases, with bit-exact float scores. */
	private static void same(Assignment a, Assignment b){
		check(a!=null && b!=null && a.assigned()==b.assigned() && a.terminal.equals(b.terminal) && a.rawHitFloor==b.rawHitFloor &&
			a.rosterRank==b.rosterRank && (a.familyId==null ? b.familyId==null : a.familyId.equals(b.familyId)) &&
			a.admittedFamilies==b.admittedFamilies && a.dpAttempts==b.dpAttempts && a.ddlMatches==b.ddlMatches &&
			Float.floatToIntBits(a.identity)==Float.floatToIntBits(b.identity) && Float.floatToIntBits(a.blosum)==Float.floatToIntBits(b.blosum),
			"Resource routing must reproduce the independently loaded package's literal decisions and counters");
	}
	/** Single-worker observer counting selected candidates and retaining the last observed alignment/gates. */
	private static final class Counts implements FamilyTailCascade.CandidateObserver{
		/** Counts only candidates selected for alignment. */
		@Override public void rankedCandidate(int f, int hits, int raw, int ranked, boolean chosen){if(chosen){selected++;}}
		/** Retains the most recent alignment scores and effective acceptance gates. */
		@Override public void alignedCandidate(int f, int attempt, float identity, float blosum, double upper, double bUpper, float lo){
			lower=lo;
			identityUpper=upper;
			blosumUpper=bUpper;
			this.identity=identity;
			this.blosum=blosum;
		}
		/** Acceptance notifications do not contribute to this observer's alignment counters. */
		@Override public void accept(int length, String terminal, int ordinal, int hits, float identity, float blosum){}
		/** Number of selected candidates observed so far. */
		int selected;
		/** Last lower gate and last observed identity/BLOSUM scores. */
		float lower, identity, blosum;
		/** Last identity and BLOSUM direct-acceptance upper gates. */
		double identityUpper, blosumUpper;
	}
	/** A synthetic family directory, its manifest settings and the structural package bound to those files. */
	private static final class Fixture{
		/** Retains the supplied path, immutable settings text and package reference without copying files. */
		Fixture(Path root_, String settings_, FamilyScoringPackage pack_){
			root=root_;
			settings=settings_;
			pack=pack_;
		}
		/** Directory containing the synthetic component files. */
		final Path root;
		/** Tab-separated manifest settings assembled by the fixture writer. */
		final String settings;
		/** Structural package used for direct-loading control assignments. */
		final FamilyScoringPackage pack;
	}
	/** Writes ASCII fixture bytes and returns their path for immediate binding. */
	private static Path write(Path path, String text)throws Exception{return Files.write(path, text.getBytes(StandardCharsets.US_ASCII));}
	/** Operation whose invalid fixture inputs are expected to fail validation. */
	private interface Action{
		/** Executes the case; checked failures propagate to the caller. */
		void run()throws Exception;
	}
	/** Requires an argument/state exception with the expected fragment; unrelated exception types propagate. */
	private static void reject(Action action, String fragment)throws Exception{
		try{action.run();}catch(IllegalArgumentException | IllegalStateException e){
			check(e.getMessage()!=null && e.getMessage().contains(fragment), "Wrong rejection: "+e);
			return;
		}
		throw new AssertionError("Expected rejection: "+fragment);
	}
	/** Fails immediately, counting only checks that succeeded. */
	private static void check(boolean valid, String why){
		if(!valid){throw new AssertionError(why);}
		checks++;
	}
	/** Successful checks in this sequential fixture invocation. */
	private static int checks;
	/** Explicit unmeasured LONG/SHORT identities and a deliberately inherited metadata identity. */
	private static final String LONG_PIN="unmeasured", SHORT_PIN="unmeasured", METADATA_PIN="33333333333333333333";
	/** Reduced sketch alphabet and the canonical residue order used by the literal codon dictionary. */
	private static final String GROUPS="A/C/DE/FY/G/H/IV/KR/L/M/N/P/Q/S/T/W", ACIDS="ACDEFGHIKLMNPQRSTVWY";
	/** One codon for each residue in ACIDS, in the same order. */
	private static final String[] CODONS={"GCT", "TGT", "GAT", "GAA", "TTT", "GGT", "CAT", "ATT", "AAA", "CTG", "ATG", "AAT", "CCT", "CAA", "CGT", "TCT", "ACT", "GTT", "TGG", "TAT"};
}
