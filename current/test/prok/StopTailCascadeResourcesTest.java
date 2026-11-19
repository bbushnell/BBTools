package test.prok;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;

import ddl.DDLFamilyProfiles;
import ddl.ProteinKmerGeometry;
import dna.GeneticCode;
import dna.AminoAcid;
import prok.CallGenes;
import prok.StopTailCascadeSupport;
import prok.StopTailCascadeResources;
import prot.DigestSuffix;
import prot.FamilyReference;
import prot.FamilyReferenceCodec;
import prot.FamilyReferenceExport;
import prot.FamilyHeuristicScore;
import prot.FamilyMemberScalarMeans;
import prot.FamilyTailCascade;
import prot.FamilyAlignmentMetric;
import prot.FamilyHitFloorBinding;
import test.prot.FamilyHitFloorIntegrationTest;
import prot.ProteinSequence;
import prot.ProteinSpectrum;
import structures.ByteBuilder;

/** Resource binding fixtures, including a tiny real classifier for caller-margin parity.
 * @author Keqing, Fischl
 */
public final class StopTailCascadeResourcesTest{

	/** Creates a fresh evidence root, checks metadata and live classifier bindings, then writes retained-file inventories. */
	public static void main(String[] args) throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Expected a fresh fixture directory");}
		final Path root=Paths.get(args[0]).toAbsolutePath().normalize();
		Files.createDirectory(root);
		evidenceRoot=root;
		evidence.clear();
		evidenceLinks.clear();
		final Path bundle=Files.createDirectory(root.resolve("bundle"));
		final Path memberDir=Files.createDirectory(bundle.resolve("profiles"));
		final Path reference=write(bundle.resolve("reference.dat"), "reference fixture\n");
		final Path blacklist=write(bundle.resolve("blacklist.faa"), ">word\nACDEF\n");
		final Path bounds=write(bundle.resolve("bounds.tsv"), "bounds fixture\n");
		final Path member=write(memberDir.resolve("profile.tsv"), "profile fixture\n");
		final Path summary=write(memberDir.resolve("summary.tsv"), "summary fixture\n");
		final Path list=write(bundle.resolve("profiles.tsv"), "profiles/profile.tsv\t"+pin(member)+"\t"+pin(summary)+"\n");
		final String text="#bbgap_tail_cascade\t1\n"+
			"familyreferences\treference.dat\nfamilyreferencessha80\t"+pin(reference)+"\n"+
			"ddlprofiles\tprofiles.tsv\nddlprofilessha80\t"+pin(list)+"\n"+
			"blacklist\tblacklist.faa\nblacklistsha80\t"+pin(blacklist)+"\n"+
			"scalarbounds\tbounds.tsv\nscalarboundssha80\t"+pin(bounds)+"\n"+
			"top\t50\nk\t5\nseed\t0\nwidth\t1024\nddlreject\t6\nddlaccept\t13\n"+
			"cosine\tranking\nrankingfeatures\tcombined25\nscalarscale\t0.1\nscalargates\tt\nscalarz\t3.8906\n"+
			"identityacceptfloor\t30\nblosumacceptfloor\t0\npositional\tt\npartitions\t4\noverlap\t2\n";
		final Path manifest=write(bundle.resolve("manifest.tsv"), text);
		final StopTailCascadeResources resources=StopTailCascadeResources.read(manifest.toString());
		check(resources.marginScale==1, "Legacy manifests must preserve exact scale 1");
		check(resources.path("familyreferences").equals(reference.toString()), "Relative reference must resolve against manifest, not cwd");
		check(resources.path("ddlprofiles").equals(list.toString()), "Relative profile list must retain its declared parent");
		check(resources.sha80.equals(pin(manifest)), "Resource identity must describe exact manifest bytes");
		resources.protectOutputs(new String[]{root.resolve("output.gff").toString()});
		reject(()->resources.protectOutputs(new String[]{bundle.resolve("new.gff").toString()}), "inside stop-tail resources");
		final Path alias=root.resolve("member-alias.tsv");
		hardLink(alias, member);
		reject(()->resources.protectOutputs(new String[]{alias.toString()}), "aliases stop-tail resource");
		final Path summaryAlias=root.resolve("summary-alias.tsv");
		hardLink(summaryAlias, summary);
		reject(()->resources.protectOutputs(new String[]{summaryAlias.toString()}), "aliases stop-tail resource");
		final Path absolute=write(bundle.resolve("absolute.tsv"), text.replace("familyreferences\treference.dat", "familyreferences\t"+reference));
		check(StopTailCascadeResources.read(absolute.toString()).path("familyreferences").equals(reference.toString()), "Development absolute paths must remain explicit");
		final Path links=Files.createDirectory(bundle.resolve("links"));
		symbolicLink(links.resolve("profile.tsv"), member);
		final Path linkedSummary=copy(summary, links.resolve("summary.tsv"));
		final Path linkedList=symbolicLink(links.resolve("profiles.tsv"), list);
		final Path linkedManifest=write(bundle.resolve("linked.tsv"), text.replace("ddlprofiles\tprofiles.tsv", "ddlprofiles\tlinks/profiles.tsv"));
		reject(()->StopTailCascadeResources.read(linkedManifest.toString()), "non-symlink input");
		// Member symlinks are supported by the profile loader; its sibling summary uses the declared parent.
		final Path memberList=write(bundle.resolve("member-links.tsv"), "links/profile.tsv\t"+pin(member)+"\t"+pin(linkedSummary)+"\n");
		final Path memberManifest=write(bundle.resolve("member-links-manifest.tsv"), text
			.replace("ddlprofiles\tprofiles.tsv", "ddlprofiles\tmember-links.tsv").replace(pin(list), pin(memberList)));
		final StopTailCascadeResources linked=StopTailCascadeResources.read(memberManifest.toString());
		final Path linkedAlias=hardLink(root.resolve("linked-summary-alias.tsv"), linkedSummary);
		reject(()->linked.protectOutputs(new String[]{linkedAlias.toString()}), "aliases stop-tail resource");
		final Path missing=write(bundle.resolve("missing.tsv"), text.replace("overlap\t2\n", ""));
		reject(()->StopTailCascadeResources.read(missing.toString()), "Missing explicit classifier setting: overlap");
		final Path duplicate=write(bundle.resolve("duplicate.tsv"), text+"top\t50\n");
		reject(()->StopTailCascadeResources.read(duplicate.toString()), "duplicate classifier binding");
		final Path unknown=write(bundle.resolve("unknown.tsv"), text+"surprise\t1\n");
		reject(()->StopTailCascadeResources.read(unknown.toString()), "Unknown");
		final Path wrong=write(bundle.resolve("wrong.tsv"), text.replace("#bbgap_tail_cascade\t1", "#bbgap_tail_cascade\t99"));
		reject(()->StopTailCascadeResources.read(wrong.toString()), "Expected");
		final String v2=text.replace("#bbgap_tail_cascade\t1", "#bbgap_tail_cascade\t2")+"indexearlyreject\tf\n";
		final String blacklistPair="blacklist\tblacklist.faa\nblacklistsha80\t"+pin(blacklist)+"\n";
		final String noBlacklist=v2.replace(blacklistPair, "")+"marginscale\t1.045\n";
		final Path without=write(bundle.resolve("without-blacklist.tsv"), noBlacklist);
		final StopTailCascadeResources unfiltered=StopTailCascadeResources.read(without.toString());
		check(unfiltered.path("blacklist")==null && unfiltered.marginScale==1.045,
			"An omitted v2 blacklist must remain absent without changing the selected margin");
		checkLegacyCallerLocator(without, unfiltered);
		final Path modeOnly=write(bundle.resolve("mode-only.tsv"), noBlacklist+"alignmentmode\tpin_right_trim_left_aligned_query\n");
		reject(()->StopTailCascadeResources.read(modeOnly.toString()), "must be supplied together");
		final Path spanOnly=write(bundle.resolve("span-only.tsv"), noBlacklist+"alignmentminspan\t70\n");
		reject(()->StopTailCascadeResources.read(spanOnly.toString()), "must be supplied together");
		final Path badMode=write(bundle.resolve("bad-mode.tsv"), noBlacklist+"alignmentmode\tunknown\nalignmentminspan\t70\n");
		reject(()->StopTailCascadeResources.read(badMode.toString()), "Unknown caller alignmentmode");
		final Path legacyMetric=write(bundle.resolve("legacy-metric.tsv"), noBlacklist+"alignmentmode\tglocal_full_query\nalignmentminspan\t0\n");
		check(StopTailCascadeResources.read(legacyMetric.toString()).marginScale==1.045, "An explicit legacy metric preserves the selected margin");
		final Path oldMetric=write(bundle.resolve("v1-metric.tsv"), text+"alignmentmode\tglocal_full_query\nalignmentminspan\t0\n");
		reject(()->StopTailCascadeResources.read(oldMetric.toString()), "Unknown");
		reject(()->unfiltered.protectOutputs(new String[]{alias.toString()}), "aliases stop-tail resource");
		reject(()->unfiltered.protectOutputs(new String[]{summaryAlias.toString()}), "aliases stop-tail resource");
		final Path legacyWithout=write(bundle.resolve("legacy-without-blacklist.tsv"), text.replace(blacklistPair, ""));
		reject(()->StopTailCascadeResources.read(legacyWithout.toString()), "Missing explicit classifier setting: blacklist");
		final Path pathOnly=write(bundle.resolve("blacklist-path-only.tsv"), noBlacklist+"blacklist\tblacklist.faa\n");
		reject(()->StopTailCascadeResources.read(pathOnly.toString()), "must be supplied together");
		final Path pinOnly=write(bundle.resolve("blacklist-pin-only.tsv"), noBlacklist+"blacklistsha80\t"+pin(blacklist)+"\n");
		reject(()->StopTailCascadeResources.read(pinOnly.toString()), "must be supplied together");
		final Path emptyList=write(bundle.resolve("empty-profiles.tsv"), "");
		for(int kind=0; kind<2; kind++){
			final String empty=(kind==0 ? text : noBlacklist).replace("ddlprofiles\tprofiles.tsv", "ddlprofiles\tempty-profiles.tsv")
				.replace(pin(list), pin(emptyList));
			final Path emptyManifest=write(bundle.resolve("empty-profile-manifest-"+kind+".tsv"), empty);
			reject(()->StopTailCascadeResources.read(emptyManifest.toString()), "Profile list is empty");
		}
		final Path missingScale=write(bundle.resolve("missing-scale.tsv"), v2);
		reject(()->StopTailCascadeResources.read(missingScale.toString()), "Missing explicit classifier setting: marginscale");
		final Path duplicateScale=write(bundle.resolve("duplicate-scale.tsv"), v2+"marginscale\t1\nmarginscale\t2\n");
		reject(()->StopTailCascadeResources.read(duplicateScale.toString()), "duplicate classifier binding");
		final Path legacyScale=write(bundle.resolve("legacy-scale.tsv"), text+"marginscale\t1\n");
		reject(()->StopTailCascadeResources.read(legacyScale.toString()), "Unknown");
		final Path legacyIndex=write(bundle.resolve("legacy-index.tsv"), text+"indexearlyreject\tf\n");
		reject(()->StopTailCascadeResources.read(legacyIndex.toString()), "Unknown");
		final Path legacyGeometry=write(bundle.resolve("legacy-geometry.tsv"), text+"aabits\thydro\n");
		reject(()->StopTailCascadeResources.read(legacyGeometry.toString()), "Unknown");
		final Path missingIndex=write(bundle.resolve("missing-index.tsv"), v2.replace("indexearlyreject\tf\n", "")+"marginscale\t1\n");
		reject(()->StopTailCascadeResources.read(missingIndex.toString()), "Missing explicit classifier setting: indexearlyreject");
		final Path duplicateIndex=write(bundle.resolve("duplicate-index.tsv"), v2+"indexearlyreject\tt\nmarginscale\t1\n");
		reject(()->StopTailCascadeResources.read(duplicateIndex.toString()), "duplicate classifier binding");
		final Path invalidIndex=write(bundle.resolve("invalid-index.tsv"), v2.replace("indexearlyreject\tf", "indexearlyreject\tmaybe")+"marginscale\t1\n");
		reject(()->StopTailCascadeResources.read(invalidIndex.toString()), "explicit boolean");
		int invalidId=0;
		for(String invalid:new String[]{"-1", "NaN", "Infinity", "-Infinity", "1e309", "bad"}){
			final Path invalidScale=write(bundle.resolve("invalid-scale-"+(invalidId++)+".tsv"), v2+"marginscale\t"+invalid+"\n");
			reject(()->StopTailCascadeResources.read(invalidScale.toString()), "finite nonnegative double");
		}
		for(String value:new String[]{"0", "1", "1.06", "2"}){
			final Path scaled=write(bundle.resolve("scale-"+value+".tsv"), v2+"marginscale\t"+value+"\n");
			check(Double.compare(StopTailCascadeResources.read(scaled.toString()).marginScale, Double.parseDouble(value))==0,
				"Version 2 must retain the exact parsed double");
		}
		final Path badPin=write(bundle.resolve("bad-pin.tsv"), text.replace(pin(list), "00000000000000000000"));
		reject(()->StopTailCascadeResources.read(badPin.toString()), "Profile list digest differs");
		checkCompositeMetadata(bundle, noBlacklist);
		checkFixed100Metadata(bundle, noBlacklist);
		checkIndependentMetadata(bundle, noBlacklist);
		checkRankingMetadata(bundle, noBlacklist);
		// Changing only the manifest must fail before the placeholder reference can be parsed.
		writeBytes(manifest, (text+"# changed\n").getBytes(StandardCharsets.UTF_8));
		reject(()->resources.load(1), "Classifier manifest changed");
		checkRealLoad(Files.createDirectory(root.resolve("real")), true);
		checkRealLoad(Files.createDirectory(root.resolve("real-without-blacklist")), false);
		finish("StopTailCascadeResourcesTest PASS checks="+checks);
	}

	/** Reuse the valid legacy file to check caller dispatch without loading its placeholder models.
	 * Omitted verification and explicit SHA80 retain exact identity and optional-pin checks;
	 * explicit structural loading retains its declared identity and parsed option state.
	 */
	private static void checkLegacyCallerLocator(final Path manifest,
			final StopTailCascadeResources expected)throws Exception{
		assert(Files.isRegularFile(manifest)) : "Legacy locator controls must exercise file dispatch, not compact directories";
		final String manifestPin=pin(manifest);
		final String wrongPin=(manifestPin.charAt(0)=='0' ? "1" : "0")+manifestPin.substring(1);
		for(boolean explicit:new boolean[]{false, true}){
			final CallGenes.TailCascadeOptions options=new CallGenes.TailCascadeOptions();
			if(explicit){options.parse("tailverification", "sha80");}
			final StopTailCascadeResources loaded=options.readResources(manifest.toString());
			check(manifestPin.equals(loaded.sha80), "Legacy caller loading must retain actual manifest identity; explicit="+explicit);
			check(loaded.path("familyreferences").equals(expected.path("familyreferences")) &&
				loaded.marginScale==expected.marginScale, "Legacy file dispatch must preserve resolved payload and margin");
			check(options.verificationSpecified()==explicit && !options.structural() && options.referencePin()==null,
				"Loading a legacy file must not rewrite omitted or explicit SHA80 options");
			options.parse("tailreferencesha80", manifestPin);
			check(manifestPin.equals(options.readResources(manifest.toString()).sha80),
				"A matching optional legacy pin must remain accepted; explicit="+explicit);
			check(options.verificationSpecified()==explicit && !options.structural() && manifestPin.equals(options.referencePin()),
				"Matching-pin loading must preserve the caller's parsed verification state");
			options.parse("tailreferencesha80", wrongPin);
			reject(()->options.readResources(manifest.toString()), "tailreferencesha80 does not match the stop-tail manifest");
		}
		final CallGenes.TailCascadeOptions structural=new CallGenes.TailCascadeOptions();
		structural.parse("tailverification", "structural");
		structural.parse("tailreferencesha80", "unmeasured");
		final StopTailCascadeResources loaded=structural.readResources(manifest.toString());
		check("unmeasured".equals(loaded.sha80), "Structural legacy loading must retain its declared unmeasured identity");
		check(loaded.path("familyreferences").equals(expected.path("familyreferences")) &&
			loaded.marginScale==expected.marginScale, "Structural legacy loading must preserve resolved payload and margin");
		check(structural.verificationSpecified() && structural.structural() && "unmeasured".equals(structural.referencePin()),
			"Legacy structural loading must preserve explicit mode and declared identity");
	}

	/** Legacy early pruning removes repetitive Q200 before alignment. With the explicit
	 * exhaustive index, scale 1 gates 99/1 accept it; scale 2 gates 148/1.5 reject it. */
	private static void checkRealLoad(final Path root, final boolean withBlacklist) throws Exception{
		final byte[] query=new byte[200];
		Arrays.fill(query, (byte)'Q');
		final char[] sketch=new char[1024];
		DDLFamilyProfiles.sketch(query, 5, 0, sketch);
		final ProteinSpectrum.Spectrum counts=new ProteinSpectrum.Spectrum(2);
		counts.add(query);
		final FamilyReference family=new FamilyReference("q", 0, 1, 200, 196,
			new ProteinSequence("q", query), sketch, FamilyReferenceExport.frequencies(counts),
			0.5f, 50, 0.5f, 1.1f, 99, 1);
		final ByteBuilder reference=new ByteBuilder();
		FamilyReferenceCodec.appendHeader(new FamilyReferenceCodec.Header(FamilyReferenceCodec.SCHEMA,
			5, 1024, 0, 1, 1, 0, FamilyReferenceCodec.FORMULA_VERSION, FamilyReferenceCodec.SPECTRUM_ENCODING), reference);
		FamilyReferenceCodec.appendFiveLines(family, reference);
		final Path refs=write(root.resolve("reference.tsv"), reference.toString());
		final byte[] original=Files.readAllBytes(refs);
		final Path blacklist=write(root.resolve("blacklist.faa"), ">unused\nACDEF\n");
		final Path bounds=write(root.resolve("bounds.tsv"), "scalar\trows\tmin\tmax\tmin_family\tmax_family\n"+
			"hydrophobicity\t79565\t-20\t20\tlo\thi\nnet_charge\t79565\t-20\t20\tlo\thi\n"+
			"aromaticity\t79565\t-20\t20\tlo\thi\nfraction_small\t79565\t-20\t20\tlo\thi\npI_tail\t79565\t-20\t20\tlo\thi\n");
		final char[] positional=new char[1024];
		DDLFamilyProfiles.sketchPositional(query, 5, 0, 4, 2, positional);
		final ByteBuilder profile=new ByteBuilder().append(DDLFamilyProfiles.HEADER16).nl()
			.append("tail\tq\t1024\t1\t200\t196\tnonempty\t");
		for(int kind=0; kind<3; kind++){
			if(kind>0){profile.tab();}
			for(int i=0; i<positional.length; i++){
				if(i>0){profile.append(',');}
				profile.append(kind==0 ? (int)positional[i] : positional[i]==0 ? 0 : 1);
			}
		}
		final Path data=write(root.resolve("profile.tsv"), profile.nl().toString());
		final Path summary=write(root.resolve("summary.tsv"),
			"namespace\tk\tseed\texponent\tregisterbits\tsketch\tpartitions\toverlap\ttag_bits\tsegment_words"+
			(withBlacklist ? "\tblacklist_sha80\tblacklist_words" : "")+"\n"+
			"tail\t5\t0\t6\t16\t"+DDLFamilyProfiles.SKETCH_POSITIONAL+"\t4\t2\t"+
			DDLFamilyProfiles.positionalTagBits(4)+"\t"+DDLFamilyProfiles.positionalSegmentWords(4)+
			(withBlacklist ? "\t"+pin(blacklist)+"\t1" : "")+"\n");
		final Path profiles=write(root.resolve("profiles.tsv"), data+"\t"+pin(data)+"\t"+pin(summary)+"\n");
		final String settings="familyreferences\t"+refs+"\nfamilyreferencessha80\t"+pin(refs)+"\n"+
			"ddlprofiles\t"+profiles+"\nddlprofilessha80\t"+pin(profiles)+"\n"+
			(withBlacklist ? "blacklist\t"+blacklist+"\nblacklistsha80\t"+pin(blacklist)+"\n" : "")+
			"scalarbounds\t"+bounds+"\nscalarboundssha80\t"+pin(bounds)+"\n"+
			"top\t1\nk\t5\nseed\t0\nwidth\t1024\nddlreject\t0\nddlaccept\t1025\n"+
			"cosine\tranking\nrankingfeatures\tcombined25\nscalarscale\t0.1\nscalargates\tf\nscalarz\t0\n"+
			"identityacceptfloor\t30\nblosumacceptfloor\t0\npositional\tt\npartitions\t4\noverlap\t2\n";
		final HashMap<String, String> options=new HashMap<String, String>();
		for(String line:settings.split("\n")){
			final String[] fields=line.split("\t");
			options.put(fields[0], fields[1]);
		}
		final FamilyTailCascade.LoadedClassifier legacyDirect=FamilyTailCascade.loadCallable(new HashMap<String, String>(options), 1);
		final FamilyTailCascade.Assignment legacyExpected=legacyDirect.newWorker().assign("q", query);
		check(legacyDirect.usesEarlyRejection() && !legacyExpected.assigned() && legacyExpected.dpAttempts==0,
			"Legacy early-index policy must retain the observed pre-alignment rejection of Q200");
		final String legacyHeader=withBlacklist ? "#bbgap_tail_cascade\t1\n" : "#bbgap_tail_cascade\t2\nindexearlyreject\tt\nmarginscale\t1\n";
		final Path legacy=write(root.resolve("legacy.tsv"), legacyHeader+settings);
		final FamilyTailCascade.Assignment actual=StopTailCascadeResources.read(legacy.toString()).load(1).newWorker().assign("q", query);
		check(Arrays.equals(decision(legacyExpected), decision(actual)), "Legacy caller decision bytes must match direct loading without a margin view");
		writeBytes(root.resolve("legacy-decision.tsv"), decision(actual));
		options.put("indexearlyreject", "f");
		final FamilyTailCascade.LoadedClassifier direct=FamilyTailCascade.loadCallable(options, 1);
		final FamilyTailCascade.Assignment expected=direct.newWorker().assign("q", query);
		writeBytes(root.resolve("direct-decision.tsv"), decision(expected));
		check(!direct.usesEarlyRejection() && expected.assigned() && expected.identity==100 && expected.dpAttempts==1,
			"With early pruning disabled, identical Q200 must pass identity 100>99 after one alignment; actual="+
			new String(decision(expected), StandardCharsets.US_ASCII));
		for(double scale:new double[]{0, 1, 2}){
			final Path manifest=write(root.resolve("scale-"+scale+".tsv"), "#bbgap_tail_cascade\t2\n"+settings+
				"indexearlyreject\tf\nmarginscale\t"+scale+"\n");
			final FamilyTailCascade.LoadedClassifier classifier=StopTailCascadeResources.read(manifest.toString()).load(2);
			check(!classifier.usesEarlyRejection(), "Version 2 must carry the explicit scientific index policy into the loaded classifier");
			final FamilyTailCascade.Assignment loaded=classifier.newWorker().assign("q", query);
			final FamilyTailCascade.Assignment oracle=direct.withEvaluationMarginScale(scale).newWorker().assign("q", query);
			check(loaded.assigned()==(scale<2), "Literal scale-2 gates must reject the identical query; scales 0/1 must accept");
			check(Arrays.equals(decision(loaded), decision(oracle)), "Caller resource must preserve every direct scaled decision field");
			writeBytes(root.resolve("decision-"+scale+".tsv"), decision(loaded));
		}
		final Path early=write(root.resolve("early.tsv"), "#bbgap_tail_cascade\t2\n"+settings+"indexearlyreject\tt\nmarginscale\t1\n");
		final FamilyTailCascade.LoadedClassifier earlyClassifier=StopTailCascadeResources.read(early.toString()).load(1);
		check(earlyClassifier.usesEarlyRejection() && Arrays.equals(decision(legacyExpected), decision(earlyClassifier.newWorker().assign("q", query))),
			"Version 2 with explicit early rejection and scale 1 must retain legacy decision bytes");
		check(Arrays.equals(original, Files.readAllBytes(refs)), "Alternate caller margins must never rewrite reference bytes");
		final Path small=write(root.resolve("bin100.tsv"), "#bbgap_tail_cascade\t2\n"+settings+
			"indexearlyreject\tf\nmarginscale\t2\nminqueryaa\t100\nmaxqueryaa\t124\n");
		final Path medium=write(root.resolve("bin150.tsv"), "#bbgap_tail_cascade\t2\n"+settings+
			"indexearlyreject\tf\nmarginscale\t1\nminqueryaa\t150\nmaxqueryaa\t199\n");
		final Path composite=write(root.resolve("composite.tsv"), "#bbgap_tail_cascade\t3\n"+settings+
			"indexearlyreject\tf\nmarginscale\t1\nshort100manifest\tbin100.tsv\nshort100sha80\t"+pin(small)+
			"\nshort150manifest\tbin150.tsv\nshort150sha80\t"+pin(medium)+"\n");
		final StopTailCascadeResources combined=StopTailCascadeResources.read(composite.toString());
		try(StopTailCascadeSupport support=combined.loadSupport(2)){
			final StringBuilder dna=new StringBuilder();
			for(int n:new int[]{100, 150, 200}){
				for(int i=0; i<n; i++){dna.append("CAA");}
				dna.append("TAA");
			}
			final StopTailCascadeSupport.ContigEvidence e=support.classify(support.scan("composite", dna.toString().getBytes(StandardCharsets.US_ASCII), GeneticCode.forTable(11)));
			check(e.size()==3 && !e.assignment(0).assigned() && e.assignment(1).assigned() && e.assignment(2).assigned(),
				"The short100 scale2 must reject while independent short150/full200 scale1 views accept");
			check(Arrays.equals(decision(e.assignment(2)), decision(expected)), "Adding short children must preserve every full200 assignment field");
		}
		final Path independent=write(root.resolve("bin100-149.tsv"), "#bbgap_tail_cascade\t2\n"+settings+
			"indexearlyreject\tf\nmarginscale\t2\nminqueryaa\t100\nmaxqueryaa\t149\n");
		final Path v4=write(root.resolve("independent.tsv"), "#bbgap_tail_cascade\t4\n"+settings+
			"indexearlyreject\tf\nmarginscale\t1\nshort100manifest\tbin100-149.tsv\nshort100sha80\t"+pin(independent)+
			"\nshort150manifest\tbin150.tsv\nshort150sha80\t"+pin(medium)+"\n");
		try(StopTailCascadeSupport support=StopTailCascadeResources.read(v4.toString()).loadSupport(2)){
			final StringBuilder dna=new StringBuilder();
			for(int n:new int[]{100, 125, 149, 150, 200}){
				for(int i=0; i<n; i++){dna.append("CAA");}
				dna.append("TAA");
			}
			final StopTailCascadeSupport.ContigEvidence evidence=support.classify(support.scan("v4", dna.toString().getBytes(StandardCharsets.US_ASCII), GeneticCode.forTable(11)));
			check(evidence.size()==5 && evidence.emitsLength(125) && evidence.emitsLength(149), "Resource-driven independent mode emits the previously disabled bin");
			for(int i=0; i<5; i++){check(evidence.assignment(i).assigned()==(i>=3), "Independent child margin must not affect150+ support");}
			check(Arrays.equals(decision(evidence.assignment(4)), decision(expected)), "Version4 preserves full200 decision bytes");
		}
		checkHeuristicLoad(root, settings, direct, query, medium);
		checkCap150Load(root, settings, direct, medium);
		checkCap100Load(root, settings, direct);
		checkFittedLoad(root);
		checkFractionalLoad(root, settings, family, query);
		checkMetricLoad(root, settings, family, query, direct);
		writeBytes(small, "# changed child\n".getBytes(StandardCharsets.US_ASCII), java.nio.file.StandardOpenOption.APPEND);
		reject(()->combined.loadSupport(1), "Classifier manifest changed");
	}

	/** The fitted short child uses the same saved parent and vector as evaluation,
	 * while the existing long classifier and scan policy stay byte-for-byte equal. */
	private static void checkFittedLoad(final Path root)throws Exception{
		final Path dir=Files.createDirectory(root.resolve("fitted"));
		final byte[] query=FamilyHitFloorIntegrationTest.sequence();
		query[0]='M';
		FamilyHitFloorIntegrationTest.fixture(dir.resolve("parent"), query);
		// Exact outputs of the shared fixture helper; no recursive discovery.
		for(String name:new String[]{"references.tsv", "profiles.tsv", "summary.tsv", "list.tsv", "bounds.tsv",
			"medians.tsv", "hydro.tsv", "classifier.config"}){retain(dir.resolve("parent").resolve(name));}
		final Path parent=dir.resolve("parent/classifier.config");
		final Path cohorts=write(dir.resolve("cohorts.tsv"), "synthetic fitting provenance\n");
		final FamilyHitFloorBinding binding=new FamilyHitFloorBinding(parent, pin(parent), 1, cohorts, pin(cohorts), 1);
		final long[] initial=new long[21];
		initial[12]=1;
		final Path floors=dir.resolve("floors.tsv");
		binding.write(floors, binding.probe().withInitialSupport(1, initial), 1, initial);
		retain(floors);
		final Path compressed=dir.resolve("floors.tsv.gz");
		try(java.util.zip.GZIPOutputStream out=new java.util.zip.GZIPOutputStream(Files.newOutputStream(compressed))){
			out.write(Files.readAllBytes(floors));
		}
		retain(compressed);
		final FamilyTailCascade.LoadedClassifier oracle=binding.load(floors, pin(floors));
		check("f1".equals(oracle.newWorker().assign("literal", query).familyId) &&
			!binding.parent().newWorker().assign("six", query).assigned(),
			"Literal five-hit family passes its fitted floor and fails fixed six");
		final String settings="#bbgap_tail_fitted\t1\nparent\t"+parent+"\nparentsha80\t"+pin(parent)+
			"\ncohorts\t"+cohorts+"\ncohortssha80\t"+pin(cohorts)+"\nmarginscale\t1\n";
		final Path oldChild=root.resolve("cap100-tiny.tsv");
		final String baseline=new String(Files.readAllBytes(root.resolve("cap100-true.tsv")), StandardCharsets.US_ASCII);
		final String oldPair="short40manifest\t"+oldChild+"\nshort40sha80\t"+pin(oldChild)+"\n";
		final byte[] coding=AminoAcid.toNTs(query);
		final ByteBuilder bases=new ByteBuilder().append(coding).append("TAAATG");
		for(int i=1; i<100; i++){bases.append("CAA");}
		bases.append("TAA");
		for(Path table:new Path[]{floors, compressed}){
			final Path child=write(dir.resolve(table.getFileName()+".child"), settings+
				"floors\t"+table+"\nfloorssha80\t"+pin(table)+"\n");
			final Path composite=write(dir.resolve(table.getFileName()+".composite"), baseline.replace(oldPair,
				"short40manifest\t"+child+"\nshort40sha80\t"+pin(child)+"\n"));
			final StopTailCascadeResources resources=StopTailCascadeResources.read(composite.toString());
			reject(()->StopTailCascadeResources.read(child.toString()), "independent40..99 child slot");
			for(int threads:new int[]{1, 2}){
				try(StopTailCascadeSupport support=resources.loadSupport(threads);
					StopTailCascadeSupport control=StopTailCascadeResources.read(root.resolve("cap100-true.tsv").toString()).loadSupport(threads)){
					final prok.StopTailCandidateScan.Scan scan=support.scan("fitted", bases.toBytes(), GeneticCode.forTable(11));
					final StopTailCascadeSupport.ContigEvidence actual=support.classify(scan);
					final StopTailCascadeSupport.ContigEvidence expected=control.classify(control.scan("fitted", bases.toBytes(), GeneticCode.forTable(11)));
					check(actual.size()==expected.size(), "Changing short admission must preserve scanner candidates");
					int literal=0;
					for(int i=0; i<actual.size(); i++){
						final byte[] tail=scan.candidates.get(i).copyTail();
						final byte[] tiny=Arrays.copyOfRange(tail, Math.max(0, tail.length-99), tail.length);
						final FamilyTailCascade.Assignment fitted=actual.familySetAssignment(i, true);
						check(Arrays.equals(decision(fitted), decision(oracle.newWorker().assign("direct", tiny))),
							"Caller short assignment must equal ordinary saved-vector loading for every scanned query");
						if(Arrays.equals(tail, query)){
							literal++;
							check(fitted.assigned() && fitted.familyHitFloorAdmission && fitted.rawHitFloor==1 && fitted.rawHitFloorMax==6,
								"The actual scanner reaches the heterogeneous fitted policy, including unsupported-family fallback");
						}
						if(tail.length==100){
							check(Arrays.equals(decision(actual.familySetAssignment(i, false)), decision(expected.familySetAssignment(i, false))),
								"The unchanged long root must retain every assignment field");
						}
					}
					check(literal==1, "Exactly one literal64-aa synthetic query must reach the fitted child");
				}
			}
			for(Path input:new Path[]{parent, cohorts, table, dir.resolve("parent/profiles.tsv"), dir.resolve("parent/summary.tsv"), dir.resolve("parent/medians.tsv")}){
				final Path alias=hardLink(root.resolve(table.getFileName()+"-alias-"+input.getFileName()), input);
				reject(()->resources.protectOutputs(new String[]{alias.toString()}), "aliases stop-tail resource");
			}
			final String childText=new String(Files.readAllBytes(child), StandardCharsets.US_ASCII);
			final String[] invalid={childText.replace("parentsha80\t"+pin(parent), "parentsha80\t00000000000000000000"),
				childText.replace("floorssha80\t"+pin(table), "floorssha80\t00000000000000000000"),
				childText.replace("cohorts\t"+cohorts+"\n", ""), childText+"marginscale\t1\n", childText+"top\t49\n"};
			for(int n=0; n<invalid.length; n++){
				final Path badChild=write(dir.resolve(table.getFileName()+"-invalid-"+n), invalid[n]);
				final Path badRoot=write(dir.resolve(table.getFileName()+"-invalid-root-"+n), baseline.replace(oldPair,
					"short40manifest\t"+badChild+"\nshort40sha80\t"+pin(badChild)+"\n"));
				reject(()->StopTailCascadeResources.read(badRoot.toString()), n<2 ? "differs from its pin" : n==2 ? "Missing fitted-child" : "duplicate fitted-child");
			}
			final Path wrongMargin=write(dir.resolve(table.getFileName()+"-wrong-margin"), childText.replace("marginscale\t1\n", "marginscale\t2\n"));
			final Path wrongRoot=write(dir.resolve(table.getFileName()+"-wrong-root"), baseline.replace(oldPair,
				"short40manifest\t"+wrongMargin+"\nshort40sha80\t"+pin(wrongMargin)+"\n"));
			reject(()->StopTailCascadeResources.read(wrongRoot.toString()).loadSupport(1), "Floor-table provenance differs");
		}
	}

	/** A fractional long root must reach the shared kernel without changing its exact-K5 child. */
	private static void checkFractionalLoad(final Path root, final String settings,
			final FamilyReference family, final byte[] query) throws Exception{
		assert(query.length==200 && family.members==1) : "This fixture derives one Q200 family from checkRealLoad";
		final Path directory=Files.createDirectory(root.resolve("fractional"));
		final ProteinKmerGeometry geometry=ProteinKmerGeometry.make(4, "hydro");
		final char[] plain=new char[1024], positional=new char[1024];
		final long words=DDLFamilyProfiles.sketch(query, geometry, 0, plain, null);
		DDLFamilyProfiles.sketchPositional(query, geometry, 0, 4, 2, positional, null);
		final FamilyReference fractional=new FamilyReference("q", 0, 1, query.length, words,
			family.consensus, plain, family.spectrum, 0.5f, 50, 0.5f, 1.1f, 99, 1);
		final ByteBuilder bytes=new ByteBuilder();
		FamilyReferenceCodec.appendHeader(new FamilyReferenceCodec.Header(FamilyReferenceCodec.SCHEMA,
			4, 1024, 0, 1, 1, 0, FamilyReferenceCodec.FORMULA_VERSION, FamilyReferenceCodec.SPECTRUM_ENCODING), bytes);
		FamilyReferenceCodec.appendFiveLines(fractional, bytes);
		final Path references=write(directory.resolve("references.tsv"), bytes.toString());
		bytes.clear().append(DDLFamilyProfiles.HEADER16).nl().append("tail\tq\t1024\t1\t200\t").append(words).append("\tnonempty\t");
		for(int kind=0; kind<3; kind++){
			if(kind>0){bytes.tab();}
			for(int i=0; i<positional.length; i++){
				if(i>0){bytes.append(',');}
				bytes.append(kind==0 ? (int)positional[i] : positional[i]==0 ? 0 : 1);
			}
		}
		final Path profile=write(directory.resolve("profile.tsv"), bytes.nl().toString());
		final boolean blacklisted=settings.contains("blacklist\t");
		final Path summary=write(directory.resolve("summary.tsv"),
			"namespace\tk\tseed\texponent\tregisterbits\tsketch\tpartitions\toverlap\ttag_bits\tsegment_words\tgeometry"+
			(blacklisted ? "\tblacklist_sha80\tblacklist_words" : "")+"\n"+
			"tail\t4\t0\t6\t16\t"+DDLFamilyProfiles.SKETCH_POSITIONAL+"\t4\t2\t"+
			DDLFamilyProfiles.positionalTagBits(4)+"\t"+DDLFamilyProfiles.positionalSegmentWords(4, geometry)+"\t"+geometry.name()+
			(blacklisted ? "\t"+pin(root.resolve("blacklist.faa"))+"\t1" : "")+"\n");
		final Path profiles=write(directory.resolve("profiles.tsv"), profile+"\t"+pin(profile)+"\t"+pin(summary)+"\n");
		final String longSettings=settings.replace("k\t5\n", "k\t4\n")
			.replace("familyreferences\t"+root.resolve("reference.tsv"), "familyreferences\t"+references)
			.replace("familyreferencessha80\t"+pin(root.resolve("reference.tsv")), "familyreferencessha80\t"+pin(references))
			.replace("ddlprofiles\t"+root.resolve("profiles.tsv"), "ddlprofiles\t"+profiles)
			.replace("ddlprofilessha80\t"+pin(root.resolve("profiles.tsv")), "ddlprofilessha80\t"+pin(profiles))+
			"indexearlyreject\tf\n";
		final Path medians=root.resolve("real-medians.tsv"), hydro=root.resolve("real-hydro.tsv");
		final Path tiny=root.resolve("cap100-tiny.tsv");
		final byte[] tinyBefore=Files.readAllBytes(tiny);
		final String policy="#bbgap_tail_cascade\t7\n"+longSettings+
			"marginscale\t1\nminqueryaa\t100\nmaxqueryaa\t100\nheuristicarm\thitfourth_h020_c150\n"+
			"familymedians\t"+medians+"\nfamilymedianssha80\t"+pin(medians)+
			"\nfamilyhydro\t"+hydro+"\nfamilyhydrosha80\t"+pin(hydro)+
			"\nshort40manifest\t"+tiny+"\nshort40sha80\t"+pin(tiny)+"\n";
		final HashMap<String, String> options=new HashMap<String, String>();
		for(String line:longSettings.split("\n")){
			final String[] fields=line.split("\t");
			options.put(fields[0], fields[1]);
		}
		options.put("aabits", "hydro");
		final FamilyTailCascade.LoadedClassifier direct=FamilyTailCascade.loadCallable(options, 1);
		final FamilyTailCascade.LoadedClassifier oracle=direct.withHeuristicRanking(new FamilyHeuristicScore(direct,
			medians, pin(medians), hydro, pin(hydro), 100, 100, "hitfourth_h020_c150"));
		final StringBuilder dna=new StringBuilder("ATG");
		for(int i=1; i<100; i++){dna.append("CAA");}
		dna.append("TAA");
		final byte[] bases=dna.toString().getBytes(StandardCharsets.US_ASCII);
		final GeneticCode code=GeneticCode.forTable(11);
		for(String key:new String[]{"aabits", "aapartition"}){
			final Path manifest=write(directory.resolve(key+".tsv"), policy+key+"\thydro\n");
			final StopTailCascadeResources resources=StopTailCascadeResources.read(manifest.toString());
			for(int threads:new int[]{1, 2}){
				try(StopTailCascadeSupport mixed=resources.loadSupport(threads);
					StopTailCascadeSupport baseline=StopTailCascadeResources.read(root.resolve("cap100-true.tsv").toString()).loadSupport(threads)){
					final prok.StopTailCandidateScan.Scan scan=mixed.scan("fractional", bases, code);
					final StopTailCascadeSupport.ContigEvidence actual=mixed.classify(scan);
					final StopTailCascadeSupport.ContigEvidence expected=baseline.classify(baseline.scan("fractional", bases, code));
					check(actual.size()==1 && expected.size()==1, "One literal100-aa ORF must reach both resource compositions");
					check(Arrays.equals(decision(actual.familySetAssignment(0, false)),
						decision(oracle.newWorker().assign("oracle", scan.candidates.get(0).copyTail()))),
						"The fractional caller root must match direct geometry-aware classification for "+key);
					check(Arrays.equals(decision(actual.familySetAssignment(0, true)), decision(expected.familySetAssignment(0, true))),
						"Fractional long geometry must not leak into the unchanged exact-K5 short40 child");
				}
			}
		}
		final Path omitted=write(directory.resolve("omitted.tsv"), policy);
		// Exact K4 cannot consume the retained pentamer blacklist; that guard precedes profile loading.
		reject(()->StopTailCascadeResources.read(omitted.toString()).load(1), blacklisted ?
			"Prefix protein blacklist length 5 exceeds raw geometry span 4" : "geometry mismatch");
		final Path wrong=write(directory.resolve("wrong.tsv"), policy+"aabits\tblosum\n");
		reject(()->StopTailCascadeResources.read(wrong.toString()).load(1), "geometry mismatch");
		final Path ambiguous=write(directory.resolve("ambiguous.tsv"), policy+"aabits\thydro\naapartition\thydro\n");
		reject(()->StopTailCascadeResources.read(ambiguous.toString()).load(1), "either aabits or aapartition");
		check(Arrays.equals(tinyBefore, Files.readAllBytes(tiny)), "Fractional root loading must preserve the exact child manifest bytes");
	}

	/** Actual version7 loading binds exact100 and the carried independent short search only. */
	private static void checkCap100Load(final Path root, final String settings,
			final FamilyTailCascade.LoadedClassifier direct)throws Exception{
		final Path medians=root.resolve("real-medians.tsv"), hydro=root.resolve("real-hydro.tsv");
		final String ranking="heuristicarm\thitfourth_h020_c150\nfamilymedians\t"+medians+
			"\nfamilymedianssha80\t"+pin(medians)+"\nfamilyhydro\t"+hydro+"\nfamilyhydrosha80\t"+pin(hydro)+"\n";
		final String base=settings+"indexearlyreject\tf\nmarginscale\t1\n"+ranking;
		final Path tinyData=copy(root.resolve("profile.tsv"), root.resolve("cap100-tiny-profile.tsv"));
		final Path tinyProfiles=write(root.resolve("cap100-tiny-profiles.tsv"), tinyData+"\t"+pin(tinyData)+"\t"+pin(root.resolve("summary.tsv"))+"\n");
		final String tinySettings=base.replace("ddlprofiles\t"+root.resolve("profiles.tsv"), "ddlprofiles\t"+tinyProfiles)
			.replace("ddlprofilessha80\t"+pin(root.resolve("profiles.tsv")), "ddlprofilessha80\t"+pin(tinyProfiles));
		final Path tiny=write(root.resolve("cap100-tiny.tsv"), "#bbgap_tail_cascade\t2\n"+tinySettings+"minqueryaa\t40\nmaxqueryaa\t99\n");
		final Path lower=write(root.resolve("cap100-lower.tsv"), "#bbgap_tail_cascade\t2\n"+
			base.replace("marginscale\t1\n", "marginscale\t2\n")+"minqueryaa\t100\nmaxqueryaa\t149\n");
		final Path medium=write(root.resolve("cap100-medium.tsv"), "#bbgap_tail_cascade\t2\n"+
			base+"minqueryaa\t150\nmaxqueryaa\t199\n");
		final String tinyBinding="short40manifest\t"+tiny+"\nshort40sha80\t"+pin(tiny)+"\n";
		final String lowerBinding="short100manifest\t"+lower+"\nshort100sha80\t"+pin(lower)+"\n";
		final String policy="#bbgap_tail_cascade\t7\n"+base+"minqueryaa\t100\nmaxqueryaa\t100\n"+tinyBinding;
		final FamilyTailCascade.LoadedClassifier exact=direct.withHeuristicRanking(new FamilyHeuristicScore(direct,
			medians, pin(medians), hydro, pin(hydro), 100, 100, "hitfourth_h020_c150"));
		final FamilyTailCascade.LoadedClassifier shortSet=direct.withHeuristicRanking(new FamilyHeuristicScore(direct,
			medians, pin(medians), hydro, pin(hydro), 40, 99, "hitfourth_h020_c150"));
		final int[] lengths={40, 99, 100, 101, 150};
		final StringBuilder dna=new StringBuilder();
		for(int length:lengths){
			dna.append("ATG");
			for(int i=1; i<length; i++){dna.append("CAA");}
			dna.append("TAA");
		}
		final byte[] bases=dna.toString().getBytes(StandardCharsets.US_ASCII);
		for(boolean first:new boolean[]{false, true}){
			final Path manifest=write(root.resolve("cap100-"+first+".tsv"), policy+(first ? "" : "firststart\tf\n"));
			final StopTailCascadeResources resources=StopTailCascadeResources.read(manifest.toString());
			check(resources.minQueryAa==100 && resources.maxQueryAa==100 && resources.firstStart==first && resources.allFamilySets &&
				resources.queryLengths().equals("40-99,100"), "Version7 preserves exact100 gates and explicit composite query policy");
			for(int threads:new int[]{1, 2}){try(StopTailCascadeSupport support=resources.loadSupport(threads)){
				final prok.StopTailCandidateScan.Scan scan=support.scan("cap100", bases, GeneticCode.forTable(11));
				final StopTailCascadeSupport.ContigEvidence evidence=support.classify(scan);
				check(evidence.maximumTailLength==100 && evidence.firstStart==first && evidence.searchesAllFamilySets(),
					"The actual resource loader reaches capped100 scanner and two-set evidence");
				for(int i=0; i<scan.candidates.size(); i++){
					final byte[] tail=scan.candidates.get(i).copyTail();
					final byte[] shortTail=Arrays.copyOfRange(tail, Math.max(0, tail.length-99), tail.length);
					final FamilyTailCascade.Assignment shortResult=shortSet.newWorker().assign("oracle", shortTail);
					check(Arrays.equals(decision(shortResult), decision(evidence.familySetAssignment(i, true))),
						"Version7 short search matches the independent99-cap direct kernel");
					FamilyTailCascade.Assignment expected=shortResult;
					if(tail.length==100){
						final FamilyTailCascade.Assignment exactResult=exact.newWorker().assign("oracle", tail);
						check(Arrays.equals(decision(exactResult), decision(evidence.familySetAssignment(i, false))),
							"Version7 dispatch uses the independent direct exact100 view");
						if(exactResult.assigned() || !shortResult.assigned()){expected=exactResult;}
					}
					check(Arrays.equals(decision(expected), decision(evidence.assignment(i))), "Version7 preserves one long-first binary decision");
				}
			}}
		}
		final String[] invalid={policy.replace("minqueryaa\t100\nmaxqueryaa\t100\n", ""),
			policy.replace("minqueryaa\t100\nmaxqueryaa\t100\n", "minqueryaa\t150\nmaxqueryaa\t150\n"),
			policy.replace(tinyBinding, ""), policy+"searchall\tf\n",
			policy+lowerBinding, policy+"short150manifest\t"+medium+"\nshort150sha80\t"+pin(medium)+"\n",
			policy.replace("#bbgap_tail_cascade\t7", "#bbgap_tail_cascade\t6")};
		for(int i=0; i<invalid.length; i++){
			final Path path=write(root.resolve("cap100-invalid-"+i+".tsv"), invalid[i]);
			reject(()->StopTailCascadeResources.read(path.toString()), i<2 || i==6 ? "Unsupported caller query interval" :
				"Version7 requires");
		}
		final Path wrongPin=write(root.resolve("cap100-wrong-child-pin.tsv"),
			policy.replace("short40sha80\t"+pin(tiny), "short40sha80\t00000000000000000000"));
		reject(()->StopTailCascadeResources.read(wrongPin.toString()), "manifest digest mismatch");
		final Path forbiddenChild=write(root.resolve("cap100-nested.tsv"), policy.replace(tinyBinding,
			"short40manifest\t"+root.resolve("cap100-true.tsv")+"\nshort40sha80\t"+pin(root.resolve("cap100-true.tsv"))+"\n"));
		reject(()->StopTailCascadeResources.read(forbiddenChild.toString()), "must be version 2");
	}

	/** Actual version6 loading must bind both capped long views and preserve the separate short search. */
	private static void checkCap150Load(final Path root, final String settings,
			final FamilyTailCascade.LoadedClassifier direct, final Path medium)throws Exception{
		final Path medians=root.resolve("real-medians.tsv"), hydro=root.resolve("real-hydro.tsv");
		final String ranking="heuristicarm\thitfourth_h020_c150\nfamilymedians\t"+medians+
			"\nfamilymedianssha80\t"+pin(medians)+"\nfamilyhydro\t"+hydro+"\nfamilyhydrosha80\t"+pin(hydro)+"\n";
		final String base=settings+"indexearlyreject\tf\nmarginscale\t1\n"+ranking;
		final Path tinyData=copy(root.resolve("profile.tsv"), root.resolve("cap150-tiny-profile.tsv"));
		final Path tinyProfiles=write(root.resolve("cap150-tiny-profiles.tsv"), tinyData+"\t"+pin(tinyData)+"\t"+pin(root.resolve("summary.tsv"))+"\n");
		final String tinySettings=base.replace("ddlprofiles\t"+root.resolve("profiles.tsv"), "ddlprofiles\t"+tinyProfiles)
			.replace("ddlprofilessha80\t"+pin(root.resolve("profiles.tsv")), "ddlprofilessha80\t"+pin(tinyProfiles));
		final Path tiny=write(root.resolve("cap150-tiny.tsv"), "#bbgap_tail_cascade\t2\n"+tinySettings+"minqueryaa\t40\nmaxqueryaa\t99\n");
		final Path lower=write(root.resolve("cap150-lower.tsv"), "#bbgap_tail_cascade\t2\n"+
			base.replace("marginscale\t1\n", "marginscale\t2\n")+"minqueryaa\t100\nmaxqueryaa\t149\n");
		final String tinyBinding="short40manifest\t"+tiny+"\nshort40sha80\t"+pin(tiny)+"\n";
		final String lowerBinding="short100manifest\t"+lower+"\nshort100sha80\t"+pin(lower)+"\n";
		final String policy="#bbgap_tail_cascade\t6\n"+base+"minqueryaa\t150\nmaxqueryaa\t150\n"+tinyBinding+lowerBinding;
		final FamilyTailCascade.LoadedClassifier exact=direct.withHeuristicRanking(new FamilyHeuristicScore(direct,
			medians, pin(medians), hydro, pin(hydro), 150, 150, "hitfourth_h020_c150"));
		final FamilyTailCascade.LoadedClassifier small=direct.withEvaluationMarginScale(2).withHeuristicRanking(new FamilyHeuristicScore(direct,
			medians, pin(medians), hydro, pin(hydro), 100, 149, "hitfourth_h020_c150"));
		final FamilyTailCascade.LoadedClassifier shortSet=direct.withHeuristicRanking(new FamilyHeuristicScore(direct,
			medians, pin(medians), hydro, pin(hydro), 40, 99, "hitfourth_h020_c150"));
		final int[] lengths={40, 99, 100, 149, 150, 151, 200};
		final StringBuilder dna=new StringBuilder();
		for(int length:lengths){
			dna.append("ATG");
			for(int i=1; i<length; i++){dna.append("CAA");}
			dna.append("TAA");
		}
		final byte[] bases=dna.toString().getBytes(StandardCharsets.US_ASCII);
		for(boolean first:new boolean[]{false, true}){
			final Path manifest=write(root.resolve("cap150-"+first+".tsv"), policy+(first ? "" : "firststart\tf\n"));
			final StopTailCascadeResources resources=StopTailCascadeResources.read(manifest.toString());
			check(resources.minQueryAa==150 && resources.maxQueryAa==150 && resources.firstStart==first && resources.allFamilySets &&
				resources.queryLengths().equals("40-99,100-149,150"), "Version6 preserves exact150 gates and explicit composite query policy");
			for(int threads:new int[]{1, 2}){try(StopTailCascadeSupport support=resources.loadSupport(threads)){
				final prok.StopTailCandidateScan.Scan scan=support.scan("cap150", bases, GeneticCode.forTable(11));
				final StopTailCascadeSupport.ContigEvidence evidence=support.classify(scan);
				check(evidence.maximumTailLength==150 && evidence.firstStart==first && evidence.searchesAllFamilySets(),
					"The actual resource loader reaches capped150 scanner and two-set evidence");
				for(int i=0; i<scan.candidates.size(); i++){
					final byte[] tail=scan.candidates.get(i).copyTail();
					final byte[] shortTail=Arrays.copyOfRange(tail, Math.max(0, tail.length-99), tail.length);
					final FamilyTailCascade.Assignment shortResult=shortSet.newWorker().assign("oracle", shortTail);
					check(Arrays.equals(decision(shortResult), decision(evidence.familySetAssignment(i, true))),
						"Version6 short search matches the independent99-cap direct kernel");
					FamilyTailCascade.Assignment expected=shortResult;
					if(tail.length>=100){
						final FamilyTailCascade.Assignment longResult=(tail.length==150 ? exact : small).newWorker().assign("oracle", tail);
						check(Arrays.equals(decision(longResult), decision(evidence.familySetAssignment(i, false))),
							"Version6 dispatch uses independent direct lower/exact150 views with their own margins");
						if(longResult.assigned() || !shortResult.assigned()){expected=longResult;}
					}
					check(Arrays.equals(decision(expected), decision(evidence.assignment(i))), "Version6 preserves one long-first binary decision");
				}
			}}
		}
		final String[] invalid={policy.replace("minqueryaa\t150\nmaxqueryaa\t150\n", ""),
			policy.replace("minqueryaa\t150\nmaxqueryaa\t150\n", "minqueryaa\t200\nmaxqueryaa\t200\n"),
			policy.replace(lowerBinding, ""), policy.replace(tinyBinding, ""), policy+"searchall\tf\n",
			policy+"short150manifest\t"+medium+"\nshort150sha80\t"+pin(medium)+"\n",
			policy.replace("#bbgap_tail_cascade\t6", "#bbgap_tail_cascade\t5")};
		for(int i=0; i<invalid.length; i++){
			final Path path=write(root.resolve("cap150-invalid-"+i+".tsv"), invalid[i]);
			reject(()->StopTailCascadeResources.read(path.toString()), i<2 ? "Unsupported caller query interval" : i==6 ? "Unknown" : "Version6 requires");
		}
		final Path wrongPin=write(root.resolve("cap150-wrong-child-pin.tsv"), policy.replace("short100sha80\t"+pin(lower), "short100sha80\t00000000000000000000"));
		reject(()->StopTailCascadeResources.read(wrongPin.toString()), "manifest digest mismatch");
		final Path forbiddenChild=write(root.resolve("cap150-nested.tsv"), policy.replace(lowerBinding,
			"short100manifest\t"+root.resolve("cap150-true.tsv")+"\nshort100sha80\t"+pin(root.resolve("cap150-true.tsv"))+"\n"));
		reject(()->StopTailCascadeResources.read(forbiddenChild.toString()), "must be version 2");
	}

	/** The manifest must reach the real bound scorer; copied margin/ranking views preserve geometry. */
	private static void checkMetricLoad(Path root, String settings, FamilyReference family, byte[] query,
			FamilyTailCascade.LoadedClassifier legacy)throws Exception{
		final FamilyAlignmentMetric metric=new FamilyAlignmentMetric(FamilyAlignmentMetric.Mode.PIN_RIGHT_TRIM_LEFT_ALIGNED_QUERY, 70);
		final FamilyReferenceCodec.Header header=new FamilyReferenceCodec.Header(FamilyReferenceCodec.ALIGNMENT_SCHEMA,
			5, 1024, 0, 1, 1, 0, FamilyReferenceCodec.FORMULA_VERSION, FamilyReferenceCodec.SPECTRUM_ENCODING, 0, metric);
		final FamilyReference bound=new FamilyReference(family.id, family.rank, family.members, family.residues, family.validWords,
			family.consensus, family.ddl, family.spectrum, family.cosineLower, family.identityLower, family.blosumLower,
			family.cosineAccept, family.identityUpper, family.blosumFinal, family.consensusRaw, family.scalarMask, family.scalarBounds,
			FamilyReference.AlignmentCalibration.ENABLED);
		final ByteBuilder bytes=new ByteBuilder();
		FamilyReferenceCodec.appendHeader(header, bytes);
		FamilyReferenceCodec.appendFiveLines(bound, bytes, header);
		final Path refs=write(root.resolve("metric-reference.tsv"), bytes.toString());
		final Path medians=root.resolve("real-medians.tsv"), hydro=root.resolve("real-hydro.tsv");
		final String policy="#bbgap_tail_cascade\t2\n"+settings
			.replace("familyreferences\t"+root.resolve("reference.tsv"), "familyreferences\t"+refs)
			.replace("familyreferencessha80\t"+pin(root.resolve("reference.tsv")), "familyreferencessha80\t"+pin(refs))+
			"indexearlyreject\tf\nmarginscale\t1\nalignmentmode\tpin_right_trim_left_aligned_query\nalignmentminspan\t70\n"+
			"heuristicarm\thitfourth_h020_c150\nfamilymedians\t"+medians+"\nfamilymedianssha80\t"+pin(medians)+
			"\nfamilyhydro\t"+hydro+"\nfamilyhydrosha80\t"+pin(hydro)+"\n";
		final Path manifest=write(root.resolve("metric-manifest.tsv"), policy);
		final FamilyTailCascade.LoadedClassifier loaded=StopTailCascadeResources.read(manifest.toString()).load(2);
		check(loaded.usesHeuristicRanking(), "New metric manifest retains the explicitly calibrated common pivot");
		// Blosum62 gives R/Q=+1: retaining R40 scores840, above the clipped Q160 score800.
		// That positive mismatch prefix must stay in the optimum and fail the 99/1 gates.
		final byte[] positivePrefix=query.clone();
		Arrays.fill(positivePrefix, 0, 40, (byte)'R');
		final FamilyTailCascade.Assignment retained=loaded.newWorker().assign("positive-prefix", positivePrefix);
		check(!retained.assigned() && retained.dpAttempts==1 && "blosum_reject".equals(retained.terminal),
			"Positive R/Q substitutions must not be clipped merely to improve identity; actual="+
			new String(decision(retained), StandardCharsets.US_ASCII));
		// W/Q=-2 makes the unwanted prefix costly; clipping leaves160 exact Q/Q pairs.
		final byte[] extended=query.clone();
		Arrays.fill(extended, 0, 40, (byte)'W');
		check(!legacy.newWorker().assign("prefix", extended).assigned(), "Legacy full-query geometry penalizes the literal forty-W prefix");
		final FamilyTailCascade.Assignment clipped=loaded.newWorker().assign("prefix", extended);
		check(clipped.assigned() && clipped.identity==100 && clipped.blosum==1,
			"Manifest loading, margin view and heuristic view retain clipped alignment and aligned-query denominator; actual="+
			new String(decision(clipped), StandardCharsets.US_ASCII));
		final Path mismatch=write(root.resolve("metric-span-mismatch.tsv"), policy.replace("alignmentminspan\t70", "alignmentminspan\t28"));
		boolean failed=false;
		try{StopTailCascadeResources.read(mismatch.toString()).load(1);}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "Manifest and reference header cannot silently disagree about alignment span");
		final Path missing=write(root.resolve("metric-missing-binding.tsv"), policy.replace("alignmentmode\tpin_right_trim_left_aligned_query\nalignmentminspan\t70\n", ""));
		failed=false;
		try{StopTailCascadeResources.read(missing.toString()).load(1);}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "A new metric header cannot masquerade as an absent legacy manifest binding");
		// Synthetic gates are bound afresh to custom scores; this does not migrate real calibrated thresholds.
		final prot.ProteinSubstitutionMatrix custom=prot.ProteinSubstitutionMatrix.defaultMatrix();
		final FamilyAlignmentMetric customMetric=new FamilyAlignmentMetric(metric.mode, 70, custom);
		final FamilyReferenceCodec.Header customHeader=new FamilyReferenceCodec.Header(FamilyReferenceCodec.MATRIX_SCHEMA,
			5, 1024, 0, 1, 1, 0, FamilyReferenceCodec.FORMULA_VERSION, FamilyReferenceCodec.SPECTRUM_ENCODING, 0, customMetric);
		final ByteBuilder customBytes=new ByteBuilder();
		FamilyReferenceCodec.appendHeader(customHeader, customBytes);
		FamilyReferenceCodec.appendFiveLines(bound, customBytes, customHeader);
		final Path customRefs=write(root.resolve("custom-metric-reference.tsv"), customBytes.toString());
		final String customBase=policy.replace("familyreferences\t"+refs, "familyreferences\t"+customRefs)
			.replace("familyreferencessha80\t"+pin(refs), "familyreferencessha80\t"+pin(customRefs));
		final String portable="matrix\tcustom40\nmatrixbinding\t"+custom.binding()+"\n";
		final Path customManifest=write(root.resolve("custom-metric-manifest.tsv"), customBase+portable);
		checkCallGenesRawScoreI0Guard(root, customManifest);
		final FamilyTailCascade.LoadedClassifier customLoaded=StopTailCascadeResources.read(customManifest.toString()).load(2);
		check(customLoaded.alignmentMatrix().same(custom), "Ordinary caller manifest must retain exact calibrated custom40 binding");
		final HashMap<String, String> customOptions=new HashMap<String, String>();
		for(String line:(customBase+portable).split("\n")){
			if(line.startsWith("#") || line.startsWith("marginscale\t")){continue;}
			final String[] fields=line.split("\t", 2);
			customOptions.put(fields[0], fields[1]);
		}
		customOptions.remove("heuristicarm");
		customOptions.remove("familymedians");
		customOptions.remove("familymedianssha80");
		customOptions.remove("familyhydro");
		customOptions.remove("familyhydrosha80");
		final FamilyTailCascade.LoadedClassifier customDirect=FamilyTailCascade.loadCallable(customOptions, 1);
		final FamilyTailCascade.Assignment callerDecision=customLoaded.newWorker().assign("custom", query);
		final FamilyTailCascade.Assignment directDecision=customDirect.newWorker().assign("custom", query);
		check(callerDecision.assigned() && callerDecision.identity==100 && callerDecision.blosum==1 && callerDecision.dpAttempts==1,
			"Identical custom-bound Q200 must reach one alignment and pass with identity100 and normalized score1");
		check(callerDecision.identity==directDecision.identity && callerDecision.blosum==directDecision.blosum &&
			callerDecision.dpAttempts==directDecision.dpAttempts && callerDecision.assigned()==directDecision.assigned(),
			"Caller custom alignment evidence must match direct calibrated scoring");
		final Path noScores=write(root.resolve("custom-no-scores.tsv"), customBase+"matrix\tcustom40\n");
		reject(()->StopTailCascadeResources.read(noScores.toString()), "portable binding");
		final Path wrongScores=write(root.resolve("custom-wrong-scores.tsv"), customBase+"matrix\tblosum62\n");
		reject(()->StopTailCascadeResources.read(wrongScores.toString()).load(1), "differs from the requested callable package");
		final Path orphan=write(root.resolve("custom-orphan-binding.tsv"), customBase+"matrixbinding\t"+custom.binding()+"\n");
		reject(()->StopTailCascadeResources.read(orphan.toString()), "explicit matrix selector");
	}

	/** Actual CallGenes startup must bind omitted raw i0 and reject explicit mismatches unless overridden. */
	private static void checkCallGenesRawScoreI0Guard(Path root, Path customManifest)throws Exception{
		final Path fasta=write(root.resolve("raw-score-i0.fna"), ">raw-score-i0\nATGCAATAA\n");
		final Path out=root.resolveSibling(root.getFileName()+"-raw-score-i0.gff");
		final String[] base={"in="+fasta, "out="+out, "taxonomy=f", "cdsonly=t", "passes=1", "t=1", "overwrite=t",
			"annotation=tailcascade", "tailverification=sha80", "tailreferencesha80="+pin(customManifest), "tailreference="+customManifest,
			"tailbonusmode=formula", "tailform=BLOSUMratio", "tailscale=0.45", "tailcap=400"};
		new CallGenes(base);
		final String[] explicit=Arrays.copyOf(base, base.length+1);
		explicit[base.length]="taili0=1";
		rejectCallGenesI0(()->new CallGenes(explicit), "matrix scale=4 formula i0=1.0");
		final String[] explicitFalse=Arrays.copyOf(base, base.length+2);
		explicitFalse[base.length]="taili0=1";
		explicitFalse[base.length+1]="tailrawscorei0override=f";
		rejectCallGenesI0(()->new CallGenes(explicitFalse), "matrix scale=4 formula i0=1.0");
		final String[] override=Arrays.copyOf(base, base.length+2);
		override[base.length]="taili0=1";
		override[base.length+1]="tailrawscorei0override=t";
		new CallGenes(override);
	}

	/** Exercise optional bindings through the actual V5 loader and dispatcher, not just metadata parsing. */
	private static void checkHeuristicLoad(final Path root, final String settings,
			final FamilyTailCascade.LoadedClassifier direct, final byte[] query, final Path medium) throws Exception{
		final Path medians=write(root.resolve("real-medians.tsv"),
			"family\tmembers\tlength_p10\tlength_median\tlength_p90\tcount_under_200\nq\t1\t200\t200\t200\t0\n");
		final Path hydro=write(root.resolve("real-hydro.tsv"),
			"family_rank\tfamily_id\tN\thydrophobicity_mean\n0\tq\t1\t-3.5\n");
		final String[] names={"hydrophobicity", "net_charge", "aromaticity", "fraction_small", "pI_tail"};
		final double[] values={-3.5, 0, 0, 0, 5.5};
		final StringBuilder moments=new StringBuilder("family_rank\tfamily_id\tN");
		for(String name:names){for(String metric:new String[]{"mean", "sd", "min", "max"}){moments.append('\t').append(name).append('_').append(metric);}}
		moments.append("\n0\tq\t1");
		for(double value:values){moments.append('\t').append(value).append("\t0\t").append(value).append('\t').append(value);}
		final Path means=write(root.resolve("real-means.tsv"), moments.append('\n').toString());
		final String policy="#bbgap_tail_cascade\t2\n"+settings+"indexearlyreject\tf\nmarginscale\t1\n";
		final int[] lengths={40, 99, 100, 149, 150, 200};
		final StringBuilder dna=new StringBuilder();
		for(int length:lengths){
			for(int i=0; i<length; i++){dna.append("CAA");}
			dna.append("TAA");
		}
		final byte[] bases=dna.toString().getBytes(StandardCharsets.US_ASCII);
		for(int mode=0; mode<4; mode++){
			final boolean ranked=mode==1 || mode==3, memberMeans=mode>=2;
			final String arm=mode==3 ? "hitfourth_h020" : "hitfourth";
			final String extra=(ranked ? "heuristicarm\t"+arm+"\nfamilymedians\t"+medians+"\nfamilymedianssha80\t"+pin(medians)+"\n" : "")+
				(mode==3 ? "familyhydro\t"+hydro+"\nfamilyhydrosha80\t"+pin(hydro)+"\n" : "")+
				(memberMeans ? "familymeans\t"+means+"\nfamilymeanssha80\t"+pin(means)+"\n" : "");
			final Path tiny=write(root.resolve("real40-"+mode+".tsv"), policy+"minqueryaa\t40\nmaxqueryaa\t99\n"+extra);
			final Path short100=write(root.resolve("real100-"+mode+".tsv"), policy+"minqueryaa\t100\nmaxqueryaa\t149\n"+extra);
			final Path parent=write(root.resolve("real-v5-"+mode+".tsv"), policy.replace("#bbgap_tail_cascade\t2\n", "#bbgap_tail_cascade\t5\n")+
				"firststart\tf\nsearchall\tf\nshort40manifest\t"+tiny+"\nshort40sha80\t"+pin(tiny)+"\nshort100manifest\t"+short100+
				"\nshort100sha80\t"+pin(short100)+"\nshort150manifest\t"+medium+"\nshort150sha80\t"+pin(medium)+"\n");
			final StopTailCascadeResources resources=StopTailCascadeResources.read(parent.toString());
			check(!resources.load(1).usesHeuristicRanking(), "Child-only bindings must leave the full200 classifier ordinary");
			final byte[][] expected=new byte[lengths.length][];
			for(int i=0; i<lengths.length; i++){
				FamilyTailCascade.LoadedClassifier oracle=direct;
				if(lengths[i]<150){
					if(memberMeans){oracle=oracle.withMemberScalarMeans(new FamilyMemberScalarMeans(oracle, means, pin(means)));}
					if(ranked){oracle=oracle.withHeuristicRanking(new FamilyHeuristicScore(oracle, medians, pin(medians),
						mode==3 ? hydro : null, mode==3 ? pin(hydro) : null,
						lengths[i]<100 ? 40 : 100, lengths[i]<100 ? 99 : 149, arm));}
				}
				expected[i]=decision(oracle.newWorker().assign("q", Arrays.copyOf(query, lengths[i])));
			}
			for(int threads:new int[]{1, 2}){
				try(StopTailCascadeSupport support=resources.loadSupport(threads)){
					final StopTailCascadeSupport.ContigEvidence evidence=support.classify(support.scan("heuristic", bases, GeneticCode.forTable(11)));
					check(evidence.size()==lengths.length, "V5 must dispatch all six boundary lengths exactly once");
					for(int i=0; i<lengths.length; i++){
						check(Arrays.equals(expected[i], decision(evidence.assignment(i))),
							"Manifest and direct immutable views must agree at length="+lengths[i]+" mode="+mode+" threads="+threads);
					}
				}
			}
		}
	}

	/** Optional payloads must be explicit, versioned and protected even when stored outside the manifest directory. */
	private static void checkRankingMetadata(final Path bundle, final String v2)throws Exception{
		final Path extra=Files.createDirectory(bundle.getParent().resolve("ranking-inputs"));
		final Path medians=write(extra.resolve("medians.tsv"), "median metadata fixture\n");
		final Path hydro=write(extra.resolve("hydro.tsv"), "KD metadata fixture\n");
		final Path means=write(extra.resolve("means.tsv"), "scalar metadata fixture\n");
		final String medianPair="familymedians\t"+medians+"\nfamilymedianssha80\t"+pin(medians)+"\n";
		final String hydroPair="familyhydro\t"+hydro+"\nfamilyhydrosha80\t"+pin(hydro)+"\n";
		final String meansPair="familymeans\t"+means+"\nfamilymeanssha80\t"+pin(means)+"\n";
		final Path plain=write(bundle.resolve("rank-plain.tsv"), v2);
		final StopTailCascadeResources ordinary=StopTailCascadeResources.read(plain.toString());
		check(ordinary.heuristicArm==null && ordinary.path("familymeans")==null && ordinary.path("familymedians")==null,
			"Absent optional keys preserve the ordinary caller binding without invented payloads");
		final Path meanManifest=write(bundle.resolve("rank-means-only.tsv"), v2+meansPair);
		final StopTailCascadeResources meanResources=StopTailCascadeResources.read(meanManifest.toString());
		check(meanResources.heuristicArm==null && meanResources.path("familymeans").equals(means.toString()),
			"External means are an independent explicit option, not an implicit heuristic selection");
		final Path alias=hardLink(bundle.getParent().resolve("means-alias.tsv"), means);
		reject(()->meanResources.protectOutputs(new String[]{alias.toString()}), "aliases stop-tail resource");
		for(String key:new String[]{"familymeans", "familymedians", "familyhydro"}){
			final Path pathOnly=write(bundle.resolve(key+"-path-only.tsv"), v2+key+"\t"+means+"\n");
			reject(()->StopTailCascadeResources.read(pathOnly.toString()), "paired path and sha80");
			final Path pinOnly=write(bundle.resolve(key+"-pin-only.tsv"), v2+key+"sha80\t"+pin(means)+"\n");
			reject(()->StopTailCascadeResources.read(pinOnly.toString()), "paired path and sha80");
		}
		final Path wrongPin=write(bundle.resolve("means-wrong-pin.tsv"), v2+meansPair.replace(pin(means), "00000000000000000000"));
		reject(()->StopTailCascadeResources.read(wrongPin.toString()), "differs from its pin");
		final Path link=symbolicLink(extra.resolve("means-link.tsv"), means);
		final Path linked=write(bundle.resolve("means-linked.tsv"), v2+meansPair.replace(means.toString(), link.toString()));
		reject(()->StopTailCascadeResources.read(linked.toString()), "symlink");
		final String legacy=v2.replace("#bbgap_tail_cascade\t2", "#bbgap_tail_cascade\t1")
			.replace("indexearlyreject\tf\n", "").replace("marginscale\t1.045\n", "");
		final Path legacyMeans=write(bundle.resolve("legacy-means.tsv"), legacy+meansPair);
		reject(()->StopTailCascadeResources.read(legacyMeans.toString()), "Unknown");
		final Path unused=write(bundle.resolve("unused-medians.tsv"), v2+medianPair);
		reject(()->StopTailCascadeResources.read(unused.toString()), "explicit heuristicarm");
		final Path badArm=write(bundle.resolve("bad-arm.tsv"), v2+"heuristicarm\tsumfourth\n"+medianPair);
		reject(()->StopTailCascadeResources.read(badArm.toString()), "Unsupported caller heuristicarm");
		final Path missingKd=write(bundle.resolve("missing-kd.tsv"), v2+"heuristicarm\thitfourth_h020\n"+medianPair);
		reject(()->StopTailCascadeResources.read(missingKd.toString()), "H arms require pinned member KD");
		final Path unusedKd=write(bundle.resolve("unused-kd.tsv"), v2+"heuristicarm\thitfourth\n"+medianPair+hydroPair);
		reject(()->StopTailCascadeResources.read(unusedKd.toString()), "H0 must omit KD");
		final Path fullHeuristic=write(bundle.resolve("full-heuristic.tsv"), v2+"heuristicarm\thitfourth\n"+medianPair);
		check("hitfourth".equals(StopTailCascadeResources.read(fullHeuristic.toString()).heuristicArm),
			"The full200 view must accept explicit ranking metadata now that every calibrated caller bin supports ranking");
		final Path medium=write(bundle.resolve("rank-medium.tsv"), v2+"minqueryaa\t150\nmaxqueryaa\t199\n");
		final String parent=v2.replace("#bbgap_tail_cascade\t2", "#bbgap_tail_cascade\t4")+
			"short150manifest\t"+medium+"\nshort150sha80\t"+pin(medium)+"\n";
		for(String arm:new String[]{"hitfourth", "hitfourth_h020", "hitfourth_h020_c150"}){
			final Path child=write(bundle.resolve("rank-child-"+arm+".tsv"), v2+"minqueryaa\t100\nmaxqueryaa\t149\n"+
				"heuristicarm\t"+arm+"\n"+medianPair+(arm.equals("hitfourth") ? "" : hydroPair)+meansPair);
			final Path manifest=write(bundle.resolve("rank-composite-"+arm+".tsv"), parent+
				"short100manifest\t"+child+"\nshort100sha80\t"+pin(child)+"\n");
			final StopTailCascadeResources composite=StopTailCascadeResources.read(manifest.toString());
			check(composite.heuristicArm==null && composite.queryLengths().equals("100-149,150-199,200"),
				"An explicit child heuristic must not change the parent policy or scanner intervals");
			final Path medianAlias=hardLink(bundle.getParent().resolve("median-alias-"+arm+".tsv"), medians);
			reject(()->composite.protectOutputs(new String[]{medianAlias.toString()}), "aliases stop-tail resource");
		}
	}

	/** Checks the separate fixed100 loader and rejects intervals/children unsupported by that entry point. */
	private static void checkFixed100Metadata(final Path bundle, final String v2) throws Exception{
		final Path fixed=write(bundle.resolve("fixed100.tsv"), v2+"minqueryaa\t100\nmaxqueryaa\t100\n");
		final StopTailCascadeResources r=StopTailCascadeResources.readFixed100(fixed.toString());
		check(r.queryLengths().equals("100") && !r.hasShortBins() && r.marginScale==1.045, "Fixed100 binds its own explicit length and margin");
		reject(()->r.loadSupport(1), "Caller support requires");
		reject(()->StopTailCascadeResources.read(fixed.toString()), "Unsupported caller query interval");
		for(String range:new String[]{"", "minqueryaa\t100\nmaxqueryaa\t124\n", "minqueryaa\t200\nmaxqueryaa\t200\n"}){
			final Path wrong=write(bundle.resolve("fixed100-wrong-"+range.length()+".tsv"), v2+range);
			reject(()->StopTailCascadeResources.readFixed100(wrong.toString()), "Unsupported caller query interval");
		}
		final Path composite=write(bundle.resolve("fixed100-child.tsv"), v2.replace("#bbgap_tail_cascade\t2", "#bbgap_tail_cascade\t3")+
			"short100manifest\t"+fixed+"\nshort100sha80\t"+pin(fixed)+"\n");
		reject(()->StopTailCascadeResources.read(composite.toString()), "Unsupported caller query interval");
		reject(()->StopTailCascadeResources.readFixed100(composite.toString()), "version-2");
	}

	/** Checks version4's paired independent short bins and rejects legacy, missing or mismatched bindings. */
	private static void checkIndependentMetadata(final Path bundle, final String v2) throws Exception{
		final Path small=write(bundle.resolve("independent100.tsv"), v2+"minqueryaa\t100\nmaxqueryaa\t149\n");
		final Path medium=write(bundle.resolve("independent150.tsv"), v2+"minqueryaa\t150\nmaxqueryaa\t199\n");
		final String base=v2.replace("#bbgap_tail_cascade\t2", "#bbgap_tail_cascade\t4");
		final String smallBinding="short100manifest\t"+small+"\nshort100sha80\t"+pin(small)+"\n";
		final String mediumBinding="short150manifest\t"+medium+"\nshort150sha80\t"+pin(medium)+"\n";
		final Path manifest=write(bundle.resolve("independent-composite.tsv"), base+smallBinding+mediumBinding);
		final StopTailCascadeResources resources=StopTailCascadeResources.read(manifest.toString());
		check(resources.queryLengths().equals("100-149,150-199,200"), "Version4 reports its complete explicit dispatch policy");
		reject(()->StopTailCascadeResources.read(small.toString()), "Unsupported caller query interval");
		final Path old=write(bundle.resolve("independent-in-v3.tsv"), (base+smallBinding+mediumBinding).replace("#bbgap_tail_cascade\t4", "#bbgap_tail_cascade\t3"));
		reject(()->StopTailCascadeResources.read(old.toString()), "Unsupported caller query interval");
		final Path missing=write(bundle.resolve("independent-missing150.tsv"), base+smallBinding);
		reject(()->StopTailCascadeResources.read(missing.toString()), "requires both");
		final Path oldChild=write(bundle.resolve("independent-old-child.tsv"), v2+"minqueryaa\t100\nmaxqueryaa\t124\n");
		final Path wrong=write(bundle.resolve("independent-wrong-child.tsv"), base+mediumBinding+
			"short100manifest\t"+oldChild+"\nshort100sha80\t"+pin(oldChild)+"\n");
		reject(()->StopTailCascadeResources.read(wrong.toString()), "Unsupported caller query interval");
	}

	/** Checks legacy composite child resolution and output protection without parsing placeholder models. */
	private static void checkCompositeMetadata(final Path bundle, final String v2) throws Exception{
		final Path children=Files.createDirectory(bundle.getParent().resolve("children"));
		final String absolute=v2.replace("familyreferences\treference.dat", "familyreferences\t"+bundle.resolve("reference.dat"))
			.replace("ddlprofiles\tprofiles.tsv", "ddlprofiles\t"+bundle.resolve("profiles.tsv"))
			.replace("scalarbounds\tbounds.tsv", "scalarbounds\t"+bundle.resolve("bounds.tsv"));
		final String smallText=absolute+"minqueryaa\t100\nmaxqueryaa\t124\n";
		final Path small=write(children.resolve("100.tsv"), smallText);
		final Path medium=write(children.resolve("150.tsv"), absolute+"minqueryaa\t150\nmaxqueryaa\t199\n");
		final String base=v2.replace("#bbgap_tail_cascade\t2", "#bbgap_tail_cascade\t3");
		final String bindings="short100manifest\t"+small+"\nshort100sha80\t"+pin(small)+
			"\nshort150manifest\t"+medium+"\nshort150sha80\t"+pin(medium)+"\n";
		final Path manifest=write(bundle.resolve("composite-metadata.tsv"), base+bindings);
		final StopTailCascadeResources r=StopTailCascadeResources.read(manifest.toString());
		check(r.hasShortBins() && r.queryLengths().equals("100-124,150-199,200") && r.marginScale==1.045,
			"Composite must explicitly bind both allowed short intervals and retain the independent parent margin");
		reject(()->r.protectOutputs(new String[]{children.resolve("would-overwrite.gff").toString()}), "inside stop-tail resources");
		reject(()->StopTailCascadeResources.read(small.toString()), "must be bound as a child");
		final Path none=write(bundle.resolve("composite-empty.tsv"), base);
		reject(()->StopTailCascadeResources.read(none.toString()), "at least one");
		final Path unpaired=write(bundle.resolve("composite-unpaired.tsv"), base+"short100manifest\t"+small+"\n");
		reject(()->StopTailCascadeResources.read(unpaired.toString()), "must be supplied together");
		final Path bad=write(bundle.resolve("composite-bad-pin.tsv"), base+bindings.replace(pin(small), "00000000000000000000"));
		reject(()->StopTailCascadeResources.read(bad.toString()), "manifest digest mismatch");
		final Path swapped=write(bundle.resolve("composite-swapped.tsv"), base+"short100manifest\t"+medium+"\nshort100sha80\t"+pin(medium)+"\n");
		reject(()->StopTailCascadeResources.read(swapped.toString()), "different query-length bin");
		final Path unknown=write(bundle.resolve("composite-125.tsv"), base+bindings+"short125manifest\t"+small+"\n");
		reject(()->StopTailCascadeResources.read(unknown.toString()), "Unknown");
		final Path nested=write(children.resolve("nested.tsv"), absolute.replace("#bbgap_tail_cascade\t2", "#bbgap_tail_cascade\t3"));
		final Path cyclic=write(bundle.resolve("composite-nested.tsv"), base+"short100manifest\t"+nested+"\nshort100sha80\t"+pin(nested)+"\n");
		reject(()->StopTailCascadeResources.read(cyclic.toString()), "must be version 2");
		final Path noBounds=write(children.resolve("unbounded.tsv"), absolute);
		final Path unbounded=write(bundle.resolve("composite-unbounded.tsv"), base+"short100manifest\t"+noBounds+"\nshort100sha80\t"+pin(noBounds)+"\n");
		reject(()->StopTailCascadeResources.read(unbounded.toString()), "explicit minqueryaa and maxqueryaa");
	}

	/** Serializes the selected decision fields used by these exact byte-parity oracles. */
	private static byte[] decision(final FamilyTailCascade.Assignment a){
		assert(a!=null) : "A completed fixture classification must have an explicit decision, including rejection";
		return (a.assigned()+"\t"+a.familyId+"\t"+a.terminal+"\t"+a.rosterRank+"\t"+a.ddlMatches+"\t"+
			a.ddlPosition+"\t"+a.dpAttempt+"\t"+a.dpAttempts+"\t"+a.score+"\t"+a.identity+"\t"+a.blosum+"\t"+
			a.scoreMeasured+"\t"+a.identityMeasured+"\t"+a.blosumMeasured+"\n").getBytes(StandardCharsets.US_ASCII);
	}

	/** Operation expected to reject invalid resource bindings. */
	private interface Action{
		/** Executes the case, allowing unrelated checked failures to propagate. */
		void run() throws Exception;
	}
	/** Requires CallGenes' classifier wrapper and the specific nested raw-score intercept diagnostic. */
	private static void rejectCallGenesI0(Action action, String expected) throws Exception{
		try{action.run();}catch(IllegalArgumentException e){
			check(e.getMessage()!=null && e.getMessage().contains("Could not load stop-tail classifier"),
				"Expected CallGenes classifier wrapper: "+e);
			for(Throwable cause=e.getCause(); cause!=null; cause=cause.getCause()){
				if(cause.getMessage()!=null && cause.getMessage().contains(expected)){return;}
			}
			throw e;
		}
		throw new AssertionError("Expected CallGenes raw-score i0 diagnostic: "+expected);
	}
	/** Requires an argument/state failure with the expected fragment; other failures propagate. */
	private static void reject(Action action, String expected) throws Exception{
		boolean rejected=false;
		try{action.run();}catch(IllegalArgumentException | IllegalStateException e){
			rejected=e.getMessage()!=null && e.getMessage().contains(expected);
			if(!rejected){throw e;}
		}
		check(rejected, "Expected diagnostic: "+expected);
	}
	/** Writes UTF-8 fixture text and records its path in the retained-evidence inventory. */
	private static Path write(Path path, String text) throws Exception{return writeBytes(path, text.getBytes(StandardCharsets.UTF_8));}
	/** Keep every known creation in a deterministic inventory without scanning directories. */
	private static String relative(Path path){
		final Path absolute=path.toAbsolutePath().normalize();
		if(!absolute.startsWith(evidenceRoot)){throw new IllegalArgumentException("Fixture output escaped root: "+path);}
		return evidenceRoot.relativize(absolute).toString();
	}
	/** Records a regular nonsymlink file under the evidence root; duplicate registrations coalesce. */
	private static void retain(Path path){
		if(!Files.isRegularFile(path) || Files.isSymbolicLink(path)){throw new IllegalArgumentException("Expected regular fixture evidence: "+path);}
		evidence.add(relative(path));
	}
	/** Writes with the caller's open options, then records the resulting regular file. */
	private static Path writeBytes(Path path, byte[] bytes, java.nio.file.OpenOption... options)throws java.io.IOException{
		Files.write(path, bytes, options);
		retain(path);
		return path;
	}
	/** Copies one fixture and records the target after successful creation. */
	private static Path copy(Path source, Path target)throws java.io.IOException{
		Files.copy(source, target);
		retain(target);
		return target;
	}
	/** Creates a hard-link alias for protection checks and records it as regular evidence. */
	private static Path hardLink(Path link, Path target)throws java.io.IOException{
		Files.createLink(link, target);
		retain(link);
		return link;
	}
	/** Creates a symbolic alias and records its declared target separately from regular evidence files. */
	private static Path symbolicLink(Path link, Path target)throws java.io.IOException{
		Files.createSymbolicLink(link, target);
		evidenceLinks.put(relative(link), target.toString());
		return link;
	}
	/** Writes the success marker and deterministic file/link inventories, then reports success to stdout. */
	private static void finish(String marker)throws Exception{
		write(evidenceRoot.resolve("PASS"), marker+"\n");
		final StringBuilder links=new StringBuilder();
		for(java.util.Map.Entry<String, String> entry:evidenceLinks.entrySet()){
			links.append(entry.getKey()).append('\t').append(entry.getValue()).append('\n');
		}
		write(evidenceRoot.resolve("evidence.links.tsv"), links.toString());
		evidence.add("evidence.paths");
		Files.write(evidenceRoot.resolve("evidence.paths"), evidence, StandardCharsets.UTF_8, java.nio.file.StandardOpenOption.CREATE_NEW);
		System.out.println(marker);
	}
	/** Fresh absolute directory bound by main; helpers reject evidence outside it. */
	private static Path evidenceRoot;
	/** Sorted relative regular-file paths created or retained by this sequential fixture. */
	private static final java.util.TreeSet<String> evidence=new java.util.TreeSet<String>();
	/** Sorted symbolic-link paths mapped to their declared target text. */
	private static final java.util.TreeMap<String, String> evidenceLinks=new java.util.TreeMap<String, String>();
	/** Measures the exact fixture identity needed by legacy/exact-verification cases. */
	private static String pin(Path path){return DigestSuffix.file(path.toString());}
	/** Counts an attempted check and fails immediately on a violated expectation. */
	private static void check(boolean valid, String why){
		checks++;
		if(!valid){throw new AssertionError(why);}
	}
	/** Attempted checks in this process, including rejected checks; main does not reset the counter. */
	private static int checks;
}
