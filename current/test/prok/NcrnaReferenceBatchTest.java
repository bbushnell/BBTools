package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;

import parse.LineParser1;
import prok.BuildNcrnaReference;
import prok.DraftNcrnaScenarios;
import prok.NcrnaEvidenceCatalog;
import prok.NcrnaReferenceEvidence;
import prok.NcrnaReferenceGff;
import prok.NcrnaReferenceGff.Row;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.PrepareCmSearchBatches;
import structures.ByteBuilder;

/** Synthetic native-report import cases; no CM search or real fingerprint claim.
 * Literal fixture provenance tokens test binding equality only. Original and
 * prepared FASTA payloads deliberately do not exist: capture import must not
 * reopen them. Separate batch-writer fixtures test actual sequence preservation.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaReferenceBatchTest{

	public static void main(final String[] args){
		String out=null;
		for(final String arg:args){if(arg.startsWith("out=")){out=arg.substring(4);}}
		check(out!=null && new File(out).mkdirs(), "Use out=fresh_fixture_directory");
		run(new File(out).getAbsolutePath());
		System.err.println("NCRNA_REFERENCE_BATCH_TEST_PASS");
	}

	private static void run(final String root){
		final Case good=new Case(root, "good", "good");
		good.importEvidence();
		metric(good.prefix+".summary.tsv", 0, 3);
		metric(good.prefix+".summary.tsv", 4, 4);
		metric(good.prefix+".summary.tsv", 10, 2);
		metric(good.prefix+".summary.tsv", 11, 2);
		metric(good.prefix+".summary.tsv", 13, 6);
		metric(good.prefix+".summary.tsv", 14, 3);
		checkHits(good);
		checkHistory(good);
		checkRetention(good);
		good.saved("profile", good.dir+"/profile");
		good.saved("geometry", good.dir+"/geometry");
		metric(good.dir+"/geometry.geometry_summary.tsv", 6, 2);
		metric(good.dir+"/geometry.geometry_summary.tsv", 11, 1);
		checkPairs(good);
		good.saved("scores", good.dir+"/scores");
		checkScoreSummary(good.dir+"/scores.score_summary.tsv", "RF00001", 2, 1, 1, 1);
		final Case scores=new Case(root, "scores", "scores");
		scores.importEvidence();
		scores.saved("scores", scores.dir+"/scores");
		checkCurve(scores.dir+"/scores.floor_curve.tsv", "RF00001", 85, false,
			new long[]{1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0});
		checkCurve(scores.dir+"/scores.floor_curve.tsv", "RF00001", 80, true,
			new long[]{1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1});
		checkCurve(scores.dir+"/scores.floor_curve.tsv", "RF00050", 80, false,
			new long[]{1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0});
		checkCurve(scores.dir+"/scores.floor_curve.tsv", "RF00050", 80, true,
			new long[]{1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0});
		check(rows(scores.dir+"/scores.score_histogram.tsv")==4, "Four distinct score/span/native-inclusion histogram cells");
		checkHistogram(scores.dir+"/scores.score_histogram.tsv", "RF00001", "[85,86)", 1);
		final Case bounds=new Case(root, "bounds", "bounds");
		bounds.importEvidence();
		bounds.saved("scores", bounds.dir+"/scores");
		// 0.85*47=39.95: both 40.0-0.05 and 39.9+0.05 touch the floor.
		checkCurve(bounds.dir+"/scores.floor_curve.tsv", "RF00001", 85, false,
			new long[]{2, 1, 0, 2, 0, 1, 1, 0, 1, 1, 1, 1, 0});
		checkCurve(bounds.dir+"/scores.floor_curve.tsv", "RF00050", 85, true,
			new long[]{2, 2, 2, 0, 0, 2, 2, 0, 0, 1, 1, 0, 1});
		checkHistogram(bounds.dir+"/scores.score_histogram.tsv", "RF00001", "[84,85)", 1);
		checkHistogram(bounds.dir+"/scores.score_histogram.tsv", "RF00001", "[85,86)", 1);
		final Case upper=new Case(root, "upper_bounds", "upper_bounds");
		upper.importEvidence();
		upper.saved("scores", upper.dir+"/scores");
		// Both competing hits have39.9: no printed pass, but upper-bound equality is possible.
		checkCurve(upper.dir+"/scores.floor_curve.tsv", "RF00050", 85, true,
			new long[]{2, 2, 2, 0, 0, 2, 2, 0, 0, 1, 0, 0, 1});
		good.saved("attributes", good.dir+"/attributes");
		BuildNcrnaReference.main(good.arguments("strata", good.dir+"/strata", "attributes="+good.dir+"/attributes"));
		check(rows(good.dir+"/strata.source_links.tsv.gz")==2, "Both cross-report historical links survive strata");
		// The old GA-only/p95 constructor must not silently reinterpret bound batch evidence.
		write(good.dir+"/unused_policy.tsv", "unused\n");
		write(good.dir+"/unused_priors.tsv", "unused\n");
		write(good.dir+"/unused_links.tsv", "unused\n");
		reject(new Runnable(){public void run(){
			DraftNcrnaScenarios.main(new String[]{"in="+good.prefix, "census="+good.dir+"/census.tsv",
				"memberships="+good.dir+"/memberships.tsv", "scope="+good.dir+"/scope.tsv", "assets="+good.dir+"/assets.tsv",
				"queries="+good.dir+"/queries.tsv", "policy="+good.dir+"/unused_policy.tsv", "priors="+good.dir+"/unused_priors.tsv",
				"links="+good.dir+"/unused_links.tsv", "out="+good.dir+"/unsupported"});
		}}, "does not support the captured batch evidence schema");
		check(!new File(good.dir+"/unsupported").exists(), "Unsupported decision tool must fail before creating output");

		final Case empty=new Case(root, "empty", "empty");
		empty.importEvidence();
		metric(empty.prefix+".summary.tsv", 0, 3);
		metric(empty.prefix+".summary.tsv", 4, 0);
		metric(empty.prefix+".summary.tsv", 13, 6);
		metric(empty.prefix+".summary.tsv", 14, 6);
		empty.saved("geometry", empty.dir+"/geometry");
		metric(empty.dir+"/geometry.geometry_summary.tsv", 6, 0);
		empty.saved("scores", empty.dir+"/scores");
		check(rows(empty.dir+"/scores.score_histogram.tsv")==0, "Empty evidence has no invented score bins");
		check(rows(empty.dir+"/scores.floor_curve.tsv")==44, "Both native-inclusion strata of both models retain all eleven zero-event floors");
		checkScoreSummary(empty.dir+"/scores.score_summary.tsv", "RF00001", 0, 0, 0, 0);
		checkAttributeReuse(good);
		checkTiers(root, good, scores, bounds, empty);
		final String[] variants={"missing_model", "failed_model", "unknown_target", "duplicate_map", "foreign_map",
			"ordinal_gap", "out_of_range", "wrong_inclusion", "wrong_binding", "wrong_space"};
		final String[] messages={"Missing model completion", "failed model search", "Unknown synthetic native target", "Repeated synthetic target",
			"foreign batch", "complete increasing ordinals", "leave the captured contig", "inclusion threshold must equal original GA",
			"exact accepted preparation binding", "search space differs"};
		for(int i=0; i<variants.length; i++){
			final Case bad=new Case(root, variants[i], variants[i]);
			reject(new Runnable(){public void run(){bad.importEvidence();}}, messages[i]);
			check(!new File(bad.prefix+".summary.tsv").exists(), "Rejected batch must not publish a completed evidence summary");
		}
		checkNativeRelocation(root);
	}

	/** Moving native bytes must preserve every witness and still reject rewritten footers. */
	private static void checkNativeRelocation(final String root){
		final Case test=new Case(root, "relocation", "good");
		test.importEvidence();
		final String physical=test.dir+"/recovered_%_\u03b2";
		check(new File(test.nativeDir).renameTo(new File(physical)) && !new File(test.nativeDir).exists(),
			"Relocation fixture must make every original native read path unavailable");
		final String row="0\t"+test.nativeDir+'\t'+physical+'\n';
		final String map=test.dir+"/native_paths.tsv", recovered=test.dir+"/recovered_evidence";
		write(map, BuildNcrnaReference.NATIVE_PATH_HEADER+'\n'+row);
		BuildNcrnaReference.main(test.arguments("batchevidence", recovered, "searchspacemb=0.001", "nativepaths="+map));
		for(final String suffix:new String[]{".contigs.tsv.gz", ".hits.tsv.gz", ".sources.tsv.gz", ".historical.tsv.gz",
			".overlaps.tsv.gz", ".relations.tsv.gz", ".coverage.tsv", ".summary.tsv", ".batch_sources.tsv"}){
			sameTable(test.prefix+suffix, recovered+suffix);
		}
		metric(recovered+".summary.tsv", 4, 4);
		metric(recovered+".summary.tsv", 13, 6);
		metric(recovered+".summary.tsv", 14, 3);
		rejectRelocation(test, "missing", "", "Missing physical-native mapping");
		rejectRelocation(test, "wrong_original", "0\t"+test.nativeDir+"_wrong\t"+physical+'\n', "differs from the collector");
		rejectRelocation(test, "duplicate", row+row, "Duplicate native relocation");
		final String unused=test.dir+"/unused_physical";
		check(new File(unused).mkdir(), "Extra-map negative needs a distinct existing physical directory");
		rejectRelocation(test, "extra", row+"1\t"+test.nativeDir+"_extra\t"+unused+'\n', "with no extras");
		reject(new Runnable(){public void run(){
			BuildNcrnaReference.main(test.arguments("profile", test.dir+"/wrong_mode", "nativepaths="+map));
		}}, "applies only to batchevidence");
		final File original=new File(physical+"/RF00001/hits.tbl"), saved=new File(physical+"/RF00001/hits.saved.tbl");
		check(original.renameTo(saved), "Preserve the literal fixture footer before testing a forged read-location claim");
		try(Input in=new Input(saved.getPath()); Output out=new Output(original.getPath())){
			for(byte[] bytes=in.reader.nextLine(); bytes!=null; bytes=in.reader.nextLine()){
				out.writer.println(new String(bytes, StandardCharsets.UTF_8).replace(test.nativeDir, physical));
			}
		}
		reject(new Runnable(){public void run(){
			BuildNcrnaReference.main(test.arguments("batchevidence", test.dir+"/rewritten_footer", "searchspacemb=0.001", "nativepaths="+map));
		}}, "names a different hit table");
		check(!new File(test.dir+"/rewritten_footer.summary.tsv").exists(), "A changed execution path cannot publish completed evidence");
		check(original.renameTo(new File(physical+"/RF00001/hits.altered.tbl")) && saved.renameTo(original),
			"Restore the original native fixture while retaining the rejected alteration");
		System.err.println("NCRNA_NATIVE_RELOCATION_TEST_PASS");
	}

	/** Invalid access maps fail before any completed import summary is published. */
	private static void rejectRelocation(final Case test, final String name, final String rows, final String reason){
		final String path=test.dir+'/'+name+".paths.tsv", output=test.dir+'/'+name;
		write(path, BuildNcrnaReference.NATIVE_PATH_HEADER+'\n'+rows);
		reject(new Runnable(){public void run(){
			BuildNcrnaReference.main(test.arguments("batchevidence", output, "searchspacemb=0.001", "nativepaths="+path));
		}}, reason);
		check(!new File(output+".summary.tsv").exists(), "Invalid relocation map must not publish a completed summary: "+name);
	}

	/** Reuse requires an explicit original evidence binding and exact current keys. */
	private static void checkAttributeReuse(final Case current){
		final String origin=legacy(current, "attribute_origin", false, true, false);
		final String capture=current.dir+"/legacy_capture";
		copyCapture(current, capture, 0);
		reject(new Runnable(){public void run(){
			BuildNcrnaReference.main(current.arguments("strata", current.dir+"/implicit_reuse", "attributes="+capture));
		}}, "population/status differs");
		BuildNcrnaReference.main(current.arguments("strata", current.dir+"/explicit_reuse", "attributes="+capture, "attributesevidence="+origin));
		check(rows(current.dir+"/explicit_reuse.event_strata.tsv.gz")==4, "Recorded attributes join all current events without a new GFF capture");
		copyCapture(current, current.dir+"/replaced_capture", 1);
		reject(new Runnable(){public void run(){
			BuildNcrnaReference.main(current.arguments("strata", current.dir+"/bad_replaced", "attributes="+current.dir+"/replaced_capture", "attributesevidence="+origin));
		}}, "lacks its historical source key");
		copyCapture(current, current.dir+"/duplicate_capture", 2);
		reject(new Runnable(){public void run(){
			BuildNcrnaReference.main(current.arguments("strata", current.dir+"/bad_capture_duplicate", "attributes="+current.dir+"/duplicate_capture", "attributesevidence="+origin));
		}}, "Repeated captured");
		check(!new File(current.dir+"/bad_replaced.strata_summary.tsv").exists()
			&& !new File(current.dir+"/bad_capture_duplicate.strata_summary.tsv").exists(), "Failed key bindings cannot publish completed strata");
	}
	private static void copyCapture(final Case current, final String prefix, final int change){
		try(Input in=new Input(current.dir+"/attributes.attribute_summary.tsv"); Output out=new Output(prefix+".attribute_summary.tsv")){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				out.writer.println(new String(row, StandardCharsets.UTF_8).replace(BuildNcrnaReference.BATCH_STATUS, "EVIDENCE_COMPLETE_NOT_ADJUDICATED"));
			}
		}
		final LineParser1 p=new LineParser1('\t');
		final ByteBuilder b=new ByteBuilder();
		try(Input in=new Input(current.dir+"/attributes.attributes.tsv.gz"); Output out=new Output(prefix+".attributes.tsv.gz")){
			out.writer.println(in.reader.nextLine());
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				if(change==1){
					p.set(row);
					b.clear();
					for(int field=0; field<4; field++){
						if(field>0){b.tab();}
						if(field==2){b.append(999);}else{p.appendTerm(b, field);}
					}
					out.writer.print(b.nl());
				}else{
					out.writer.println(row);
					if(change==2){out.writer.println(row);}
				}
			}
		}
	}

	/** Literal threshold/competition/prior oracles. Payloads are absent throughout. */
	private static void checkTiers(final String root, final Case good, final Case scores, final Case bounds, final Case empty){
		final String none=good.dir+"/no_priors.tsv";
		write(none, PRIOR_HEADER+'\n');
		for(final String model:new String[]{"RF00001", "RF00050"}){
			check(new File(good.table(model)).renameTo(new File(good.table(model)+"_not_needed")), "Tier construction must not reopen native alignment reports");
		}
		tiers(good, "tiers90", "90", none, new long[]{1, 3, 0, 0, 0, 0, 0}, 1, 3, 0);
		tiers(good, "tiers100", "100", none, new long[]{1, 1, 0, 2, 0, 0, 0}, 1, 1, 0);
		scores.saved("attributes", scores.dir+"/attributes");
		tiers(scores, "tiers90", "90", none, new long[]{2, 0, 2, 0, 0, 0, 0}, 2, 0, 0);
		tiers(scores, "tiers80", "80", none, new long[]{2, 1, 0, 1, 0, 0, 0}, 2, 1, 0);
		BuildNcrnaReference.main(new String[]{"config="+scores.dir+"/tiers90/rethreshold.config",
			"out="+scores.dir+"/from_config80", "marginalpercent=80", "pigz=f", "bgzip=f"});
		sameTable(scores.dir+"/tiers80/event_counts.tsv", scores.dir+"/from_config80/event_counts.tsv");
		tiers(scores, "tiers75", "75", none, new long[]{1, 3, 0, 0, 0, 0, 0}, 1, 3, 0);
		bounds.saved("attributes", bounds.dir+"/attributes");
		tiers(bounds, "tiers85", "85", none, new long[]{1, 1, 0, 2, 0, 0, 0}, 1, 1, 0);
		empty.saved("attributes", empty.dir+"/attributes");
		tiers(empty, "tiers90", "90", none, new long[]{0, 0, 0, 0, 0, 0, 0}, 0, 0, 0);
		metric(empty.dir+"/tiers90/source_family_counts.tsv", 11, 1);
		final Case tie=new Case(root, "native_tie", "tie_only");
		tie.importEvidence();
		tie.saved("attributes", tie.dir+"/attributes");
		tiers(tie, "tiers90", "90", none, new long[]{2, 1, 1, 0, 0, 0, 0}, 2, 1, 0);
		checkGaTie(tie.dir+"/tiers90/event_provenance.tsv.gz");
		final Case shallow=new Case(root, "shallow_scan", "shallow_scan");
		shallow.importEvidence();
		shallow.saved("attributes", shallow.dir+"/attributes");
		reject(new Runnable(){public void run(){
			BuildNcrnaReference.main(shallow.arguments("tiers", shallow.dir+"/unsupported75",
				"attributes="+shallow.dir+"/attributes", "priors="+none, "marginalpercent=90"));
		}}, "reporting floor");
		check(!new File(shallow.dir+"/unsupported75").exists(), "A shallower scan cannot promise75%re-threshold coverage");
		final Case partial=new Case(root, "partial", "partial");
		partial.importEvidence();
		partial.saved("attributes", partial.dir+"/attributes");
		tiers(partial, "tiers90", "90", none, new long[]{1, 0, 3, 0, 0, 0, 0}, 1, 0, 0);
		tiers(partial, "tiers80", "80", none, new long[]{1, 2, 0, 1, 0, 0, 0}, 1, 2, 0);
		final Case duplicates=new Case(root, "duplicates", "duplicates");
		duplicates.importEvidence();
		duplicates.saved("attributes", duplicates.dir+"/attributes");
		tiers(duplicates, "tiers90", "90", none, new long[]{1, 3, 0, 0, 0, 0, 0}, 1, 2, 0);
		final Case sameFamily=new Case(root, "same_family_models", "same_family_models");
		sameFamily.importEvidence();
		sameFamily.saved("attributes", sameFamily.dir+"/attributes");
		tiers(sameFamily, "tiers90", "90", none, new long[]{3, 1, 0, 0, 0, 0, 0}, 2, 1, 0);
		checkConcordance(root, none);
		final String prior=good.dir+"/accepted_priors.tsv";
		write(prior, PRIOR_HEADER+'\n'+prior(good, 1, "shared", 300, 60, 90)+prior(good, 3, "zero", 80, 10, 20));
		tiers(good, "tiers_priors", "90", prior, new long[]{1, 1, 0, 0, 1, 1, 0}, 3, 1, 2);
		final NcrnaReferenceGff codec=new NcrnaReferenceGff();
		final Row hitless=codec.read(good.dir+"/tiers_priors/input_3.gff.gz", "3").get(0);
		check(hitless.scores.missing && !hitless.mask && hitless.model.equals("RF00050"), "Hitless protected prior retains its original assignment and explicit unavailable scores");
		for(final Row row:codec.read(good.dir+"/tiers_priors/input_1.gff.gz", "1")){
			check(row.locus.start!=1 && row.locus.start!=50, "A prior cannot acquire widened bounds or a wrong-family masking exemption");
		}
		checkReportDetails(good, partial);
		System.err.println("NCRNA_REFERENCE_TIERS_TEST_PASS");
	}
	/** Independent loci cannot be swallowed by a stronger broad bridge, even when
	 * native report order changes. Weak witnesses never choose geometry. */
	private static void checkConcordance(final String root, final String priors){
		final String[] variants={"concordance_legit", "concordance_alt", "concordance_chain",
			"concordance_bridge", "concordance_bridge_reverse", "concordance_weak", "concordance_uncertain", "concordance_bridge_source_warning"};
		final NcrnaReferenceGff codec=new NcrnaReferenceGff();
		for(final String variant:variants){
			final Case test=new Case(root, variant, variant);
			test.importEvidence();
			test.saved("attributes", test.dir+"/attributes");
			final boolean weak=variant.equals("concordance_weak"), uncertain=variant.equals("concordance_uncertain");
			final boolean sourceWarning=variant.equals("concordance_bridge_source_warning");
			final boolean one=variant.equals("concordance_legit") || variant.equals("concordance_alt");
			final long[] counts={sourceWarning ? 2 : weak || uncertain ? 3 : 4, sourceWarning ? 2 : 0, weak ? 1 : 0, uncertain ? 1 : 0, 0, 0, 0};
			tiers(test, "tiers90", "90", priors, counts, one || sourceWarning ? 2 : 3, sourceWarning ? 1 : 0, 0);
			if(sourceWarning){
				reportReasons(test.dir+"/tiers90/report.reasons.tsv", new long[]{2, 0, 1, 0, 1});
				reportPhysical(test.dir+"/tiers90/report.physical.tsv", "first", new long[]{2, 1, 0, 1, 0, 0, 0, 0});
			}
			final java.util.ArrayList<Row> rows=codec.read(test.dir+"/tiers90/input_1.gff.gz", "1");
			check(rows.size()==(one ? 1 : 2), "Literal same-gene/tandem locus count: "+variant);
			if(variant.equals("concordance_legit")){
				check(rows.get(0).locus.start==5 && rows.get(0).locus.stop==104 && rows.get(0).scores.size()==3,
					"Legitimate cross-model support keeps the best original interval and every model witness");
			}else if(variant.equals("concordance_alt")){
				check(rows.get(0).locus.start==11 && rows.get(0).locus.stop==110,
					"Alternative same-model bounds retain the best native-included representative, without widening");
			}else if(variant.equals("concordance_chain")){
				check(rows.get(0).locus.start==1 && rows.get(1).locus.start==41,
					"A-B and B-C agreement cannot merge incompatible A-C bounds");
			}else{
				check(rows.get(0).locus.start==1 && rows.get(0).locus.stop==100
					&& rows.get(1).locus.start==201 && rows.get(1).locus.stop==300,
					"Tandem copy intervals survive a high-scoring bridge or weak witness regardless of report order");
				for(final Row row:rows){check(row.scores.size()==2, "Each copy retains its own support plus the bridge/weak witness");}
			}
		}
	}
	private static void sameTable(final String expected, final String actual){
		try(Input a=new Input(expected); Input b=new Input(actual)){
			for(byte[] row=a.reader.nextLine(); row!=null; row=a.reader.nextLine()){
				check(java.util.Arrays.equals(row, b.reader.nextLine()), "Decoded table differs from its expected records: "+actual);
			}
			check(b.reader.nextLine()==null, "Compared table must not acquire extra rows: "+actual);
		}
	}
	private static void checkGaTie(final String path){
		final LineParser1 p=new LineParser1('\t');
		int found=0;
		try(Input in=new Input(path)){
			in.reader.nextLine();
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				if(!p.termEquals("1", 1) || p.parseInt(3)!=201){continue;}
				check(p.parseDouble(16)==100 && p.termEquals("false", 17) && p.termEquals("true", 18)
					&& p.termEquals("MARGINAL", 22), "A native-excluded printed-GA tie must not become confident; disagreement stays visible");
				found++;
			}
		}
		check(found==1, "The isolated native-GA tie must remain independently identifiable");
	}
	private static String prior(final Case test, final int input, final String target, final int length, final int start, final int stop){
		final String fasta=input==1 ? test.dir+"/original_\u03b1.fa" : test.dir+"/original"+input+".fa";
		return "fixture\t/missing/fixture.review.tsv\t"+input+"\tprior"+input+"\tg"+input+'\t'+fasta+'\t'+input+'\t'+fasta+'\t'+target+'\t'+length
			+"\tRF00050\tACCEPTED_NOVEL_FMN\t"+start+'\t'+stop+"\t0\t0\t+\tPRIOR_REVIEWED_WHOLE_INTERVAL_NOT_NEW_APTAMER_BOUNDS\tfixture_review\tfixture_review\tfixture\n";
	}
	private static void tiers(final Case test, final String name, final String floor, final String priors,
		final long[] expected, final int scored, final int masked, final int protectedCount){
		final String out=test.dir+'/'+name;
		BuildNcrnaReference.main(test.arguments("tiers", out, "attributes="+test.dir+"/attributes", "priors="+priors, "marginalpercent="+floor));
		final LineParser1 p=new LineParser1('\t');
		long events=0;
		for(final long n:expected){events+=n;}
		int totals=0;
		try(Input in=new Input(out+"/event_counts.tsv")){
			in.reader.nextLine();
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				if(!p.termEquals("ALL_QUERIED", 0)){continue;}
				check(p.parseLong(1)==events && p.terms()==9, "All captured events survive the tier partition");
				for(int i=0; i<expected.length; i++){check(p.parseLong(i+2)==expected[i], "Literal tier oracle "+name+" state="+i);}
				totals++;
			}
		}
		check(totals==1 && rows(out+"/event_provenance.tsv.gz")==events, "Every raw event remains available for later re-thresholding");
		int observedScored=0, observedMasks=0, observedPriors=0;
		final NcrnaReferenceGff codec=new NcrnaReferenceGff();
		for(int input=1; input<=3; input++){
			for(final Row row:codec.read(out+"/input_"+input+".gff.gz", Integer.toString(input))){
				check(row.scores!=null, "Every primary annotation carries explicit score/GA/percent measurements or missing values");
				if(row.mask){observedMasks++;}else{observedScored++;}
			}
		}
		try(Input in=new Input(out+"/summary.tsv")){
			in.reader.nextLine();
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				observedPriors+=p.parseInt(3);
			}
		}
		check(observedScored==scored && observedMasks==masked && observedPriors==protectedCount, "Physical annotations and protected priors must match the literal oracle: "+name);
		check(rows(out+"/reference.tsv")==3, "Every input retains a reference manifest entry, including empty views");
		check(rows(out+"/source_provenance.tsv.gz")==1, "Unsupported and supported original source annotations both remain explicit provenance");
		BuildNcrnaReference.main(new String[]{"mode=tierreport", "in="+out, "out="+out+"/report",
			"samplemodel=RF00050", "examplespergroup=2", "pigz=f", "bgzip=f", "t=1"});
		metric(out+"/report.summary.tsv", 0, 3);
		metric(out+"/report.summary.tsv", 1, events);
		metric(out+"/report.summary.tsv", 2, scored);
		metric(out+"/report.summary.tsv", 3, masked);
		metric(out+"/report.summary.tsv", 4, protectedCount);
	}
	/** Literal reasons and prior-status oracles, independent of report counters. */
	private static void checkReportDetails(final Case good, final Case partial){
		reportReasons(good.dir+"/tiers90/report.reasons.tsv", new long[]{3, 0, 3, 0, 0});
		reportPhysical(good.dir+"/tiers90/report.physical.tsv", "first", new long[]{0, 2, 0, 2, 0, 0, 0, 0});
		reportReasons(partial.dir+"/tiers80/report.reasons.tsv", new long[]{2, 1, 0, 1, 0});
		reportPhysical(partial.dir+"/tiers80/report.physical.tsv", "first", new long[]{0, 2, 1, 0, 1, 0, 0, 0});
		metric(good.dir+"/tiers_priors/report.summary.tsv", 5, 1);
		metric(good.dir+"/tiers_priors/report.summary.tsv", 6, 1);
		reportPhysical(good.dir+"/tiers_priors/report.physical.tsv", "second", new long[]{3, 0, 0, 0, 0, 2, 1, 1});
		check(rows(good.dir+"/tiers90/report.examples.tsv")==2 && rows(good.dir+"/tiers90/report.example_overlaps.tsv")==2,
			"Bounded cross-family examples preserve two original overlaps and never join the same contig name across inputs");
		final LineParser1 p=new LineParser1('\t');
		int found=0;
		try(Input in=new Input(good.dir+"/tiers90/report.events.tsv")){
			in.reader.nextLine();
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				if(!p.termEquals("MODEL", 0) || !p.termEquals("RF00050", 2) || !p.termEquals("MARGINAL", 3)){continue;}
				check(p.parseLong(4)==1 && p.parseLong(7)==1 && p.parseLong(8)==1,
					"One marginal event has BOTH cross-family and source-warning causes; incidences are not additive");
				found++;
			}
		}
		check(found==1, "Expected one independently identified marginal model row");
		reject(new Runnable(){public void run(){
			BuildNcrnaReference.main(new String[]{"mode=tierreport", "in="+good.dir+"/tiers90",
				"out="+good.dir+"/unknown_sample", "samplemodel=RF99999"});
		}}, "not in the saved model population");
		final File original=new File(good.dir+"/tiers90/family_counts.tsv"), saved=new File(good.dir+"/tiers90/family_counts.original.tsv");
		check(original.renameTo(saved), "Preserve the literal family accounting before a negative fixture");
		final ByteBuilder b=new ByteBuilder();
		try(Input in=new Input(saved.getPath()); Output out=new Output(original.getPath())){
			out.writer.println(in.reader.nextLine());
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				if(!p.termEquals("first", 0)){
					out.writer.println(row);
					continue;
				}
				b.clear();
				for(int i=0; i<p.terms(); i++){
					if(i>0){b.tab();}
					if(i==9){b.append(p.parseLong(i)+1);}else{p.appendTerm(b, i);}
				}
				out.writer.print(b.nl());
			}
		}
		reject(new Runnable(){public void run(){
			BuildNcrnaReference.main(new String[]{"mode=tierreport", "in="+good.dir+"/tiers90", "out="+good.dir+"/wrong_physical"});
		}}, "Physical annotation/prior counts disagree");
		check(!new File(good.dir+"/wrong_physical.summary.tsv").exists(), "A mismatched physical population cannot publish a completed report");
		check(original.renameTo(new File(good.dir+"/tiers90/family_counts.altered.tsv")) && saved.renameTo(original), "Restore original fixture accounting without discarding failure evidence");
		System.err.println("NCRNA_REFERENCE_TIER_REPORT_TEST_PASS");
	}
	private static void reportReasons(final String path, final long[] expected){
		final LineParser1 p=new LineParser1('\t');
		int found=0;
		try(Input in=new Input(path)){
			in.reader.nextLine();
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				if(!p.termEquals("ALL", 0)){continue;}
				for(int i=0; i<expected.length; i++){check(p.parseLong(i+3)==expected[i], "Literal raw marginal reason count at "+i);}
				found++;
			}
		}
		check(found==1, "Exactly one whole-population reason record");
	}
	private static void reportPhysical(final String path, final String family, final long[] expected){
		final LineParser1 p=new LineParser1('\t');
		int found=0;
		try(Input in=new Input(path)){
			in.reader.nextLine();
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				if(!p.termEquals(family, 0)){continue;}
				for(int i=0; i<expected.length; i++){check(p.parseLong(i+1)==expected[i], "Literal physical mask/prior count at "+i+" for "+family);}
				found++;
			}
		}
		check(found==1, "Exactly one physical accounting row per declared family");
	}

	/** Literal multisets distinguish duplicate loss, changed scores and native flags. */
	private static void checkRetention(final Case current){
		final String old=legacy(current, "old", false, true, false);
		BuildNcrnaReference.main(current.arguments("retention", current.dir+"/retention", "compare="+old));
		retentionCounts(current.dir+"/retention.retention_summary.tsv",
			new long[]{3, 4, 2, 1, 2, 3, 3, 2, 1, 1, 1, 1, 1, 1, 0, 0});
		check(rows(current.dir+"/retention.retention_differences.tsv.gz")==4,
			"Duplicate loss, two new intervals and changed reverse-strand score need four coordinate groups");
		final String oldTie=legacy(current, "old_tie", true, true, false);
		BuildNcrnaReference.main(current.arguments("retention", current.dir+"/retention_tie", "compare="+oldTie));
		retentionCounts(current.dir+"/retention_tie.retention_summary.tsv",
			new long[]{4, 4, 3, 1, 1, 4, 3, 2, 2, 1, 1, 0, 2, 1, 1, 0});
		// Both old duplicates have printed GA_pass=true, but their actual native flags are '?'.
		final String oldExcluded=legacy(current, "old_excluded", false, true, true);
		BuildNcrnaReference.main(current.arguments("retention", current.dir+"/retention_excluded", "compare="+oldExcluded));
		retentionCounts(current.dir+"/retention_excluded.retention_summary.tsv",
			new long[]{3, 4, 2, 1, 2, 1, 3, 1, 0, 2, 1, 1, 1, 1, 0, 1});
		final String wrongCommand=legacy(current, "wrong_control", false, false, false);
		reject(new Runnable(){public void run(){
			BuildNcrnaReference.main(current.arguments("retention", current.dir+"/rejected_retention", "compare="+wrongCommand));
		}}, "requires exactly one --cut_ga");
		check(!new File(current.dir+"/rejected_retention.retention_summary.tsv").exists(),
			"A non-GA native control cannot publish a completed retention report");
	}

	/** Creates small legacy saved tables and their actual native witnesses, without searches. */
	private static String legacy(final Case current, final String name, final boolean tie, final boolean cutGa,
		final boolean excluded){
		final String dir=current.dir+'/'+name+"%_\u03b1";
		check(new File(dir).mkdir(), "Fresh legacy fixture directory");
		final String prefix=dir+"/evidence", first=dir+"/input1.tbl", second=dir+"/input2.tbl";
		final int count=tie ? 4 : 3, firstCount=count-1;
		final String footer="# Program: cmsearch\n# Version: 1.1.5 (fixture)\n# Option settings: cmsearch "
			+(cutGa ? "--cut_ga " : "")+"models.cm input.fa\n# [ok]\n";
		final String heading="# synthetic legacy native-format fixture\n# physical hit lines start at three\n";
		final String flag=excluded ? "?" : "!";
		final int duplicateBits=excluded ? 20 : 21;
		write(first, heading+hit("shared", "first", "RF00001", 1, 100, "+", 20, flag)
			+hit("shared", "first", "RF00001", 1, 100, "+", duplicateBits, flag)
			+(tie ? hit("shared", "first", "RF00001", 201, 300, "+", 20, "!") : "")+footer);
		write(second, heading+hit("shared", "second", "RF00050", 110, 11, "-", 31, "!")+footer);
		write(prefix+".summary.tsv", BuildNcrnaReference.SUMMARY_HEADER+"\n3\t3\t500\t3\t"+count
			+"\t3\t2\t1\t1\t0\t0\t0\t2\t6\t4\tEVIDENCE_COMPLETE_NOT_ADJUDICATED\n");
		write(prefix+".coverage.tsv", BuildNcrnaReference.COVERAGE_HEADER+"\nRF00001\t15.1\t100\t100\t0\t"
			+firstCount+'\t'+firstCount+"\t3\nRF00050\t15.1\t100\t100\t0\t1\t1\t3\n");
		try(Input in=new Input(current.prefix+".sources.tsv.gz"); Output out=new Output(prefix+".sources.tsv.gz")){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){out.writer.println(row);}
		}
		write(prefix+".overlaps.tsv.gz", BuildNcrnaReference.OVERLAP_HEADER+'\n');
		try(Output out=new Output(prefix+".hits.tsv.gz")){
			out.writer.println(BuildNcrnaReference.HIT_HEADER);
			legacyHit(out, 1, first, 3, 1, 100, '+', "RF00001", "first", 20);
			legacyHit(out, 1, first, 4, 1, 100, '+', "RF00001", "first", duplicateBits);
			if(tie){legacyHit(out, 1, first, 5, 201, 300, '+', "RF00001", "first", 20);}
			legacyHit(out, 2, second, 3, 11, 110, '-', "RF00050", "second", 31);
		}
		return prefix;
	}

	private static void legacyHit(final Output out, final int input, final String path, final int line,
		final int start, final int stop, final char strand, final String model, final String family, final int bits){
		final ByteBuilder b=new ByteBuilder();
		b.append(input).tab();
		BuildNcrnaReference.field(b, path).tab().append(line).append("\tshared\t");
		b.append(start).tab().append(stop).tab().append(strand).tab().append(model).append("\t15.1\t");
		b.append(family).append("\tTARGET_CLASS\t1\t100\t100\t0\t").append(bits);
		out.writer.print(b.append("\t1e-5\ttrue\tVALID_INTERVAL\tfalse\tfalse\n"));
	}

	private static void retentionCounts(final String path, final long[] expected){
		int found=0, models=0;
		final long[] sum=new long[expected.length];
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(path)){
			in.reader.nextLine();
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				check(p.terms()==18, "Retention summary schema");
				if(p.termEquals("ALL_QUERIED", 0)){
					for(int i=0; i<expected.length; i++){
						check(p.parseLong(i+1)==expected[i] && sum[i]==expected[i], "Literal retention oracle at metric "+i);
					}
					found++;
				}else{
					for(int i=0; i<expected.length; i++){sum[i]+=p.parseLong(i+1);}
					models++;
				}
			}
		}
		check(found==1 && models==2, "One whole-population and two per-model retention summaries");
	}

	/** Literal curve oracle: a reported weak bridge is not a qualifying competitor. */
	private static void checkCurve(final String path, final String model, final int floor,
		final boolean included, final long[] expected){
		int found=0;
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(path)){
			in.reader.nextLine();
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				if(!p.termEquals(model, 0) || p.parseInt(2)!=floor || !p.termEquals(Boolean.toString(included), 3)){continue;}
				check(p.terms()==17, "Floor curve schema");
				for(int i=0; i<expected.length; i++){check(p.parseLong(4+i)==expected[i], "Literal floor oracle "+model+" floor="+floor+" metric="+i);}
				found++;
			}
		}
		check(found==1, "Each model/floor/native-inclusion row occurs exactly once");
	}

	private static void checkScoreSummary(final String path, final String model,
		final long total, final long included, final long excluded, final long disagreements){
		int found=0;
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(path)){
			in.reader.nextLine();
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				if(!p.termEquals(model, 0)){continue;}
				check(p.terms()==7 && p.parseLong(2)==total && p.parseLong(3)==included
					&& p.parseLong(4)==excluded && p.parseLong(5)==disagreements,
					"Native inclusion remains distinct from printed score/GA agreement");
				found++;
			}
		}
		check(found==1, "Explicit summary row even for a zero-event model");
	}

	private static void checkHistogram(final String path, final String model, final String bin, final long expected){
		long count=0;
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(path)){
			in.reader.nextLine();
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				if(p.termEquals(model, 0) && p.termEquals(bin, 6)){count+=p.parseLong(7);}
			}
		}
		check(count==expected, "Exact decimal histogram boundary: "+model+' '+bin);
	}

	/** Two model tables use line3 in one input; both keys and the tied '?' must survive. */
	private static void checkHits(final Case test){
		final LineParser1 p=new LineParser1('\t');
		int first=0, bridge=0, rounded=0, reverse=0;
		try(Input in=new Input(test.prefix+".hits.tsv.gz")){
			check(new String(in.reader.nextLine(), StandardCharsets.UTF_8).equals(BuildNcrnaReference.BATCH_HIT_HEADER), "Batch hit schema is explicit");
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				check(p.terms()==24 && p.termEquals("shared", 3) && p.termEquals(BuildNcrnaReference.BOUND_GA_BASIS, 23), "Restored target/binding basis");
				if(p.termEquals("1", 0) && p.termEquals("RF00001", 7) && p.parseLong(2)==3){
					check(p.termEquals(test.table("RF00001"), 1) && p.termEquals("true", 17) && p.termEquals("true", 22), "First model's line3");
					first++;
				}else if(p.termEquals("1", 0) && p.termEquals("RF00050", 7)){
					check(p.parseLong(2)==3 && p.termEquals(test.table("RF00050"), 1), "Second model's independent line3");
					bridge++;
				}else if(p.termEquals("1", 0)){
					check(p.parseDouble(15)==20 && p.termEquals("false", 17) && p.termEquals("false", 22), "Printed GA tie must retain excluded native flag");
					rounded++;
				}else{
					check(p.termEquals("2", 0) && p.termEquals("-", 6) && p.parseInt(4)==11 && p.parseInt(5)==110
						&& p.termEquals("i2c1", 21), "Same token in second genome remains its reverse-strand event");
					reverse++;
				}
			}
		}
		check(first==1 && bridge==1 && rounded==1 && reverse==1, "All four distinct native events retained once");
	}

	private static void checkHistory(final Case test){
		try(Input in=new Input(test.prefix+".historical.tsv.gz")){
			check(new String(in.reader.nextLine(), StandardCharsets.UTF_8).equals(BuildNcrnaReference.BATCH_HISTORICAL_HEADER), "Batch history schema is explicit");
			final LineParser1 p=new LineParser1('\t');
			p.set(in.reader.nextLine());
			check(p.terms()==20 && p.parseLong(18)==3 && p.termEquals(test.table("RF00050"), 19), "Best witness needs its actual report as well as line3");
			check(in.reader.nextLine()==null, "One historical source feature");
		}
	}

	private static void checkPairs(final Case test){
		final LineParser1 p=new LineParser1('\t');
		int pairs=0;
		try(Input in=new Input(test.dir+"/geometry.model_pairs.tsv.gz")){
			in.reader.nextLine();
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				check(p.terms()==17 && !p.termEquals(p.parseString(1), 3), "Each cross-report pair retains both source paths");
				pairs++;
			}
		}
		check(pairs==2, "Bridge has two actual overlap edges");
		try(Input in=new Input(test.dir+"/geometry.bridges.tsv.gz")){
			in.reader.nextLine();
			p.set(in.reader.nextLine());
			check(p.terms()==17 && p.termEquals(test.table("RF00050"), 2) && p.parseLong(3)==3
				&& p.termEquals(test.table("RF00001"), 4) && p.parseLong(5)==3
				&& p.termEquals(test.table("RF00001"), 6) && p.parseLong(7)==4, "Bridge and both disjoint witnesses retain real table/line keys");
			check(in.reader.nextLine()==null, "One deterministic bridge witness");
		}
	}

	private static final class Case{
		Case(final String root, final String name, final String variant_){
			variant=variant_;
			dir=root+'/'+name;
			prepared=dir+"/prepared";
			nativeDir=dir+"/native";
			prefix=dir+"/evidence";
			check(new File(prepared).mkdirs() && new File(nativeDir).mkdir(), "Fresh fixture directories");
			final String a=dir+"/original_\u03b1.fa", b=dir+"/original2.fa", c=dir+"/original3.fa";
			final String members="1\t"+a+"\t10\t300\t1\n2\t"+b+"\t10\t120\t1\n3\t"+c+"\t10\t80\t1\n";
			write(dir+"/census.tsv", NcrnaEvidenceCatalog.CENSUS_HEADER+'\n'+members);
			write(prepared+"/members.tsv", NcrnaEvidenceCatalog.CENSUS_HEADER+'\n'+members);
			write(dir+"/source.gff", "##gff-version 3\nshared\tfixture\tncRNA\t60\t90\t.\t+\t.\tID=source;Dbxref=RFAM:RF00001"
				+(variant.equals("concordance_bridge_source_warning") ? ";pseudo=true" : "")+'\n');
			write(dir+"/memberships.tsv", NcrnaEvidenceCatalog.MEMBERSHIP_HEADER+'\n'
				+"fixture\t1\t1\t"+a+"\tg1\ta1\t1\tPROK\t"+dir+"/source.gff\tNA\n"
				+"fixture\t2\t2\t"+b+"\tg2\ta2\t2\tPROK\tNA\tNA\n"
				+"fixture\t3\t3\t"+c+"\tg3\ta3\t3\tPROK\tNA\tNA\n");
			final boolean concordance=variant.startsWith("concordance_");
			write(dir+"/scope.tsv", NcrnaReferenceEvidence.SCOPE_HEADER+"\nRF00001\tfirst\tPROK\tfixture\tfirst\tTARGET_CLASS\nRF00050\tsecond\tPROK\tfixture\t"
				+(variant.equals("same_family_models") || concordance ? "first" : "second")+"\tTARGET_CLASS\n");
			final boolean boundary=variant.equals("bounds") || variant.equals("upper_bounds");
			final boolean partial=variant.equals("partial"), duplicate=variant.equals("duplicates");
			final boolean sameFamily=variant.equals("same_family_models");
			final boolean low=variant.equals("scores") || partial;
			final boolean weakBridge=low || variant.equals("tie_only");
			final int firstGa=boundary ? 47 : 20;
			final double firstLower=boundary ? 35.25 : variant.equals("shallow_scan") ? 16 : 15;
			final double secondLower=variant.equals("shallow_scan") ? 24 : 22.5;
			final double firstBits=partial ? 17 : variant.equals("upper_bounds") ? 39.9 : boundary ? 40 : 20;
			write(dir+"/assets.tsv", NcrnaReferenceEvidence.ASSET_HEADER+"\nRF00001\t15.1\t100\t"+firstGa+'\t'+SHA+"\t"+model("RF00001")
				+"\nRF00050\t15.1\t100\t30\t"+SHA+"\t"+model("RF00050")+'\n');
			write(dir+"/queries.tsv", NcrnaReferenceEvidence.QUERY_HEADER+"\nRF00001\tfirst\t100\nRF00050\tsecond\t100\n");
			write(prepared+"/summary.tsv", "batch_id\ttask_id\tgenomes\tnucleotide_bases\tfasta_records\tstatus\n0\t0\t3\t500\t3\tPREPARED_NOT_SEARCHED\n");
			final String map=(variant.equals("foreign_map") ? "1" : "0")+"\ti1c1\t1\t1\tshared\t300\n"
				+(variant.equals("ordinal_gap") ? "0\ti2c2\t2\t2" : "0\ti2c1\t2\t1")+"\tshared\t120\n0\ti3c1\t3\t1\tzero\t80\n";
			write(prepared+"/targets.tsv", PrepareCmSearchBatches.TARGET_HEADER+'\n'+map
				+(variant.equals("duplicate_map") ? "0\ti1c1\t1\t1\tshared\t300\n" : ""));
			final String binding="0\t"+prepared+"/batch.fa\t"+SHA+'\t'+prepared+"/targets.tsv\t"+SHA+'\t'
				+prepared+"/members.tsv\t"+SHA+"\t3\t500\t3\tPREPARATION_BOUND\n";
			write(prepared+"/binding.tsv", BINDING+'\n'+binding);
			write(nativeDir+"/input.tsv", BINDING+'\n'+(variant.equals("wrong_binding") ? binding.replace(SHA, OTHER_SHA) : binding));
			write(nativeDir+"/exit.txt", "0\n");
			final boolean empty=variant.equals("empty");
			write(nativeDir+"/completion.tsv", COMPLETION+"\n0\tRF00001\t3\t500\t1000\t"+(empty ? 0 : 2)+'\t'+firstGa+'\t'+firstLower+"\tCOMPLETE\n"
				+(variant.equals("missing_model") ? "" : "0\tRF00050\t3\t500\t1000\t"+(empty ? 0 : 2)+"\t30\t"+secondLower+"\tCOMPLETE\n"));
			final String firstTarget=variant.equals("unknown_target") ? "i99c1" : "i1c1";
			final String firstFlag=boundary || partial ? "?" : "!";
			final double secondBits=boundary ? 39.9 : low ? 17 : 20;
			String firstHit=hit(firstTarget, "first", "RF00001", 1, 100, "+", firstBits, firstFlag);
			String secondHit=hit("i1c1", "first", "RF00001", duplicate ? 1 : 201,
				duplicate ? 100 : 300, "+", secondBits, duplicate ? "!" : "?");
			String bridge=hit("i1c1", "second", "RF00050", sameFamily ? 1 : 50,
				sameFamily ? 100 : 250, "+", weakBridge ? 24 : 40, weakBridge ? "?" : "!");
			if(concordance){
				final boolean chain=variant.equals("concordance_chain"), alt=variant.equals("concordance_alt");
				final boolean legit=variant.equals("concordance_legit");
				final boolean weak=variant.equals("concordance_weak"), uncertain=variant.equals("concordance_uncertain");
				final int secondStart=chain ? 41 : alt || legit ? 11 : 201;
				final int bridgeStart=chain ? 21 : legit ? 5 : 1;
				final int bridgeStop=chain ? 120 : legit ? 104 : alt ? 100 : 300;
				final double aBits=chain ? 30 : 20, bBits=alt ? 30 : 22;
				final double bridgeBits=weak ? 24 : uncertain ? 27 : chain || alt || legit ? 40 : 90;
				firstHit=hit("i1c1", "first", "RF00001", 1, 100, "+", aBits, "!");
				secondHit=hit("i1c1", "first", "RF00001", secondStart, secondStart+99, "+", bBits, "!");
				bridge=hit("i1c1", "second", "RF00050", bridgeStart, bridgeStop, "+", bridgeBits, weak || uncertain ? "?" : "!");
			}
			final String reverse=hit("i2c1", "second", "RF00050", variant.equals("out_of_range") ? 130 : 110, 11, "-", 30, "!");
			final String forward=variant.equals("concordance_bridge_reverse") ? secondHit+firstHit : firstHit+secondHit;
			modelOutput("RF00001", "first", firstGa, firstLower, empty ? "" : forward, empty ? 0 : 2);
			modelOutput("RF00050", "second", 30, secondLower, empty ? "" : bridge+reverse, empty ? 0 : 2);
			write(dir+"/collector.tsv", BuildNcrnaReference.BATCH_COLLECTOR_HEADER+"\n0\t"+prepared+'\t'+nativeDir+"\t2\t"+(empty ? 0 : 4)+'\n');
		}
		void modelOutput(final String accession, final String name, final double ga, final double lower, final String hits, final int count){
			final String point=nativeDir+'/'+accession;
			check(new File(point).mkdir(), "Fresh native model point");
			write(point+"/exit.txt", variant.equals("failed_model") && accession.equals("RF00050") ? "1\n" : "0\n");
			final double inclusion=variant.equals("wrong_inclusion") ? ga-1 : ga;
			final String command="cmsearch --cpu 1 --default --nohmmonly --noali --notextw -Z "+(variant.equals("wrong_space") ? "0.002" : "0.001")
				+" -T "+lower+" --incT "+inclusion+" --tformat fasta --tblout "+table(accession)+" -o "+point+"/search.out "+model(accession)+' '+prepared+"/batch.fa";
			final String printedHits=variant.equals("partial") && accession.equals("RF00001") ? hits.replace("cm 1 100 1 100 +", "cm 1 80 1 100 +") : hits;
			write(table(accession), "# synthetic native-format fixture\n# physical hit lines start at three\n"+printedHits
				+"# Program: cmsearch\n# Version: 1.1.5 (fixture)\n# Option settings: "+command+"\n# [ok]\n");
			write(point+"/search.out", "Query: "+name+" [CLEN=100]\nAccession: "+accession
				+"\nTarget sequences: 3 (1000 residues searched)\nTotal CM hits reported: "+count+"\n[ok]\n");
		}
		String model(final String accession){return dir+'/'+accession+".cm";}
		String table(final String accession){return nativeDir+'/'+accession+"/hits.tbl";}
		void importEvidence(){BuildNcrnaReference.main(arguments("batchevidence", prefix, "searchspacemb=0.001"));}
		void saved(final String mode, final String output){BuildNcrnaReference.main(arguments(mode, output));}
		String[] arguments(final String mode, final String output, final String... extra){
			final String[] args=new String[11+extra.length];
			args[0]="mode="+mode;
			args[1]="in="+(mode.equals("batchevidence") ? dir+"/collector.tsv" : prefix);
			args[2]="out="+output;
			args[3]="census="+dir+"/census.tsv";
			args[4]="memberships="+dir+"/memberships.tsv";
			args[5]="scope="+dir+"/scope.tsv";
			args[6]="assets="+dir+"/assets.tsv";
			args[7]="queries="+dir+"/queries.tsv";
			args[8]="t=1";
			args[9]="pigz=f";
			args[10]="bgzip=f";
			System.arraycopy(extra, 0, args, 11, extra.length);
			return args;
		}
		final String dir, prepared, nativeDir, prefix, variant;
	}

	private static String hit(final String target, final String name, final String model, final int from, final int to,
		final String strand, final double bits, final String included){
		return target+" - "+name+' '+model+" cm 1 100 "+from+' '+to+' '+strand+" no 1 0.5 0.0 "+bits+" 1e-5 "+included+" fixture\n";
	}
	private static void metric(final String path, final int column, final long expected){
		try(Input in=new Input(path)){
			in.reader.nextLine();
			final LineParser1 p=new LineParser1('\t');
			p.set(in.reader.nextLine());
			check(p.parseLong(column)==expected && in.reader.nextLine()==null, "Literal summary oracle: "+path+" column "+column);
		}
	}
	private static int rows(final String path){
		int count=0;
		try(Input in=new Input(path)){
			in.reader.nextLine();
			while(in.reader.nextLine()!=null){count++;}
		}
		return count;
	}
	private static void write(final String path, final String text){try(Output out=new Output(path)){out.writer.print(text.getBytes(StandardCharsets.UTF_8));}}
	private static void reject(final Runnable action, final String fragment){
		try{action.run();}
		catch(RuntimeException failure){
			for(Throwable cause=failure; cause!=null; cause=cause.getCause()){
				if(cause.getMessage()!=null && cause.getMessage().contains(fragment)){return;}
			}
			throw new AssertionError("Wrong rejection for "+fragment, failure);
		}
		throw new AssertionError("Expected rejection: "+fragment);
	}
	private static void check(final boolean condition, final String message){if(!condition){throw new AssertionError(message);}}
	private static final String SHA="00000000000000000001", OTHER_SHA="00000000000000000002";
	private static final String BINDING="batch_id\tfasta\tfasta_sha80\tmapping\tmapping_sha80\tmembers\tmembers_sha80\tgenomes\tnucleotide_bases\tfasta_records\tstatus";
	private static final String COMPLETION="batch_id\tmodel\tFASTA_records\tinput_bases\tsearched_residues\treported_hits\toriginal_GA_bits\treport_bits\tstatus";
	private static final String PRIOR_HEADER="source_id\tdecision_file\tdecision_line\toriginal_locus\toriginal_genome\toriginal_fasta\tinput_id\tresolved_fasta\tliteral_target\tcontig_length\tsource_family\tassignment\tsource_start1\tsource_stop1\tcore_start1\tcore_stop1\tstrand\tinterval_semantics\tsource_authority\trow_authority\toriginal_row";
}
