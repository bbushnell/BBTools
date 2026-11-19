package test.prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Objects;

import ddl.DDLFamilyProfiles;
import prot.DigestSuffix;
import prot.FamilyAlignmentMetric;
import prot.FamilyHeuristicScore;
import prot.FamilyHitFloorVector;
import prot.FamilyHitHistogram;
import prot.FamilyReference;
import prot.FamilyReferenceCodec;
import prot.FamilyTailCascade;
import prot.FamilyTailCascade.Assignment;
import structures.ByteBuilder;
import prot.ProteinSequence;

/** Literal admission, binding, observation and unchanged-control checks.
 * The three reference profiles have exactly two, five and zero query matches.
 * Mutable assertion counts belong to one isolated test JVM, not concurrent calls.
 * @author Brian, Keqing, Fischl
 */
public final class FamilyHitFloorIntegrationTest{

	/*--------------------------------------------------------------*/
	/*----------------      Admission checks        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates synthetic files under a fresh root and checks floor-policy behavior.
	 * Keeps two independently loaded indexes resident to test binding identity;
	 * the public launcher supplies a 4-GiB default heap for this mode.
	 * @param args one output directory whose parent already exists
	 * @throws Exception if fixture files cannot be created or loaded
	 */
	public static void main(final String[] args)throws Exception{
		check(args.length==1, "Supply a fresh Dori fixture directory");
		final Path root=Paths.get(args[0]).toAbsolutePath();
		Files.createDirectory(root);
		final byte[] query=sequence();
		final FamilyTailCascade.LoadedClassifier base=fixture(root.resolve("parent"), query);
		final FamilyHitFloorVector probe=base.newFamilyHitFloorProbe();
		final FamilyTailCascade.LoadedClassifier allOne=base.withFamilyHitFloors(probe);
		check(allOne.usesFamilyHitFloors() && !allOne.usesFamilyMeanAdmission() && !base.usesFamilyHitFloors(),
			"Family-specific floors have a distinct policy identity and leave the fixed parent unchanged");
		final Assignment parent=base.newWorker().assign("parent", query);
		check(!parent.assigned() && parent.admittedFamilies==0 && parent.rawHitFloor==6,
			"The five-hit best family cannot pass the fixed-six control");
		final Assignment open=allOne.newWorker().assign("probe", query);
		check(open.assigned() && open.familyId.equals("f1") && open.ddlMatches==5 && open.admittedFamilies==2 &&
			open.familyHitFloorAdmission && !open.familyMeanAdmission && open.rawHitFloor==1 && open.rawHitFloorMax==1,
			"Both positive-hit families enter the pool and the five-hit family wins the unchanged ranking/alignment");
		final long[] fittingWinners=new long[21];
		fittingWinners[2]=1;
		fittingWinners[12]=1;
		final FamilyHitFloorVector supported=probe.withInitialSupport(2, fittingWinners);
		final FamilyTailCascade.LoadedClassifier mixed=base.withFamilyHitFloors(supported.withFloors(new byte[]{3, 5, 6}));
		final Assignment admitted=mixed.newWorker().assign("mixed", query);
		check(admitted.familyId.equals("f1") && admitted.admittedFamilies==1 &&
			admitted.rawHitFloor==3 && admitted.rawHitFloorMax==6,
			"Heterogeneous floors apply before shortlist selection; the two-hit family fails its own floor");
		final Assignment none=base.withFamilyHitFloors(supported.withFloors(new byte[]{3, 6, 6}))
			.newWorker().assign("none", query);
		check(!none.assigned() && none.admittedFamilies==0 && none.terminal.equals("ddl_reject"),
			"A maximum of five exceeding the minimum floor three is insufficient when every family's own floor fails");
		same(parent, base.withFamilyHitFloors(supported.withFloors(new byte[]{6, 6, 6})).newWorker().assign("six", query));
		reject(()->mixed.withEvaluationMarginScale(.5), "Attached floors cannot silently acquire different margins");
		reject(()->base.withEvaluationMarginScale(.5).withFamilyHitFloors(supported),
			"The original vector cannot transfer to a changed parent margin");
		reject(()->mixed.withFamilyMeanAdmission(), "Mean-length and family-specific policies cannot be combined");
		reject(()->mixed.withLengthScaledShortAdmission(),
			"Query-length and family-specific policies cannot be combined");
		final FamilyTailCascade.LoadedClassifier other=fixture(root.resolve("other"), query);
		reject(()->other.withFamilyHitFloors(supported),
			"An independently loaded same-byte roster/index is not this resident binding");
		// Public admission and observation checks do not load historical caller-gap diagnostics.
		observe(base, mixed, query);
		System.out.println("FAMILY_HIT_FLOOR_INTEGRATION_TEST_PASS checks="+checks+" scientific_runs=0");
	}

	/** Checks raw and terminal observations, normalization and callback reentry.
	 * Classifiers are borrowed; this method owns its mutable workers and histograms.
	 * A rejected query must not corrupt counters or subsequent assignments.
	 */
	private static void observe(final FamilyTailCascade.LoadedClassifier base,
		final FamilyTailCascade.LoadedClassifier mixed, final byte[] query){
		final FamilyTailCascade.Worker worker=mixed.newWorker();
		final FamilyHitHistogram hist=new FamilyHitHistogram(3, true);
		final FamilyTailCascade.AssignmentObserver terminal=
			(length, state, family, hits, identity, blosum)->hist.finishQuery(family, hits);
		final FamilyTailCascade.RawHitObserver raw=hist::beginQuery;
		worker.assignObserve("accepted", query, terminal, raw);
		final byte[] unknown=new byte[64];
		Arrays.fill(unknown, (byte)'X');
		worker.assignObserve("zero", unknown, terminal, raw);
		final long[] expectedRaw=new long[21];
		expectedRaw[0]=expectedRaw[2]=expectedRaw[7]=expectedRaw[12]=1;
		expectedRaw[14]=2;
		final long[] expectedWinners=new long[21];
		expectedWinners[12]=1;
		check(hist.completedQueries()==2 && hist.acceptedQueries()==1 &&
			Arrays.equals(hist.copyRawBins(), expectedRaw) && Arrays.equals(hist.copyWinnerBins(), expectedWinners),
			"Raw observations include rejected and below-floor families; all-zero queries contribute only derived zeros");
		final FamilyHitHistogram rejected=new FamilyHitHistogram(3, true);
		base.newWorker().assignObserve("fixed_reject", query,
			(l, s, f, h, i, b)->rejected.finishQuery(f, h), rejected::beginQuery);
		check(rejected.completedQueries()==1 && rejected.acceptedQueries()==0 &&
			rejected.copyRawBins()[2]==1 && rejected.copyRawBins()[12]==1,
			"The fixed-six control observes uncensored hits even when admission rejects the entire pool");
		final byte[] terminated=Arrays.copyOf(query, query.length+1);
		terminated[query.length]='*';
		same(worker.assign("normalized", query), worker.assign("terminal_stop", terminated));
		reject(()->worker.assignObserve("short", Arrays.copyOf(query, 39), terminal, raw),
			"Unsupported queries must not become measured zero-hit records");
		check(hist.completedQueries()==2, "A pre-index failure must not increment the measured-query denominator");
		reject(()->worker.assignObserve("outer", query, (l, s, f, h, i, b)->{},
			(t, s, n)->worker.assign("nested", query)), "Raw callbacks cannot reenter mutable query scratch");
		reject(()->worker.assignObserve("outer_terminal", query, (l, s, f, h, i, b)->worker.assign("nested", query)),
			"Terminal callbacks cannot replace the result before observation finishes");
		same(mixed.newWorker().assign("independent", query), worker.assign("after_callback_failure", query));
	}

	/*--------------------------------------------------------------*/
	/*----------------       Shared fixture         ----------------*/
	/*--------------------------------------------------------------*/

	/** Writes a three-family synthetic package and loads it through the public API.
	 * The default profiles retain two, five and zero occupied query registers.
	 * @param root fresh directory beneath an existing parent
	 * @param raw amino-acid query used for consensus and positional profiles
	 * @return a new classifier with the fixture's member-length ranking policy
	 * @throws Exception if fixture files or classifier validation fail
	 */
	public static FamilyTailCascade.LoadedClassifier fixture(final Path root, final byte[] raw)throws Exception{
		return fixture(root, raw, new int[]{2, 5, 0});
	}
	/** Builds the same package with explicit occupied-register counts per family.
	 * The fitting CLI fixture uses this overload to exercise a failed fixed-six
	 * endpoint. The supplied counts are read in family order.
	 * @param root fresh directory beneath an existing parent
	 * @param raw query with enough occupied registers for every requested count
	 * @param hits three nonnegative counts, in family order
	 * @return a newly loaded classifier bound to the written synthetic package
	 * @throws Exception if a count cannot be realized or package creation fails
	 */
	public static FamilyTailCascade.LoadedClassifier fixture(final Path root, final byte[] raw,
		final int[] hits)throws Exception{
		check(hits!=null && hits.length==3 && hits[0]>=0 && hits[1]>=0 && hits[2]>=0,
			"The three-family fixture needs explicit nonnegative self-hit counts");
		Files.createDirectory(root);
		final char[] sketch=new char[1024];
		DDLFamilyProfiles.sketchPositional(raw, 5, 0, 4, 2, sketch);
		final FamilyAlignmentMetric metric=new FamilyAlignmentMetric(
			FamilyAlignmentMetric.Mode.PIN_RIGHT_TRIM_LEFT_ALIGNED_QUERY, 28);
		final FamilyReferenceCodec.Header header=new FamilyReferenceCodec.Header(FamilyReferenceCodec.ALIGNMENT_SCHEMA,
			5, 1024, 0, 3, 3, 0, FamilyReferenceCodec.FORMULA_VERSION,
			FamilyReferenceCodec.SPECTRUM_ENCODING, 3.8906, metric);
		final ByteBuilder refs=FamilyReferenceCodec.appendHeader(header, new ByteBuilder());
		final ByteBuilder profiles=new ByteBuilder().append(DDLFamilyProfiles.HEADER16).nl();
		for(int f=0; f<3; f++){
			final char[] row=sketch.clone();
			int available=0;
			for(int bin=0; bin<row.length; bin++){
				if(row[bin]!=0 && ++available>hits[f]){row[bin]=0;}
			}
			check(available>=hits[f], "The fixture query must supply at least five distinct positional registers");
			final String id="f"+f;
			final float[] spectrum=new float[400];
			spectrum[0]=1;
			FamilyReferenceCodec.appendFiveLines(new FamilyReference(id, f, 30, 1920, 1800,
				new ProteinSequence(id, raw), row, spectrum,
				0, 30, 0, 1.1f, 99, 1, raw, 0, new double[10],
				FamilyReference.AlignmentCalibration.ENABLED), refs, header);
			profiles.append("tail\t").append(id).append("\t1024\t30\t1920\t1800\tnonempty\t");
			for(int field=0; field<3; field++){
				if(field>0){profiles.tab();}
				for(int bin=0; bin<row.length; bin++){
					if(bin>0){profiles.append(',');}
					profiles.append(field==0 ? (int)row[bin] : row[bin]==0 ? 0 : 1);
				}
			}
			profiles.nl();
		}
		final Path ref=root.resolve("references.tsv"), data=root.resolve("profiles.tsv"),
			summary=root.resolve("summary.tsv"),
			list=root.resolve("list.tsv"), bounds=root.resolve("bounds.tsv"),
			medians=root.resolve("medians.tsv"), hydro=root.resolve("hydro.tsv");
		Files.write(ref, refs.toBytes());
		Files.write(data, profiles.toBytes());
		write(summary, "namespace\tk\tseed\texponent\tregisterbits\tsketch\tpartitions\toverlap\ttag_bits\tsegment_words\n"+
			"tail\t5\t0\t6\t16\t"+DDLFamilyProfiles.SKETCH_POSITIONAL+"\t4\t2\t"+
			DDLFamilyProfiles.positionalTagBits(4)+"\t"+DDLFamilyProfiles.positionalSegmentWords(4)+"\n");
		write(list, data+"\t"+pin(data)+"\t"+pin(summary)+"\n");
		final StringBuilder scalar=new StringBuilder("scalar\trows\tmin\tmax\tmin_family\tmax_family\n");
		for(String name:new String[]{"hydrophobicity", "net_charge", "aromaticity", "fraction_small", "pI_tail"}){
			scalar.append(name).append("\t79565\t-20\t20\tlo\thi\n");
		}
		write(bounds, scalar.toString());
		final StringBuilder ms=new StringBuilder("family\tmembers\tlength_p10\tlength_median\tlength_p90\tcount_under_200\n");
		final StringBuilder hs=new StringBuilder("family_rank\tfamily_id\tN\thydrophobicity_mean\n");
		for(int f=0; f<3; f++){
			ms.append('f').append(f).append("\t30\t64\t64\t64\t30\n");
			hs.append(f).append("\tf").append(f).append("\t30\t1.8\n");
		}
		write(medians, ms.toString());
		write(hydro, hs.toString());
		final HashMap<String, String> options=new HashMap<String, String>();
		for(String entry:new String[]{"top=50", "k=5", "seed=0", "width=1024", "ddlreject=6", "ddlaccept=1025",
			"positional=t", "partitions=4", "overlap=2", "indexearlyreject=f", "cosine=ranking", "rankingfeatures=combined25",
			"scalargates=t", "scalarz=3.8906", "identityacceptfloor=30", "blosumacceptfloor=0.145", "scalarscale=0.1",
			"minmemberaa=40", "maxmemberaa=99", "alignmentmode=pin_right_trim_left_aligned_query", "alignmentminspan=28"}){
			final int split=entry.indexOf('=');
			options.put(entry.substring(0, split), entry.substring(split+1));
		}
		options.put("familyreferences", ref.toString());
		options.put("familyreferencessha80", pin(ref));
		options.put("ddlprofiles", list.toString());
		options.put("ddlprofilessha80", pin(list));
		options.put("scalarbounds", bounds.toString());
		options.put("scalarboundssha80", pin(bounds));
		final FamilyTailCascade.LoadedClassifier base=FamilyTailCascade.loadCallable(options, 1);
		final StringBuilder recipe=new StringBuilder();
		for(String key:new java.util.TreeSet<String>(options.keySet())){
			recipe.append(key).append('=').append(options.get(key)).append('\n');
		}
		recipe.append("heuristic=hitfourth_h020_c150\nminqueryaa=40\nmaxqueryaa=99\nlengthmedians=").append(medians)
			.append("\nlengthmedianssha80=").append(pin(medians)).append("\nfamilyhydro=").append(hydro)
			.append("\nfamilyhydrosha80=").append(pin(hydro)).append('\n');
		write(root.resolve("classifier.config"), recipe.toString());
		return base.withHeuristicRanking(new FamilyHeuristicScore(
			base, medians, pin(medians), hydro, pin(hydro), 40, 99, "hitfourth_h020_c150"));
	}
	/** Returns a new deterministic 64-residue query over the canonical amino alphabet. */
	public static byte[] sequence(){
		final byte[] alphabet="ACDEFGHIKLMNPQRSTVWY".getBytes(StandardCharsets.US_ASCII), raw=new byte[64];
		for(int i=0; i<raw.length; i++){raw[i]=alphabet[(i*7+i/20)%alphabet.length];}
		return raw;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Test utilities         ----------------*/
	/*--------------------------------------------------------------*/

	/** Compares the assignment and measured stage fields used by the unchanged control. */
	private static void same(final Assignment a, final Assignment b){
		check(Objects.equals(a.familyId, b.familyId) && a.terminal.equals(b.terminal) && a.rosterRank==b.rosterRank &&
			a.ddlMatches==b.ddlMatches && a.rawHitFloor==b.rawHitFloor && a.rawHitFloorMax==b.rawHitFloorMax &&
			a.admittedFamilies==b.admittedFamilies && a.dpAttempts==b.dpAttempts &&
			a.identity==b.identity && a.blosum==b.blosum &&
			a.score==b.score && a.alignmentState==b.alignmentState,
			"The control must retain its ordinary assignment and measured stage values");
	}
	/** Writes the fixture text as ASCII; all paths are inside its fresh root. */
	private static void write(final Path path, final String text)throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII));
	}
	/** Returns the digest suffix required by the fixture package binding. */
	private static String pin(final Path path){return DigestSuffix.file(path.toString());}
	/** An operation expected to reject an incompatible fixture policy. */
	private interface Action{
		/** Executes the operation, allowing checked construction failures to propagate. */
		void run()throws Exception;
	}
	/** Counts only IllegalArgumentException as the intended policy rejection. */
	private static void reject(final Action action, final String why){
		try{
			action.run();
		}catch(IllegalArgumentException expected){
			checks++;
			return;
		}catch(Exception unexpected){throw new AssertionError("Unexpected fixture failure", unexpected);}
		throw new AssertionError(why);
	}
	/** Fails loudly on an unmet fixture contract; counts successful checks only. */
	private static void check(final boolean valid, final String why){
		if(!valid){throw new AssertionError(why);}
		checks++;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	/** Process-local count of successful assertions and expected policy rejections. */
	private static int checks;
}
