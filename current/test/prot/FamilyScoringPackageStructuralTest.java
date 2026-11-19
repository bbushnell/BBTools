package test.prot;

import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Objects;

import ddl.DDLFamilyProfiles;
import ddl.DDLPentamerBlacklist;
import prot.FamilyAlignmentMetric;
import prot.FamilyCalibration.VerificationMode;
import prot.FamilyReference;
import prot.FamilyReferenceCodec;
import prot.FamilyScoringPackage;
import prot.FamilyTailCascade;
import prot.FamilyTailCascade.Assignment;
import prot.ProteinSequence;
import structures.ByteBuilder;

/** File-backed structural verification is explicit and preserves classifier behavior.
 * Wrong inherited pins prove that nested readers do not silently hash; malformed
 * data and observable changes must still fail. Synthetic inputs only.
 * @author Brian Bushnell, Keqing
 */
public final class FamilyScoringPackageStructuralTest{

	/** Creates synthetic files in a new directory, then checks reloads, replacement profiles and input rejection. */
	public static void main(String[] args)throws Exception{
		check(args.length==1, "Supply a fresh Dori fixture directory");
		final Path root=Paths.get(args[0]).toAbsolutePath();
		Files.createDirectory(root);
		final byte[] alphabet="ACDEFGHIKLMNPQRSTVWY".getBytes(StandardCharsets.US_ASCII), query=new byte[64];
		for(int i=0; i<query.length; i++){query[i]=alphabet[(i*7+i/20)%alphabet.length];}
		final Path data=root.resolve("package");
		final char[][] profileRows=fixture(data, query);
		final Path config=data.resolve("classifier.config"), list=data.resolve("list.tsv");
		final FamilyScoringPackage initial=new FamilyScoringPackage(config, VerificationMode.STRUCTURAL);
		final FamilyTailCascade.LoadedClassifier baseline=initial.load(1);
		final Assignment expected=baseline.newWorker().assign("known", query);
		check(expected.assigned() && expected.familyId.equals("f1") && expected.ddlMatches==5 && expected.dpAttempts==1,
			"Literal five-hit family must reach the unchanged alignment/acceptance path");
		replacement(data, baseline, query, profileRows, expected);
		initial.verify();
		final FamilyScoringPackage structural=new FamilyScoringPackage(config, VerificationMode.STRUCTURAL);
		check(structural.verificationMode()==VerificationMode.STRUCTURAL &&
			structural.resolvedOptions().get("familyreferencessha80").equals(INHERITED),
			"Structural mode keeps the declared identity rather than inventing a fresh digest");
		final FamilyTailCascade.LoadedClassifier loaded=structural.load(1);
		same(expected, loaded.newWorker().assign("structural", query));
		final byte[] unknown=new byte[query.length];
		Arrays.fill(unknown, (byte)'X');
		same(baseline.newWorker().assign("unknown", unknown), loaded.newWorker().assign("unknown_structural", unknown));
		same(expected, structural.loadWithOverrides(1, 50, 5).newWorker().assign("override", query));
		same(baseline.withEvaluationMarginScale(.985).newWorker().assign("scaled", query),
			loaded.withEvaluationMarginScale(.985).newWorker().assign("scaled_structural", query));
		structural.verify();
		reject(()->structural.bindLengthRanking(loaded, null, null, null, null), "Length-model");
		// Every bound leaf matters, including payloads found through the profile list.
		for(String name:new String[]{"classifier.config", "references.tsv", "list.tsv", "profiles.tsv",
				"summary.tsv", "bounds.tsv", "medians.tsv", "hydro.tsv"}){
			final FamilyScoringPackage frozen=new FamilyScoringPackage(config, VerificationMode.STRUCTURAL);
			frozen.load(1);
			final Path file=data.resolve(name);
			final byte[] before=Files.readAllBytes(file);
			try{
				final byte[] changed=Arrays.copyOf(before, before.length+1);
				changed[before.length]='\n';
				Files.write(file, changed);
				reject(frozen::verify, file.toString());
			}finally{Files.write(file, before);}
		}
		final String inheritedConfig=text(config), inheritedSummary=text(data.resolve("summary.tsv"));
		write(config, inheritedConfig+"verification=structural\n");
		new FamilyScoringPackage(config, VerificationMode.STRUCTURAL).load(1);
		String unmeasuredConfig=inheritedConfig;
		for(String key:new String[]{"familyreferencessha80", "ddlprofilessha80", "scalarboundssha80",
				"lengthmedianssha80", "familyhydrosha80"}){
			unmeasuredConfig=unmeasuredConfig.replace(key+"="+INHERITED+"\n", "");
		}
		write(config, unmeasuredConfig);
		same(expected, new FamilyScoringPackage(config, VerificationMode.STRUCTURAL).load(1).newWorker().assign("unmeasured", query));
		write(config, inheritedConfig);
		write(config, inheritedConfig.replace("familyreferencessha80="+INHERITED, "familyreferencessha80=invalid"));
		reject(()->new FamilyScoringPackage(config, VerificationMode.STRUCTURAL), "Inherited package component");
		write(config, inheritedConfig+"blacklist="+data.resolve("absent.faa")+"\nblacklistsha80="+INHERITED+"\n");
		reject(()->new FamilyScoringPackage(config, VerificationMode.STRUCTURAL), "complete explicit blacklist policy");
		write(config, inheritedConfig);
		write(data.resolve("summary.tsv"), inheritedSummary.replace("tail\t5\t0\t", "tail\t6\t0\t"));
		reject(()->new FamilyScoringPackage(config, VerificationMode.STRUCTURAL).load(1), "recipe mismatch");
		write(data.resolve("summary.tsv"), inheritedSummary);
		final Path profiles=data.resolve("profiles.tsv");
		final String profileText=text(profiles);
		write(profiles, profileText+profileText.split("\n", -1)[1]+"\n");
		reject(()->new FamilyScoringPackage(config, VerificationMode.STRUCTURAL).load(1), "Duplicate profile");
		write(profiles, profileText);
		final FamilyScoringPackage beforeLoad=new FamilyScoringPackage(config, VerificationMode.STRUCTURAL);
		write(config, inheritedConfig+"# changed after binding\n");
		reject(()->beforeLoad.load(1), config.toString());
		write(config, inheritedConfig);
		new FamilyScoringPackage(config, VerificationMode.STRUCTURAL).load(1);
		System.out.println("FAMILY_SCORING_PACKAGE_STRUCTURAL_TEST_PASS checks="+checks+
			" structural_reload_parity=true inherited_pins=true transitive_stability=true default_hash_path_not_executed=true");
	}

	/** Literal three-family layout with two, five and zero hits; no digest API is called.
	 * The reference/scalar schema follows the accepted hit-floor fixture, but this
	 * builder writes inherited placeholder pins and never invokes its hashing loader. */
	private static char[][] fixture(Path root, byte[] raw)throws Exception{
		Files.createDirectory(root);
		final char[] sketch=new char[1024];
		DDLFamilyProfiles.sketchPositional(raw, 5, 0, 4, 2, sketch);
		final FamilyAlignmentMetric metric=new FamilyAlignmentMetric(FamilyAlignmentMetric.Mode.PIN_RIGHT_TRIM_LEFT_ALIGNED_QUERY, 28);
		final FamilyReferenceCodec.Header header=new FamilyReferenceCodec.Header(FamilyReferenceCodec.ALIGNMENT_SCHEMA,
			5, 1024, 0, 3, 3, 0, FamilyReferenceCodec.FORMULA_VERSION, FamilyReferenceCodec.SPECTRUM_ENCODING, 3.8906, metric);
		final ByteBuilder refs=FamilyReferenceCodec.appendHeader(header, new ByteBuilder());
		final ByteBuilder profiles=new ByteBuilder().append(DDLFamilyProfiles.HEADER16).nl();
		final int[] hits={2, 5, 0};
		final char[][] generated=new char[3][];
		for(int f=0; f<3; f++){
			final char[] row=sketch.clone();
			int available=0;
			generated[f]=row;
			for(int bin=0; bin<row.length; bin++){if(row[bin]!=0 && ++available>hits[f]){row[bin]=0;}}
			check(available>=hits[f], "Fixture needs enough distinct registers for the literal admission oracle");
			final String id="f"+f;
			final float[] spectrum=new float[400];
			spectrum[0]=1;
			FamilyReferenceCodec.appendFiveLines(new FamilyReference(id, f, 30, 1920, 1800, new ProteinSequence(id, raw), row, spectrum,
				0, 30, 0, 1.1f, 99, 1, raw, 0, new double[10], FamilyReference.AlignmentCalibration.ENABLED), refs, header);
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
		Files.write(root.resolve("references.tsv"), refs.toBytes());
		Files.write(root.resolve("profiles.tsv"), profiles.toBytes());
		write(root.resolve("summary.tsv"), "namespace\tk\tseed\texponent\tregisterbits\tsketch\tpartitions\toverlap\ttag_bits\tsegment_words\n"+
			"tail\t5\t0\t6\t16\t"+DDLFamilyProfiles.SKETCH_POSITIONAL+"\t4\t2\t"+
			DDLFamilyProfiles.positionalTagBits(4)+"\t"+DDLFamilyProfiles.positionalSegmentWords(4)+"\n");
		write(root.resolve("list.tsv"), root.resolve("profiles.tsv")+"\t"+INHERITED+"\t"+INHERITED+"\n");
		final StringBuilder scalar=new StringBuilder("scalar\trows\tmin\tmax\tmin_family\tmax_family\n");
		for(String name:new String[]{"hydrophobicity", "net_charge", "aromaticity", "fraction_small", "pI_tail"}){
			scalar.append(name).append("\t79565\t-20\t20\tlo\thi\n");
		}
		write(root.resolve("bounds.tsv"), scalar.toString());
		final StringBuilder medians=new StringBuilder("family\tmembers\tlength_p10\tlength_median\tlength_p90\tcount_under_200\n");
		final StringBuilder hydro=new StringBuilder("family_rank\tfamily_id\tN\thydrophobicity_mean\n");
		for(int f=0; f<3; f++){
			medians.append('f').append(f).append("\t30\t64\t64\t64\t30\n");
			hydro.append(f).append("\tf").append(f).append("\t30\t1.8\n");
		}
		write(root.resolve("medians.tsv"), medians.toString());
		write(root.resolve("hydro.tsv"), hydro.toString());
		final StringBuilder recipe=new StringBuilder("top=50\nk=5\nseed=0\nwidth=1024\nddlreject=5\nddlaccept=1025\n"+
			"positional=t\npartitions=4\noverlap=2\nindexearlyreject=f\ncosine=ranking\nrankingfeatures=combined25\n"+
			"scalargates=t\nscalarz=3.8906\nidentityacceptfloor=30\nblosumacceptfloor=0.145\nscalarscale=0.1\n"+
			"minmemberaa=40\nmaxmemberaa=99\nalignmentmode=pin_right_trim_left_aligned_query\nalignmentminspan=28\n"+
			"heuristic=hitfourth_h020_c150\nminqueryaa=40\nmaxqueryaa=99\n");
		for(String[] binding:new String[][]{{"familyreferences", "references.tsv"}, {"ddlprofiles", "list.tsv"}, {"scalarbounds", "bounds.tsv"},
			{"lengthmedians", "medians.tsv"}, {"familyhydro", "hydro.tsv"}}){
			recipe.append(binding[0]).append('=').append(root.resolve(binding[1])).append('\n')
				.append(binding[0]).append("sha80=").append(INHERITED).append('\n');
		}
		write(root.resolve("classifier.config"), recipe.toString());
		return generated;
	}

	/** A no-filter rebuild is an exact control; blocking every query word yields
	 * zero hits without changing the parent or its calibration/heuristic metadata. */
	private static void replacement(Path root, FamilyTailCascade.LoadedClassifier parent, byte[] query,
			char[][] rows, Assignment expected)throws Exception{
		final String[] ids={"f0", "f1", "f2"};
		final int[] ranks={0, 1, 2};
		check(parent.blacklist()==null && parent.matchesProfiles(rows), "Baseline parity inspects the actual indexed registers");
		final char[][] mutable=new char[rows.length][];
		for(int i=0; i<rows.length; i++){mutable[i]=rows[i].clone();}
		final FamilyTailCascade.LoadedClassifier control=parent.withReplacementProfiles(ids, ranks, mutable, null, 2);
		same(expected, control.newWorker().assign("same_rows", query));
		same(parent.withEvaluationMarginScale(.985).newWorker().assign("parent_scaled", query),
			parent.withEvaluationMarginScale(.985).withReplacementProfiles(ids, ranks, rows, null, 2).newWorker().assign("control_scaled", query));
		Arrays.fill(mutable[1], (char)0);
		check(!control.matchesProfiles(mutable) && control.matchesProfiles(rows), "Derived index must not borrow mutable caller rows");
		same(expected, control.newWorker().assign("after_mutation", query));
		reject(()->parent.withReplacementProfiles(new String[]{"f1", "f0", "f2"}, ranks, rows, null, 1), "identity/rank/width");
		reject(()->parent.withReplacementProfiles(ids, new int[]{0, 2, 1}, rows, null, 1), "identity/rank/width");
		final char[][] wrong=rows.clone();
		wrong[0]=new char[512];
		reject(()->parent.withReplacementProfiles(ids, ranks, wrong, null, 1), "identity/rank/width");
		final java.util.TreeSet<String> words=new java.util.TreeSet<String>();
		for(int i=0; i+5<=query.length; i++){words.add(new String(query, i, 5, StandardCharsets.US_ASCII));}
		final StringBuilder fasta=new StringBuilder();
		int ordinal=0;
		for(String word:words){fasta.append('>').append(++ordinal).append('\n').append(word).append('\n');}
		final Path file=root.resolve("all_query_words.faa");
		write(file, fasta.toString());
		final DDLPentamerBlacklist filter=DDLPentamerBlacklist.loadStructural(file, "unmeasured", words.size(), 5,
			DDLPentamerBlacklist.MatchMode.LITERAL, "synthetic_all_query_words");
		// Wiring-only negative control intentionally keeps unfiltered matching rows.
		// Empty reference rows alone cannot detect a missing query-side filter; this
		// control does not claim that those rows were produced with the filter.
		final FamilyTailCascade.LoadedClassifier queryControl=parent.withReplacementProfiles(ids, ranks, rows, filter, 2);
		check(queryControl.matchesProfiles(rows), "Query-filter control retains the known five-hit profile");
		final Assignment queryRejected=queryControl.newWorker().assign("query_filter_control", query);
		check(!queryRejected.assigned() && queryRejected.admittedFamilies==0 && queryRejected.dpAttempts==0 &&
			queryRejected.terminal.equals("ddl_reject"), "Query filtering must remove matches even when reference registers remain");
		final char[][] empty=new char[3][1024];
		final FamilyTailCascade.LoadedClassifier blocked=parent.withReplacementProfiles(ids, ranks, empty, filter, 2);
		check(blocked.blacklist()==filter && blocked.matchesProfiles(empty), "Query/profile policy retains the exact immutable filter object");
		final Assignment rejected=blocked.newWorker().assign("blocked", query);
		check(!rejected.assigned() && rejected.admittedFamilies==0 && rejected.dpAttempts==0 && rejected.terminal.equals("ddl_reject"),
			"Excluding every raw word rejects before alignment; it cannot fabricate recalibrated acceptance");
		check(parent.blacklist()==null && parent.matchesProfiles(rows), "Replacement view must preserve the unfiltered parent");
		same(expected, parent.newWorker().assign("parent_unchanged", query));
	}

	/** Compare all public assignment fields, including later-added observable counters. */
	private static void same(Assignment expected, Assignment actual)throws Exception{
		int fields=0;
		for(Field field:Assignment.class.getFields()){
			if(Modifier.isStatic(field.getModifiers())){continue;}
			check(Objects.deepEquals(field.get(expected), field.get(actual)), "Structural loading changed Assignment."+field.getName());
			fields++;
		}
		check(fields>0, "Assignment parity must inspect its public result fields");
	}
	/** Reads an ASCII fixture so a rejection case can restore its exact original text. */
	private static String text(Path path)throws Exception{return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);}
	/** Replaces one synthetic ASCII fixture; callers retain any text needed for restoration. */
	private static void write(Path path, String value)throws Exception{Files.write(path, value.getBytes(StandardCharsets.US_ASCII));}
	/** An operation expected to reject malformed or changed fixture inputs. */
	private interface Action{
		/** Executes the operation, allowing its checked diagnostic to reach the rejection helper. */
		void run()throws Exception;
	}
	/** Requires an exception with the specified message fragment; assertion failures propagate unchanged. */
	private static void reject(Action action, String reason)throws Exception{
		try{action.run();}catch(Exception failure){
			check(failure.getMessage()!=null && failure.getMessage().contains(reason), "Wrong rejection for "+reason+": "+failure);
			return;
		}
		throw new AssertionError("Malformed or changed input was accepted: "+reason);
	}
	/** Counts an attempted check and fails immediately with its fixture-specific explanation. */
	private static void check(boolean valid, String why){
		checks++;
		if(!valid){throw new AssertionError(why);}
	}
	/** Attempted checks in this process; this sequential fixture does not reset the counter. */
	private static int checks;
	/** Deliberately incorrect inherited identity, proving structural loading does not recompute a digest. */
	private static final String INHERITED="00000000000000000000";
}
