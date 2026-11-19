package test.prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import ddl.DDLFamilyHoldout;
import fileIO.ByteFile;
import prot.DigestSuffix;
import prot.FamilyCalibrate;
import prot.FamilyEval;
import prot.FamilyLengthScore;
import prot.FamilyScoringPackage;
import prot.FamilyTailCascade;
import shared.Shared;

/** Builds a literal two-family fixture and tests paired CLI and ordinary parity.
 * Owns the shared parent setup used by the input-only validation producer.
 * No private exporter, real panel or checkpoint is consumed.
 * Mutable assertion counts and amino-input settings belong to an isolated test JVM.
 * @author Sayu, Keqing, Fischl
 */
public final class FamilyPairedEvalTest{
	/*--------------------------------------------------------------*/
	/*----------------       Evaluation checks      ----------------*/
	/*--------------------------------------------------------------*/

	/** Builds a fresh parent and checks literal outcomes, CLI modes and worker parity.
	 * Creates files only beneath the supplied fresh absolute root. Construction
	 * enables the process-wide amino-input setting; use an isolated test JVM.
	 */
	public static void main(String[] args) throws Exception{
		check(args.length==1, "Require a fresh absolute fixture directory");
		final Path fixture=Paths.get(args[0]).normalize();
		check(fixture.isAbsolute() && !Files.exists(fixture) && !Files.isSymbolicLink(fixture),
			"Paired fixture root must be fresh, absolute and nonsymlinked");
		Files.createDirectory(fixture);
		prepareFixture(fixture);
		final Path root=fixture.resolve("paired_eval");
		final Inputs inputs=writeDefaultHashInputs(fixture, root);
		final Path config=inputs.config, model=inputs.model, medians=inputs.medians, panel=inputs.panel;
		final byte[][] sequences=inputs.sequences;
		final FamilyScoringPackage source=new FamilyScoringPackage(config);
		final FamilyTailCascade.LoadedClassifier base=source.load(2);
		final FamilyLengthScore binding=source.bindLengthRanking(base, model, hash(model), medians, hash(medians));
		final FamilyLengthScore.Mode[] modes={
			FamilyLengthScore.Mode.OFF, FamilyLengthScore.Mode.MULTIPLY, FamilyLengthScore.Mode.ADD};
		final long[] expected=new long[3];
		for(int a=0; a<3; a++){
			final FamilyTailCascade.Worker worker=base.withLengthRanking(binding, modes[a]).newWorker();
			for(int i=0; i<sequences.length; i++){
				final FamilyTailCascade.Assignment result=worker.assign("q"+i, sequences[i]);
				check(result.assigned()==(i%3!=0),
					"Literal fixture outcome differs: arm="+a+" query="+i+" terminal="+result.terminal);
				if(result.assigned()){
					expected[a]++;
					check(result.familyId.equals("f0"), "Both rankings recover the callable f0 family");
					check(result.dpAttempt==(a==0 ? 1 : 2),
						"Learned family-fraction ranking must try rejected f1 before accepted f0: arm="+a+" dp="+result.dpAttempt);
				}
			}
			check(expected[a]==8, "Each arm must retain eight acceptances and four rejections");
		}
		check(expected[0]>0 && expected[0]<12, "Fixture must exercise accepted and rejected queries");
		for(boolean positive:new boolean[]{false, true}){
			for(int threads:new int[]{1, 4}){
				final Path out=root.resolve((positive ? "positive" : "negative")+threads);
				FamilyEval.main(arguments(config, model, medians, panel, out, positive, threads));
				final List<String> rows=Files.readAllLines(out.resolve("summary.tsv"), StandardCharsets.US_ASCII);
				check(rows.size()==4, "Paired CLI writes three arm rows for one panel");
				checkSummary(rows);
				for(int a=0; a<3; a++){
					final String[] row=rows.get(a+1).split("\t");
					check(Long.parseLong(row[4])==12 && Long.parseLong(row[5])==expected[a],
						"Paired counts match independent single-arm scoring");
				}
				check(Files.readAllLines(out.resolve("oracle_membership.tsv")).size()==9,
					"All eight membership masks are retained");
				check(Files.readAllLines(out.resolve("paired_transitions.tsv")).size()==3,
					"Both ON comparisons have paired outcomes");
				check(Files.readAllLines(out.resolve("PASS")).get(0).contains("queries=12 arms=3"),
					"Completion counts every query once, not three times");
			}
			for(String name:new String[]{"summary.tsv", "paired_transitions.tsv",
				"oracle_membership.tsv", "family_false_accepts.tsv.gz"}){
				check(lines(root.resolve((positive ? "positive" : "negative")+1).resolve(name))
					.equals(lines(root.resolve((positive ? "positive" : "negative")+4).resolve(name))),
					"Worker count preserves "+name);
			}
		}
		final Path ordinary=root.resolve("ordinary");
		FamilyEval.main(new String[]{"in="+config, "neg="+panel, "out="+ordinary, "t=1"});
		check(Long.parseLong(Files.readAllLines(ordinary.resolve("summary.tsv")).get(1).split("\t")[3])==expected[0],
			"Default OFF result matches ordinary FamilyEval");
		final String[] good=arguments(config, model, medians, panel, root.resolve("invalid"), false, 1);
		FamilyEval.main(append(replace(good, "out", root.resolve("default_structural").toString()), "negcounts=12"));
		check(Files.readAllLines(root.resolve("default_structural/evaluation_policy.tsv"), StandardCharsets.US_ASCII)
			.contains("verification\tstructural"), "Omitted verification uses structural mode with inherited labels");
		FamilyEval.main(structural(good, root.resolve("structural"), 12));
		check(Files.readAllLines(root.resolve("structural/evaluation_policy.tsv"), StandardCharsets.US_ASCII)
			.contains("identity_semantics\tinherited_unverified"),
			"Structural paired CLI reports inherited, unverified identities");
		FamilyEval changed=new FamilyEval(structural(good, root.resolve("changed"), 12));
		Files.setLastModifiedTime(model,
			java.nio.file.attribute.FileTime.fromMillis(Files.getLastModifiedTime(model).toMillis()+5000));
		rejectAction(changed::process, "Structural input changed");
		reject(append(good, "pos="+panel), "Paired scoring cannot mix selection negatives and positive evaluation");
		reject(structural(good, root.resolve("wrong_count"), 11), "Structural paired counts must match observations");
		reject(append(good, "marginscale=1"), "Ambiguous global and per-arm margins fail");
		reject(append(good, "pairedscales=1,NaN,1"), "Nonfinite paired scale fails");
		reject(replace(append(good, "verification=sha80"), "lengthmodelsha80", hash(medians)), "Wrong model pin fails");
		check(!Files.exists(root.resolve("invalid/PASS")), "Rejected policy does not publish PASS");
		System.out.println("FAMILY_PAIRED_EVAL_TEST_PASS checks="+checks);
	}

	/*--------------------------------------------------------------*/
	/*----------------       Shared fixture         ----------------*/
	/*--------------------------------------------------------------*/

	/** Writes and calibrates the one shared five-member parent under an existing root.
	 * Both inputs/ and package/ must be absent. The returned paths and arguments
	 * describe this actual construction, for the input-only producer's manifest.
	 * Like the standalone tests, this sets the current JVM's amino-input mode.
	 */
	static ParentFixture prepareFixture(Path root) throws Exception{
		if(!root.isAbsolute() || !Files.isDirectory(root) || Files.isSymbolicLink(root)){
			throw new IllegalArgumentException("Fixture parent must be an existing absolute nonsymlink directory: "+root);
		}
		final Path inputs=root.resolve("inputs"), output=root.resolve("package");
		if(Files.exists(inputs) || Files.isSymbolicLink(inputs) ||
			Files.exists(output) || Files.isSymbolicLink(output)){
			throw new IllegalArgumentException("Fixture inputs and package must be fresh: "+root);
		}
		final ParentFixture fixture=new ParentFixture(inputs, output);
		Shared.AMINO_IN=true;
		FamilyCalibrate.main(fixture.constructionArgs);
		return fixture;
	}

	/** Literal population shared by paired evaluation and default-hash validation.
	 * The three genome labels come from the retained five-member fixture. They
	 * are data, not a request to run the historical exporter's split algorithm.
	 */
	static final class ParentFixture{
		/** Creates the exact two input chunks, bounds and preparation recipe. */
		private ParentFixture(Path root, Path output) throws Exception{
			Files.createDirectory(root);
			final Path families=root.resolve("families.tsv"), chunks=root.resolve("chunks.tsv");
			final Path consensus=root.resolve("consensus.faa"), bounds=root.resolve("bounds.tsv");
			final Path negative=root.resolve("negative.faa"), panels=root.resolve("panels.tsv");
			manifest=root.resolve("members.config");
			lengths=root.resolve("lengths.tsv");
			final String a=new String(repeat('A', 100), StandardCharsets.US_ASCII);
			final String a200=new String(repeat('A', 200), StandardCharsets.US_ASCII);
			write(families, DDLFamilyHoldout.FAMILY_HEADER+"\n0\tf0\t3\n1\tf1\t2\n");
			write(consensus, ">f0\n"+a+"\n>f1\n"+a+"\n");
			final String[] members={
				"0\tf0\tgenome0|a\tgenome0\t0\t"+a+"\n1\tf1\tgenome0|b\tgenome0\t1\t"+a+
					"\n0\tf0\tgenome0|long\tgenome0\t2\t"+a200+"\n",
				"0\tf0\tgenome2|a\tgenome2\t0\t"+a+"\n1\tf1\tgenome4|a\tgenome4\t0\t"+a+"\n"};
			final StringBuilder parts=new StringBuilder(DDLFamilyHoldout.CHUNK_HEADER).append('\n');
			for(int i=0; i<members.length; i++){
				final Path part=root.resolve("part"+i+".tsv");
				write(part, DDLFamilyHoldout.MEMBER_HEADER+"\n"+members[i]);
				parts.append("b0").append(i).append("\t0\t").append(part).append('\t').append(hash(part))
					.append(i==0 ? "\t3\t400\n" : "\t2\t200\n");
			}
			write(chunks, parts.toString());
			write(negative, ">fixture_control\n"+new String(repeat('C', 100), StandardCharsets.US_ASCII)+"\n");
			write(panels, "panel\tpath\tsha80\tqueries\nk5\t"+negative+"\t"+hash(negative)+"\t1\n");
			final StringBuilder scalar=new StringBuilder("scalar\trows\tmin\tmax\tmin_family\tmax_family\n");
			for(String name:new String[]{"hydrophobicity", "net_charge", "aromaticity", "fraction_small", "pI_tail"}){
				scalar.append(name).append("\t79565\t-20\t20\tlo\thi\n");
			}
			write(bounds, scalar.toString());
			write(lengths, "family\tmembers\tlength_p10\tlength_median\tlength_p90\tcount_under_200\n"+
				"f1\t2\t100\t100\t100\t2\nf0\t3\t100\t100\t200\t2\n");
			write(manifest, "schema=family_preparation_v1\nfamilies="+families+"\nfamiliessha80="+hash(families)+
				"\nchunks="+chunks+"\nchunkssha80="+hash(chunks)+"\nconsensus="+consensus+
				"\nconsensussha80="+hash(consensus)+"\nnegativepanels="+panels+"\nnegativepanelssha80="+hash(panels)+
				"\ncosinepanel=k5\nexpectedfamilies=2\nexpectedmembers=5\nbuckets=1\nbatches=2\nmaxreferenceaa=1000\n"+
				"maxbucketmembers=16\nmaxbucketbytes=16384\nmaxpreparedbytes=1048576\nmaxoutputbytes=1048576\n"+
				"retention=metric_specific\nuppertop=1\nminnegativeaa=100\nmaxnegativeaa=200\nk=5\nseed=0\nwidth=1024\n"+
				"matrix=blosum62\npositional=t\npartitions=4\noverlap=2\nscalarbounds="+bounds+
				"\nscalarboundssha80="+hash(bounds)+"\nscalarscale=0.1\nblosumfloor=0\nscalarz=3.8906\n"+
				"monomer=t\ntop=50\nddlreject=8\nddlaccept=1025\ncosine=ranking\nrankingfeatures=combined25\n"+
				"scalargates=t\nidentityacceptfloor=30\nblosumacceptfloor=0\nindexearlyreject=f\n");
			constructionArgs=new String[]{"prepare="+manifest, "out="+output, "verification=sha80",
				"t=2", "rounds=1", "targetfpr=1", "perclass=t"};
		}

		/** Literal input recipe and unmodified source lengths retained by the manifest. */
		final Path manifest, lengths;
		/** Actual seven arguments passed to FamilyCalibrate, owned by this fixture. */
		final String[] constructionArgs;
	}

	/** The retrieval fixture has only two matching registers; lower its TEST-ONLY reject gate explicitly.
	 * Resolve relative components against the original recipe, preserving their pins and original files.
	 */
	static Inputs writeDefaultHashInputs(Path fixture, Path root) throws Exception{
		Files.createDirectory(root);
		final Path config=derivedConfig(fixture.resolve("package/classifier.config"), root.resolve("classifier.config"));
		final Path model=root.resolve("family_fraction.bbnet"), medians=root.resolve("medians.tsv");
		// f0 has four occupied profile registers, f1 two; an A query matches two in either.
		// Increasing family-fraction weight makes ON try f1 first, then recover through f0.
		write(model,
			"##bbnet\n#version 1\n#concise\n#dense\n#density 1\n#blocksize 1\n#seed 1\n#layers 2\n#dims 3 1\n#coding decimal\nC4 SIG 0 0 4 0\n");
		// Keep derived-test medians separate from the parent fixture's source lengths.
		write(medians,
			"family\tmembers\tlength_p10\tlength_median\tlength_p90\tcount_under_200\nf1\t2\t100\t100\t100\t2\nf0\t3\t100\t100\t200\t2\n");
		final Path panel=root.resolve("queries.faa");
		final StringBuilder text=new StringBuilder();
		final byte[][] sequences=new byte[12][];
		for(int i=0; i<12; i++){
			sequences[i]=repeat(i%3==0 ? 'C' : 'A', 100+i%2);
			text.append(">q").append(i).append('\n')
				.append(new String(sequences[i], StandardCharsets.US_ASCII)).append('\n');
		}
		write(panel, text.toString());
		return new Inputs(config, model, medians, panel, sequences);
	}
	/** Paths and query arrays shared with the input-only producer; not defensive copies. */
	static final class Inputs{
		/** Retains the supplied paths and array references without further file access. */
		Inputs(Path config, Path model, Path medians, Path panel, byte[][] sequences){
			this.config=config;
			this.model=model;
			this.medians=medians;
			this.panel=panel;
			this.sequences=sequences;
		}
		/** Derived classifier, length model, median table and twelve-query FASTA paths. */
		final Path config, model, medians, panel;
		/** Query bytes in FASTA order, used by the independent single-arm assertions. */
		final byte[][] sequences;
	}
	/** Copies a parent recipe, resolves its resources and changes only the test DDL reject gate. */
	private static Path derivedConfig(Path original, Path out) throws Exception{
		final StringBuilder recipe=new StringBuilder("# paired-test-only source=")
			.append(original).append(" sha80=").append(hash(original)).append('\n');
		int changed=0;
		for(String line:Files.readAllLines(original, StandardCharsets.US_ASCII)){
			if(line.isEmpty() || line.startsWith("#")){continue;}
			final String[] pair=line.split("=", 2);
			check(pair.length==2, "Fixture recipe must contain key=value");
			if(pair[0].equals("ddlreject")){
				check(pair[1].equals("8"), "Expected the original fixture's reject8 gate");
				pair[1]="1";
				changed++;
			}
			if(Arrays.asList("familyreferences", "ddlprofiles", "scalarbounds", "blacklist").contains(pair[0])){
				pair[1]=original.getParent().resolve(pair[1]).normalize().toString();
			}
			recipe.append(pair[0]).append('=').append(pair[1]).append('\n');
		}
		check(changed==1, "Exactly one test-only DDL gate must change");
		write(out, recipe.toString());
		return out;
	}

	/*--------------------------------------------------------------*/
	/*----------------       CLI test utilities     ----------------*/
	/*--------------------------------------------------------------*/

	/** Builds one paired CLI invocation with explicit model and median bindings. */
	private static String[] arguments(Path config, Path model, Path medians, Path panel,
		Path out, boolean pos, int threads){
		return new String[]{"in="+config, "paired=t", "lengthmodel="+model, "lengthmodelsha80="+hash(model),
			"lengthmedians="+medians,
			"lengthmedianssha80="+hash(medians), "pairedscales=1,1,1", (pos ? "pos=" : "neg=")+panel,
			"out="+out, "t="+threads};
	}
	/** Copies arguments for structural verification with an explicit expected query count. */
	private static String[] structural(String[] args, Path out, int count){
		final ArrayList<String> result=new ArrayList<String>();
		for(String arg:args){
			if(arg.startsWith("out=")){
				result.add("out="+out);
			}else if(arg.startsWith("lengthmodelsha80=") || arg.startsWith("lengthmedianssha80=")){
				result.add(arg.replaceFirst("=.*", "=unmeasured"));
			}else{result.add(arg);}
		}
		result.add("verification=structural");
		result.add("negsha80=unmeasured");
		result.add("negcounts="+count);
		return result.toArray(new String[0]);
	}
	/** Returns a new argument array containing the original values followed by extra values. */
	private static String[] append(String[] values, String... extra){
		final String[] result=Arrays.copyOf(values, values.length+extra.length);
		System.arraycopy(extra, 0, result, values.length, extra.length);
		return result;
	}
	/** Copies arguments and replaces the first exact key prefix; fails if the key is absent. */
	private static String[] replace(String[] values, String key, String value){
		final String[] result=values.clone();
		final String prefix=key+"=";
		for(int i=0; i<result.length; i++){
			if(result[i].startsWith(prefix)){
				result[i]=prefix+value;
				return result;
			}
		}
		throw new AssertionError(key);
	}
	/** Requires the declared stage columns and conservation of queries across all stages. */
	private static void checkSummary(List<String> rows){
		final String[] header=rows.get(0).split("\t");
		check(header.length==11+METRIC_STAGES.length, "Paired summary labels every metric stage");
		for(int i=0; i<METRIC_STAGES.length; i++){
			check(header[11+i].equals(METRIC_STAGES[i]), "Paired metric stage label "+i+" differs");
		}
		for(int r=1; r<rows.size(); r++){
			final String[] row=rows.get(r).split("\t");
			long stages=0;
			check(row.length==header.length, "Paired summary row/header field count differs");
			for(int i=0; i<METRIC_STAGES.length; i++){stages+=Long.parseLong(row[11+i]);}
			check(stages==Long.parseLong(row[4]), "Paired metric stages do not conserve queries");
		}
	}
	/** Requires a CLI policy rejection; other exceptions remain visible to the caller. */
	private static void reject(String[] args, String why) throws Exception{
		boolean failed=false;
		try{
			FamilyEval.main(args);
		}catch(IllegalArgumentException e){failed=true;}
		check(failed, why);
	}
	/** A deferred operation for a structural-input rejection check. */
	private interface Action{
		/** Executes the operation without suppressing checked failures. */
		void run()throws Exception;
	}
	/** Requires an argument or structural-state rejection from the deferred operation. */
	private static void rejectAction(Action action, String why)throws Exception{
		boolean failed=false;
		try{
			action.run();
		}catch(IllegalArgumentException|IllegalStateException e){failed=true;}
		check(failed, why);
	}
	/** Decodes text through ByteFile so compressed and plain outputs compare by content. */
	private static List<String> lines(Path p){
		final ArrayList<String> out=new ArrayList<String>();
		for(byte[] row:ByteFile.toLines(p.toString())){
			out.add(new String(row, StandardCharsets.US_ASCII));
		}
		return out;
	}
	/** Returns a new constant-residue sequence of the requested length. */
	private static byte[] repeat(char c, int n){
		final byte[] b=new byte[n];
		Arrays.fill(b, (byte)c);
		return b;
	}
	/** Returns the digest suffix consumed by explicit fixture binding checks. */
	private static String hash(Path p){return DigestSuffix.file(p.toString());}
	/** Creates an ASCII fixture file, refusing to overwrite an existing path. */
	private static void write(Path p, String text) throws Exception{
		Files.write(p, text.getBytes(StandardCharsets.US_ASCII), java.nio.file.StandardOpenOption.CREATE_NEW);
	}
	/** Counts an attempted assertion and fails loudly if its contract is not met. */
	private static void check(boolean valid, String why){
		checks++;
		if(!valid){throw new AssertionError(why);}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	/** Process-local assertion count; standalone execution uses a fresh JVM. */
	private static int checks;
	/** Exact summary stage order after the eleven leading metadata/count columns. */
	private static final String[] METRIC_STAGES={"ddl_reject", "ddl_accept", "identity_accept", "blosum_accept",
		"identity_reject", "blosum_reject", "no_hit", "scalar_reject", "alignment_unavailable",
		"input_too_short", "no_alignment", "span_too_short"};
}
