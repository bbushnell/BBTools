package benchmark;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Locale;

import benchmark.CallGenesBenchmarkSupport.InputVerification;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;

import parse.LineParser1;
import parse.Parse;
import parse.PreParser;

import shared.Shared;

import structures.ByteBuilder;

/** Extracts scorer metadata from saved paired-caller tables; never reads sequence payloads.
 * @author Sayu
 */
public final class CallGenesSavedMetadata{
	/** Parse explicit saved-table bindings; close the diagnostic stream on every exit. */
	public static void main(String[] args)throws Exception{
		final PreParser pp=new PreParser(args, CallGenesSavedMetadata.class, false);
		try{
			final HashMap<String, String> options=new HashMap<String, String>();
			for(String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				require(pair.length==2 && !pair[1].isEmpty() && options.put(pair[0].toLowerCase(Locale.ROOT), pair[1])==null,
					"Require unique key=value arguments");
			}
			require(options.size()==7 && options.keySet().containsAll(Arrays.asList(
				"accepted", "raw", "retained", "summary", "out", "policy", "expectedgenomes")),
				"Require accepted/raw/retained/summary/out/policy/expectedgenomes");
			Shared.setThreads(1);
			extract(Paths.get(options.get("accepted")), Paths.get(options.get("raw")), Paths.get(options.get("retained")),
			Paths.get(options.get("summary")), Paths.get(options.get("out")), options.get("policy"),
			Integer.parseInt(options.get("expectedgenomes")));
		}finally{
			Shared.closeStream(pp.outstream);
		}
	}

	/** The cohort, aliases and row order are the saved paired-caller ABI. */
	public static void extract(Path accepted, Path raw, Path retained, Path summary, Path output,
	String policy, int expectedGenomes)throws Exception{
		require(expectedGenomes>0, "expectedgenomes must be positive");
		require(policy!=null && policy.matches("[A-Za-z0-9][A-Za-z0-9_.-]*"), "Unsafe policy prefix");
		require(output.isAbsolute() && !Files.exists(output, LinkOption.NOFOLLOW_LINKS), "Fresh absolute output directory required");
		final Path out=output.normalize().getParent().toRealPath().resolve(output.getFileName());
		require(!Files.exists(out, LinkOption.NOFOLLOW_LINKS), "Output already exists");
		final InputVerification stability=new InputVerification("structural");
		for(Path input:new Path[]{
			accepted, raw, retained, summary
		}){
			stability.pin(input, "unmeasured");
		}
		final LinkedHashMap<String, Input> roster=new LinkedHashMap<String, Input>();
		read(accepted, CallGenesBenchmarkSupport.INPUT_HEADER, 13, p->{
			require(p.length(0)>0, "Empty accepted genome");
			final Input in=new Input(p.parseByteArray(0), p.parseByteArray(1), p.parseByteArray(3));
			require(roster.put(new String(in.genome, StandardCharsets.UTF_8), in)==null, "Duplicate accepted genome");
		});
		require(roster.size()==expectedGenomes, "Accepted cohort size differs from expectedgenomes");
		read(retained, CallGenesBenchmarkSupport.RETAINED_GFF_HEADER, 3, p->{
			final Input in=roster.get(text(p, 0));
			require(in!=null && in.nativeCalls==null, "Foreign or duplicate retained genome");
			require(Paths.get(text(p, 1)).isAbsolute() && Paths.get(text(p, 2)).isAbsolute(), "Retained calls must be absolute paths");
			in.nativeCalls=p.parseByteArray(1);
			in.formulaCalls=p.parseByteArray(2);
		});
		for(Input in:roster.values()){
			require(in.nativeCalls!=null, "Retained inventory differs from accepted cohort");
		}
		final ByteBuilder grades=new ByteBuilder("genome\tarm\t").append(String.join("\t", CallGenesBenchmarkSupport.GRADE_NAMES)).nl();
		final String scoreHeader=CallGenesBenchmarkSupport.RAW_SCORE_HEADER.trim();
		read(raw, scoreHeader, 19, p->{
			final int arm=arm(p);
			if(arm<0 || !p.termEquals("ALL", 2)){
				return;
			}
			final Input in=roster.get(text(p, 1));
			require(in!=null, "Raw grades name a foreign genome");
			final int bit=1<<arm;
			require((in.rawSeen&bit)==0, "Duplicate ALL-band grade");
			in.rawSeen|=bit;
			grades.appendTerm(p, 1).tab().append(alias(arm));
			for(int i=11; i<19; i++){
				nonnegative(p, i);
				grades.tab().appendTerm(p, i);
			}
			grades.nl();
		});
		final ByteBuilder inputs=new ByteBuilder(CallGenesBenchmarkSupport.SAVED_SCORE_INPUT_HEADER).nl();
		final byte[] gradesPath=out.resolve("grades.tsv").toString().getBytes(StandardCharsets.UTF_8);
		for(Input in:roster.values()){
			require(in.rawSeen==3, "Missing paired ALL-band grade");
			for(int arm=0; arm<2; arm++){
				inputs.append(policy).append('_').append(arm==0?"native":"formula").tab().append(in.genome).tab()
				.append(in.fasta).tab().append(in.reference).tab().append(arm==0?in.nativeCalls:in.formulaCalls).tab()
				.append(gradesPath).tab().append(alias(arm)).nl();
			}
		}
		final boolean[] controls=new boolean[2];
		final ByteBuilder expected=new ByteBuilder(CallGenesBenchmarkSupport.SAVED_SCORE_EXPECTED_HEADER).nl();
		read(summary, scoreHeader, 19, p->{
			final int arm=arm(p);
			if(arm<0 || !p.termEquals("ALL", 1) || !p.termEquals("ALL", 2)){
				return;
			}
			require(!controls[arm], "Duplicate ALL/ALL summary policy");
			controls[arm]=true;
			expected.append(policy).append('_').appendTerm(p, 0).tab().append(expectedGenomes);
			for(int i=3; i<9; i++){
				nonnegative(p, i);
				expected.tab().appendTerm(p, i);
			}
			expected.nl();
		});
		require(controls[0] && controls[1], "Missing paired ALL/ALL summary controls");
		stability.verifyComplete();
		Files.createDirectory(out);
		write(out.resolve("grades.tsv"), grades);
		write(out.resolve("inputs.tsv"), inputs);
		write(out.resolve("expected.tsv"), expected);
		stability.verifyComplete();
		write(out.resolve("PASS"), new ByteBuilder("SAVED_CALLER_METADATA_PASS genomes=").append(expectedGenomes)
		.append(" policies=2 caller_runs=0 classifier_runs=0\n"));
	}

	/** Own selected field copies from the borrowed parser; rawSeen is the two-policy presence mask. */
	private static final class Input{
		Input(byte[] genome_, byte[] fasta_, byte[] reference_){
			genome=genome_;
			fasta=fasta_;
			reference=reference_;
		}
		final byte[] genome, fasta, reference;
		byte[] nativeCalls, formulaCalls;
		int rawSeen;
	}

	/** The parser and its input line are borrowed only for this callback. */
	private interface Row{
		void accept(LineParser1 parser)throws Exception;
	}

	/** Check the exact ABI and reuse one byte parser; retained fields must be copied by the callback. */
	private static void read(Path path, String header, int fields, Row consumer)throws Exception{
		final ByteFile in=ByteFile.makeByteFile(path.toString(), false);
		final LineParser1 parser=new LineParser1('\t');
		try{
			require(Arrays.equals(in.nextLine(), header.getBytes(StandardCharsets.UTF_8)), "Unexpected header: "+path);
			for(byte[] line; (line=in.nextLine())!=null; ){
				parser.set(line);
				require(parser.terms()==fields, "Unexpected column count: "+path);
				consumer.accept(parser);
			}
		}finally{
			require(!in.close(), "Failed reading metadata: "+path);
		}
	}

	/** Decode only a selected lookup/path field as UTF-8, leaving emitted fields as original bytes. */
	private static String text(LineParser1 parser, int term){
		parser.setBounds(term);
		return new String(parser.line(), parser.a(), parser.b()-parser.a(), StandardCharsets.UTF_8);
	}

	private static int arm(LineParser1 parser){
		return parser.termEquals("native", 0)?0:parser.termEquals("formula", 0)?1:-1;
	}

	private static String alias(int arm){
		return arm==0?"native":"bonus400";
	}

	/** Validate a borrowed unsigned-decimal span with checked arithmetic; preserve its spelling. */
	private static void nonnegative(LineParser1 parser, int term){
		parser.setBounds(term);
		require(parser.a()<parser.b(), "Empty count");
		long value=0;
		for(int i=parser.a(); i<parser.b(); i++){
			final int digit=parser.line()[i]-'0';
			require(digit>=0 && digit<=9 && value<=(Long.MAX_VALUE-digit)/10, "Invalid or overflowing nonnegative count");
			value=value*10+digit;
		}
	}

	/** Publish one fresh byte table and wait for its writer to close before reporting success. */
	private static void write(Path path, ByteBuilder bytes){
		final ByteStreamWriter out=new ByteStreamWriter(path.toString(), false, false, false);
		out.start();
		out.print(bytes);
		require(!out.poisonAndWait(), "Metadata output write failed: "+path);
	}

	private static void require(boolean valid, String message){
		if(!valid){
			throw new IllegalArgumentException(message);
		}
	}
}
