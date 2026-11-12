package prok;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Locale;

import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaReferenceGff.Row;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/** Materializes explicit reviewed disposition claims without making biological
 * decisions. Derived from A_SampleByteFile; native parsing, reusable buffers and
 * checked stream closure. Rule approval/authority fields are required declarations,
 * not independent verification of human approval or biological truth.
 *
 * All decision rows survive in a companion. Only resolved scored/dont_care rows
 * enter GFFs; no threshold, model ranking, merging or coordinate repair occurs.
 * Each input gets a GFF, including empty views. A completed reference.tsv is
 * atomically published last; a failed partial directory has no completed manifest.
 * Single-threaded metadata processing; no FASTA, model or caller execution.
 * Optional rebase mode preserves existing line content while replacing one
 * declared evidence pointer; witness copying and its provenance are external.
 * @author Brian Bushnell, Raiden
 */
public final class WriteNcrnaReference{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final PreParser pp=new PreParser(args, WriteNcrnaReference.class, false);
		try{new WriteNcrnaReference(pp.args, pp.outstream).process(timer);}
		finally{Shared.closeStream(pp.outstream);}
	}

	private WriteNcrnaReference(final String[] args, final PrintStream stream){
		outstream=stream;
		Shared.setThreads(1);
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("mode")){
				require("write".equals(b) || "rebase".equals(b), "mode must be write or rebase");
				rebase=b.equals("rebase");
			}else if(a.equals("evidencefrom")){
				evidenceFrom=b;
			}else if(a.equals("evidenceto")){
				evidenceTo=b;
			}else if(a.equals("inputs")){
				inputsPath=b;
			}else if(a.equals("contigs")){
				contigsPath=b;
			}else if(a.equals("rules")){
				rulesPath=b;
			}else if(!parser.parse(arg, a, b)){
				throw new IllegalArgumentException("Unknown parameter: "+arg);
			}
		}
		decisionsPath=parser.in1;
		require(decisionsPath!=null && parser.out1!=null, "Supply in= and out=new_directory");
		require(!parser.overwrite && !parser.append, "Reference materialization requires a fresh directory; overwrite/append unsupported");
		if(rebase){
			require(inputsPath==null && contigsPath==null && rulesPath==null && evidenceFrom!=null && evidenceTo!=null,
				"Rebase requires in=reference_manifest evidencefrom=old_pointer evidenceto=existing_witness out=new_directory");
			require(!evidenceFrom.equals(evidenceTo) && new File(evidenceFrom).isAbsolute() && new File(evidenceTo).isAbsolute()
				&& new File(evidenceTo).isFile(), "Rebase pointers must be distinct absolute paths; the destination witness must already exist");
			require(parser.in2==null && parser.out2==null && parser.maxReads<0 && parser.samplerate==1,
				"Rebasing must retain the complete unpaired reference manifest");
			require(Tools.testInputFiles(false, true, decisionsPath, evidenceTo), "Rebase manifest and destination witness must be readable");
		}else{
			require(inputsPath!=null && contigsPath!=null && rulesPath!=null && evidenceFrom==null && evidenceTo==null,
				"Supply in=dispositions inputs=roster contigs=observations rules=approved_rules out=new_directory");
			require(Tools.testInputFiles(false, true, decisionsPath, inputsPath, contigsPath, rulesPath), "Materializer inputs must be readable");
		}
		output=new File(parser.out1).getAbsoluteFile();
		require(rawSafe(output.getPath()), "Output directory must contain no tab/CR/LF for the literal reference manifest");
		// File.mkdir can replace unrepresentable names under an ASCII filesystem
		// locale. Validate the native path before any output side effect.
		outputPath=output.toPath();
		require(!output.exists(), "Reference output directory must not exist: "+output);
		companion=new File(output, "dispositions.tsv").getPath();
	}

	/*--------------------------------------------------------------*/
	/*----------------           Metadata           ----------------*/
	/*--------------------------------------------------------------*/

	private void loadRules(){
		try(Input in=new Input(rulesPath)){
			RrnaResourceIO.header(in.reader, RULE_HEADER, rulesPath);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				parser.set(line);
				require(parser.terms()==5, "Rule rows require five fields");
				final String id=decoded(0), family=decoded(1), status=decoded(2), approval=decoded(3), authority=decoded(4);
				require(present(id) && !family.isEmpty() && present(authority), "Rules need explicit IDs, family fields and authority pointers");
				require(approval.equals("APPROVED"), "Rule table must record APPROVED authority; a scored locus alone is not approval: "+id);
				final int kind=kind(status);
				require(kind>=0, "Unknown rule disposition: "+status);
				require(kind>MASK || present(family), "A scored/masked rule requires a named family");
				require(!rules.contains(id), "Duplicate rule ID: "+id);
				rules.put(id, new Rule(id, family, kind));
				ruleRows.add(line);
			}
		}
		require(rules.size()>0, "An explicit approved rule table is required");
	}

	private void loadInputs(){
		final ObjectSet<String> fastas=new ObjectSet<String>(String.class);
		try(Input in=new Input(inputsPath)){
			RrnaResourceIO.header(in.reader, INPUT_HEADER, inputsPath);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				parser.set(line);
				require(parser.terms()==2, "Input roster requires input_id and literal FASTA path");
				final String id=literal(0), fasta=literal(1);
				require(present(id) && new File(fasta).isAbsolute() && rawSafe(id) && rawSafe(fasta),
					"Roster needs nonempty IDs and absolute paths representable in the literal UTF-8 manifest");
				require(!inputs.contains(id) && fastas.add(fasta), "Input IDs and exact FASTA paths must be unique");
				final Genome genome=new Genome(id, fasta, new File(output, "input_"+(order.size()+1)+".gff.gz").getPath());
				inputs.put(id, genome);
				order.add(genome);
			}
		}
		require(!order.isEmpty(), "Reference roster must contain at least one input");
	}

	/** Saved contig observations use the evidence producer's escaped TSV framing. */
	private void loadContigs(){
		Genome genome=null;
		try(Input in=new Input(contigsPath)){
			RrnaResourceIO.header(in.reader, CONTIG_HEADER, contigsPath);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				parser.set(line);
				require(parser.terms()==4, "Contig observations require four fields");
				if(genome==null || !same(0, genome.idEncoded)){genome=inputs.get(decoded(0));}
				require(genome!=null && equalText(1, genome.fasta, genome.fastaEncoded), "Contig observation must bind to the exact roster input/FASTA");
				final String target=decoded(2);
				final long length=parser.parseLong(3);
				require(present(target) && length>=0 && length<=Integer.MAX_VALUE, "Observed contig identity/length is invalid");
				require(!genome.targets.contains(target), "Duplicate literal target within an input: "+target);
				genome.targets.put(target, new Target(target, (int)length));
			}
		}
		for(final Genome input:order){require(input.targets.size()>0, "Missing contig observations for input: "+input.id);}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Processing          ----------------*/
	/*--------------------------------------------------------------*/

	private void process(final Timer timer){
		if(rebase){
			rebase();
			timer.stop();
			outstream.println("Rebased declared reference evidence pointers. Time: "+timer);
			return;
		}
		loadRules();
		loadInputs();
		loadContigs();
		try{Files.createDirectory(outputPath);}
		catch(IOException failure){throw new IllegalStateException("Could not create a fresh reference directory: "+output, failure);}
		try(Output out=new Output(new File(output, "rules.tsv").getPath())){
			out.writer.println(RULE_HEADER);
			for(final byte[] row:ruleRows){out.writer.println(row);}
		}
		readDecisions();
		for(final Genome genome:order){
			Collections.sort(genome.identities, IDENTITY_ORDER);
			Identity previous=null;
			for(final Identity current:genome.identities){
				if(previous!=null && IDENTITY_ORDER.compare(previous, current)==0){
					throw new IllegalArgumentException("Duplicate or contradictory physical dispositions in input "+genome.id
						+" at source lines "+previous.line+" and "+current.line+"; resolve explicitly, never choose one");
				}
				if(current.row!=null){genome.emitted.add(current.row);}
				previous=current;
			}
		}
		for(final Genome genome:order){codec.write(genome.gff, genome.id, genome.emitted);}
		writeSummary();
		publishManifest();
		timer.stop();
		outstream.println("Serialized explicit dispositions for "+order.size()+" declared inputs. Time: "+timer);
		outstream.println("Recorded approval/authority are claims, not independently verified biological truth. No admission rules were inferred.");
	}

	/** Streams an existing manifest and each GFF once. Scientific fields and input
	 * identity remain unchanged; witness copying and biological approval are external.
	 * Publish only after all rewritten GFFs and the provenance summary close.
	 */
	private void rebase(){
		try{Files.createDirectory(outputPath);}
		catch(IOException failure){throw new IllegalStateException("Could not create fresh rebase directory: "+output, failure);}
		final File pending=new File(output, "reference.partial.tsv"), complete=new File(output, "reference.tsv");
		final ObjectSet<String> ids=new ObjectSet<String>(String.class), fastas=new ObjectSet<String>(String.class);
		long features=0;
		int count=0;
		try(Input in=new Input(decisionsPath); Output out=new Output(pending.getPath())){
			RrnaResourceIO.header(in.reader, GradeNcrnaReference.MANIFEST_HEADER, decisionsPath);
			out.writer.println(GradeNcrnaReference.MANIFEST_HEADER);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				parser.set(line);
				require(parser.terms()==3, "Reference manifest requires input_id, fasta and gff");
				final String id=literal(0), fasta=literal(1), source=literal(2);
				require(present(id) && rawSafe(id) && rawSafe(fasta) && new File(fasta).isAbsolute()
					&& new File(source).isAbsolute() && ids.add(id) && fastas.add(fasta),
					"Rebase manifest needs unique exact input IDs/FASTA paths and absolute GFF paths");
				final String destination=new File(output, "input_"+(++count)+".gff.gz").getPath();
				features=Math.addExact(features, codec.rebase(source, destination, id, evidenceFrom, evidenceTo));
				buffer.clear().append(id.getBytes(StandardCharsets.UTF_8)).tab().append(fasta.getBytes(StandardCharsets.UTF_8)).tab();
				out.writer.print(buffer.append(destination.getBytes(StandardCharsets.UTF_8)).nl());
			}
		}
		require(count>0, "An empty manifest cannot publish a rebased reference set");
		try(Output out=new Output(new File(output, "rebase.tsv").getPath())){
			out.writer.println("inputs\tfeature_rows\tsource_manifest\tevidence_from\tevidence_to\tqualification");
			buffer.clear().append(count).tab().append(features).tab();
			BuildNcrnaReference.field(buffer, decisionsPath).tab();
			BuildNcrnaReference.field(buffer, evidenceFrom).tab();
			BuildNcrnaReference.field(buffer, evidenceTo).append("\tDECLARED_PATH_CHANGE_NO_BIOLOGICAL_OR_WITNESS_COPY_VALIDATION\n");
			out.writer.print(buffer);
		}
		try{Files.move(pending.toPath(), complete.toPath(), StandardCopyOption.ATOMIC_MOVE);}
		catch(IOException failure){throw new IllegalStateException("Could not publish rebased reference manifest", failure);}
		outstream.println("Inputs rebased: "+count+"; feature evidence links changed: "+features);
	}

	private void readDecisions(){
		Genome genome=null;
		Rule rule=null;
		Target target=null;
		try(Input in=new Input(decisionsPath); Output out=new Output(companion)){
			RrnaResourceIO.header(in.reader, DECISION_HEADER, decisionsPath);
			out.writer.println(DECISION_HEADER);
			long number=1;
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				number++;
				parser.set(line);
				check(parser.terms()==13, "Disposition rows require thirteen fields", number);
				final Genome previous=genome;
				if(genome==null || !same(0, genome.idEncoded)){genome=inputs.get(decoded(0));}
				check(genome!=null, "Disposition input is absent from the exact roster", number);
				if(previous!=genome){target=null;}
				if(rule==null || !same(8, rule.idEncoded)){rule=rules.get(decoded(8));}
				check(rule!=null, "Missing explicit approved rule for disposition", number);
				check(equalText(5, rule.family, rule.familyEncoded) && parser.termEquals(STATES[rule.kind], 7),
					"Disposition family/status conflicts with its approved rule", number);
				final boolean resolved=parser.termEquals("RESOLVED", 6);
				check(resolved || parser.termEquals("UNRESOLVED", 6), "Family state must be explicit RESOLVED or UNRESOLVED", number);
				check(!resolved || present(rule.family), "A resolved family cannot be the absent-family marker", number);
				check(nonempty(11) && nonempty(12), "Every disposition needs decision and witness pointers", number);
				check(rule.kind==SCORED || nonempty(10), "A non-scored disposition needs an explicit reason", number);
				if(target==null || !same(1, target.encoded)){target=genome.targets.get(decoded(1));}
				final long start=coordinate(2), stop=coordinate(3);
				final char strand=parser.termEquals("+", 4) ? '+' : parser.termEquals("-", 4) ? '-' : '.';
				final boolean valid=target!=null && start>0 && stop>=start && stop<=target.length && strand!='.';
				final boolean emitting=rule.kind==SCORED || rule.kind==MASK;
				Row row=null;
				if(emitting){
					check(resolved && valid, "Emitted loci need resolved families and valid observed contig bounds/strand", number);
					final String reason=rule.kind==MASK ? decoded(10) : null;
					final String model=parser.termEquals(".", 9) ? null : decoded(9);
					row=codec.create(genome.id, new Locus(target.name, rule.family, (int)start, (int)stop, strand),
						rule.kind==MASK, reason, model, companion+":"+number);
				}
				if(resolved && valid){genome.identities.add(new Identity(target, rule.family, (int)start, (int)stop, strand, number, row));}
				genome.counts[rule.kind]++;
				genome.rows++;
				out.writer.println(line);
			}
		}
	}

	private void writeSummary(){
		try(Output out=new Output(new File(output, "summary.tsv").getPath())){
			out.writer.println("input_id\trows\tscored\tdont_care\tunresolved\tremoved\tcontext\tqualification");
			for(final Genome genome:order){
				long sum=0;
				for(final long count:genome.counts){sum+=count;}
				require(sum==genome.rows && genome.emitted.size()==genome.counts[SCORED]+genome.counts[MASK],
					"Every disposition must be conserved; only scored/masked rows enter GFFs");
				BuildNcrnaReference.field(buffer.clear(), genome.id).tab().append(genome.rows);
				for(final long count:genome.counts){buffer.tab().append(count);}
				out.writer.print(buffer.append("\tDECLARED_APPROVED_RULES_NOT_INDEPENDENTLY_VERIFIED\n"));
			}
		}
	}

	/** Same-directory atomic publication follows the existing NcrnaChildProcess
	 * ready-marker idiom. No final manifest exists before every payload closes.
	 */
	private void publishManifest(){
		final File pending=new File(output, "reference.partial.tsv"), complete=new File(output, "reference.tsv");
		try(Output out=new Output(pending.getPath())){
			out.writer.println(GradeNcrnaReference.MANIFEST_HEADER);
			for(final Genome genome:order){
				buffer.clear().append(genome.id.getBytes(StandardCharsets.UTF_8)).tab()
					.append(genome.fasta.getBytes(StandardCharsets.UTF_8)).tab().append(genome.gff.getBytes(StandardCharsets.UTF_8)).nl();
				out.writer.print(buffer);
			}
		}
		try{Files.move(pending.toPath(), complete.toPath(), StandardCopyOption.ATOMIC_MOVE);}
		catch(IOException failure){throw new IllegalStateException("Could not atomically publish completed reference manifest", failure);}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Parsing            ----------------*/
	/*--------------------------------------------------------------*/

	private String decoded(final int column){
		parser.setBounds(column);
		return decoder.decode(parser.line(), parser.a(), parser.b());
	}
	private String literal(final int column){
		parser.setBounds(column);
		final String text=new String(parser.line(), parser.a(), parser.b()-parser.a(), StandardCharsets.UTF_8);
		require(text.indexOf('\uFFFD')<0, "Invalid UTF-8 changes roster identity");
		return text;
	}
	private long coordinate(final int column){return parser.termEquals(".", column) ? Long.MIN_VALUE : parser.parseLong(column);}
	private boolean nonempty(final int column){
		parser.setBounds(column);
		if(parser.a()==parser.b() || parser.termEquals(".", column)){return false;}
		for(int i=parser.a(); i<parser.b(); i++){
			if(parser.line()[i]=='%' || parser.line()[i]<0){return present(decoded(column));}
		}
		return true;
	}
	private boolean same(final int column, final byte[] expected){
		parser.setBounds(column);
		if(parser.b()-parser.a()!=expected.length){return false;}
		for(int i=0; i<expected.length; i++){if(parser.line()[parser.a()+i]!=expected[i]){return false;}}
		return true;
	}
	private boolean equalText(final int column, final String expected, final byte[] encoded){
		return same(column, encoded) || decoded(column).equals(expected);
	}
	private static byte[] encoded(final String text){return BuildNcrnaReference.field(new ByteBuilder(), text).toBytes();}
	private static boolean present(final String value){return !value.isEmpty() && !value.equals(".");}
	private static boolean rawSafe(final String value){return value.indexOf('\t')<0 && value.indexOf('\r')<0 && value.indexOf('\n')<0;}
	private static int kind(final String status){
		for(int i=0; i<STATES.length; i++){if(STATES[i].equals(status)){return i;}}
		return -1;
	}
	private static void require(final boolean valid, final String reason){if(!valid){throw new IllegalArgumentException(reason);}}
	private static void check(final boolean valid, final String reason, final long line){
		if(!valid){throw new IllegalArgumentException(reason+" at disposition line "+line);}
	}

	/*--------------------------------------------------------------*/
	/*----------------             Types            ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Rule{
		Rule(final String id, final String family_, final int kind_){
			family=family_;
			kind=kind_;
			idEncoded=encoded(id);
			familyEncoded=encoded(family);
		}
		final String family;
		final int kind;
		final byte[] idEncoded, familyEncoded;
	}
	private static final class Target{
		Target(final String name_, final int length_){
			name=name_;
			length=length_;
			encoded=encoded(name);
		}
		final String name;
		final int length;
		final byte[] encoded;
	}
	private static final class Genome{
		Genome(final String id_, final String fasta_, final String gff_){
			id=id_;
			fasta=fasta_;
			gff=gff_;
			idEncoded=encoded(id);
			fastaEncoded=encoded(fasta);
		}
		final String id, fasta, gff;
		final byte[] idEncoded, fastaEncoded;
		final ObjectMap<String, Target> targets=new ObjectMap<String, Target>(16, String.class, Target.class);
		final ArrayList<Identity> identities=new ArrayList<Identity>();
		final ArrayList<Row> emitted=new ArrayList<Row>();
		final long[] counts=new long[STATES.length];
		long rows;
	}
	private static final class Identity{
		Identity(final Target target_, final String family_, final int start_, final int stop_, final char strand_, final long line_, final Row row_){
			target=target_;
			family=family_;
			start=start_;
			stop=stop_;
			strand=strand_;
			line=line_;
			row=row_;
		}
		final Target target;
		final String family;
		final int start, stop;
		final char strand;
		final long line;
		final Row row;
	}
	private static final Comparator<Identity> IDENTITY_ORDER=(a, b)->{
		int c=a.target.name.compareTo(b.target.name);
		if(c==0){c=a.family.compareTo(b.family);}
		if(c==0){c=Integer.compare(a.start, b.start);}
		if(c==0){c=Integer.compare(a.stop, b.stop);}
		return c==0 ? Character.compare(a.strand, b.strand) : c;
	};

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String decisionsPath, inputsPath, contigsPath, rulesPath, companion;
	private boolean rebase;
	private String evidenceFrom, evidenceTo;
	private final File output;
	private final Path outputPath;
	private final PrintStream outstream;
	private final LineParser1 parser=new LineParser1('\t');
	private final NcrnaHistoricalEvidence.Attributes decoder=new NcrnaHistoricalEvidence.Attributes();
	private final NcrnaReferenceGff codec=new NcrnaReferenceGff();
	private final ByteBuilder buffer=new ByteBuilder();
	private final ObjectMap<String, Rule> rules=new ObjectMap<String, Rule>(16, String.class, Rule.class);
	private final ObjectMap<String, Genome> inputs=new ObjectMap<String, Genome>(16, String.class, Genome.class);
	private final ArrayList<Genome> order=new ArrayList<Genome>();
	private final ArrayList<byte[]> ruleRows=new ArrayList<byte[]>();
	private static final int SCORED=0, MASK=1;
	private static final String[] STATES={"scored", "dont_care", "unresolved", "removed", "context"};
	public static final String INPUT_HEADER="input_id\tfasta";
	public static final String CONTIG_HEADER="input_id\tfasta\tliteral_target\tobserved_length";
	public static final String RULE_HEADER="rule_id\tfamily\tdisposition\tapproval\tauthority";
	public static final String DECISION_HEADER="input_id\tliteral_target\tstart\tstop\tstrand\treference_family\tfamily_state\tdisposition\trule_id\tmodel\treason\tdecision\twitnesses";
}
