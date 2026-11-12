package prok;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
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
import prok.NcrnaReferenceCalls.Call;
import prok.NcrnaReferenceCalls.Rule;
import prok.NcrnaReferenceMatcher.Result;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/** Compares explicit hypothetical dispositions against the same saved calls.
 * Derived from A_SampleByteFile. Native text owners and reusable byte parsers
 * perform metadata-only work; no caller, predicate evaluation or biological
 * arbitration runs. Outputs are diagnostic tables, never reference GFFs/manifests.
 *
 * An input-local decision_id declares correspondence between the two scenarios.
 * It is separate from physical coordinates and is not proof that the decisions
 * concern the same biological locus. Original provenance rows remain visible.
 * Scenario rows are retained; calls and matcher work are confined to one input.
 * @author Brian Bushnell, Raiden
 */
public final class CompareNcrnaScenarios{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final PreParser pp=new PreParser(args, CompareNcrnaScenarios.class, false);
		try{new CompareNcrnaScenarios(pp.args, pp.outstream).process(timer);}
		finally{Shared.closeStream(pp.outstream);}
	}

	private CompareNcrnaScenarios(final String[] args, final PrintStream stream){
		outstream=stream;
		Shared.setThreads(1);
		final Parser common=new Parser();
		common.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("calls")){callsPath=b;}
			else if(a.equals("contigs")){contigsPath=b;}
			else if(a.equals("familymap")){familyPath=b;}
			else if(a.equals("overlap")){overlap=Float.parseFloat(b);}
			else if(!common.parse(arg, a, b)){throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		beforePath=common.in1;
		afterPath=common.in2;
		require(beforePath!=null && afterPath!=null && callsPath!=null && contigsPath!=null && familyPath!=null && common.out1!=null,
			"Supply in=before in2=after calls=manifest contigs=observations familymap=projection out=new_directory");
		require(!common.overwrite && !common.append, "Scenario reports require a fresh directory");
		require(overlap>0 && overlap<=1, "Reciprocal scored overlap must be finite and in (0,1]");
		require(Tools.testInputFiles(false, true, beforePath, afterPath, callsPath, contigsPath, familyPath), "Scenario inputs must be readable");
		output=new File(common.out1).getAbsoluteFile().toPath();
		require(!Files.exists(output), "Scenario output directory must not exist: "+output);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Metadata           ----------------*/
	/*--------------------------------------------------------------*/

	private Rule[] loadFamilies(){
		final ArrayList<Rule> rules=new ArrayList<Rule>();
		try(Input in=new Input(familyPath)){
			RrnaResourceIO.header(in.reader, GradeNcrnaReference.RULE_HEADER, familyPath);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				parser.set(line);
				require(parser.terms()==4, "Family projection needs four fields");
				final boolean include=parser.termEquals("true", 3) || parser.termEquals("t", 3);
				require(include || parser.termEquals("false", 3) || parser.termEquals("f", 3), "Family include must be true/t or false/f");
				final Rule rule=new Rule(literal(0), literal(1), literal(2), include);
				rules.add(rule);
				family(rule.family).included|=include;
			}
		}
		return rules.toArray(new Rule[rules.size()]);
	}

	/** One common manifest supplies both scenarios' exact input/FASTA/call bindings. */
	private void loadInputs(){
		final ObjectSet<String> fastas=new ObjectSet<String>(String.class), files=new ObjectSet<String>(String.class);
		try(Input in=new Input(callsPath)){
			RrnaResourceIO.header(in.reader, GradeNcrnaReference.MANIFEST_HEADER, callsPath);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				parser.set(line);
				require(parser.terms()==3, "Call manifest needs input_id, fasta and gff");
				final String id=literal(0), fasta=literal(1), gff=literal(2);
				require(present(id) && new File(fasta).isAbsolute() && new File(gff).isAbsolute(), "Manifest needs an input identity and absolute paths");
				require(!inputs.contains(id) && fastas.add(fasta) && files.add(gff), "Duplicate input ID, FASTA path or call GFF");
				final Genome genome=new Genome(id, fasta, gff);
				inputs.put(id, genome);
				genomes.add(genome);
			}
		}
		require(!genomes.isEmpty(), "A common input population is required, including empty-call inputs");
	}

	private void loadContigs(){
		Genome genome=null;
		try(Input in=new Input(contigsPath)){
			RrnaResourceIO.header(in.reader, WriteNcrnaReference.CONTIG_HEADER, contigsPath);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				parser.set(line);
				require(parser.terms()==4, "Contig observations need four fields");
				if(genome==null || !same(0, genome.encoded)){genome=inputs.get(decoded(0));}
				require(genome!=null && (same(1, genome.fastaEncoded) || decoded(1).equals(genome.fasta)), "Contigs must bind to the exact input/FASTA");
				final String name=decoded(2);
				final long length=parser.parseLong(3);
				require(present(name) && length>=0 && length<=Integer.MAX_VALUE && !genome.targets.contains(name), "Invalid or repeated observed contig");
				genome.targets.put(name, new Target(name, (int)length));
			}
		}
		for(final Genome input:genomes){require(input.targets.size()>0, "Missing contig observations for "+input.id);}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Scenarios           ----------------*/
	/*--------------------------------------------------------------*/

	private void readScenario(final String path, final int side, final Output retained){
		Genome genome=null;
		Family family=null;
		Target target=null;
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, HEADER, path);
			long number=1;
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				number++;
				parser.set(line);
				require(parser.terms()==14, "Scenario rows require decision_id followed by the thirteen disposition fields");
				final Genome previous=genome;
				if(genome==null || !same(1, genome.encoded)){genome=inputs.get(decoded(1));}
				require(genome!=null, "Scenario input is absent from the common call manifest");
				if(genome!=previous){target=null;}
				final String key=decoded(0);
				require(present(key), "A stable input-local decision_id is required; correspondence is never guessed");
				Entry entry=genome.entries.get(key);
				if(entry==null){
					entry=new Entry(key);
					genome.entries.put(key, entry);
					genome.order.add(entry);
				}
				require(entry.rows[side]==null, "Duplicate decision_id in one scenario/input: "+key);
				if(family==null || !same(6, family.encoded)){family=family(decoded(6));}
				final int state=state(8);
				final boolean resolved=parser.termEquals("RESOLVED", 7);
				require(resolved || parser.termEquals("UNRESOLVED", 7), "Family state must be explicitly RESOLVED or UNRESOLVED");
				require(!resolved || present(family.name), "A resolved family must be named");
				require(nonempty(9) && nonempty(12) && nonempty(13), "Scenario rule, decision and witness pointers are required, without approval claims");
				require(state==SCORED || nonempty(11), "Non-scored dispositions need an explicit reason");
				final String targetName;
				if(target!=null && same(2, target.encoded)){targetName=target.name;}
				else{
					targetName=decoded(2);
					target=genome.targets.get(targetName);
				}
				final long start=coordinate(3), stop=coordinate(4);
				final char strand=parser.termEquals("+", 5) ? '+' : parser.termEquals("-", 5) ? '-' : '.';
				final boolean valid=target!=null && start>0 && stop>=start && stop<=target.length && strand!='.';
				final boolean scored=state==SCORED || state==MASK;
				require(!scored || (resolved && valid && family.included), "Scored/masked scenario rows need resolved included families and valid observed bounds/strand");
				final Model model=family.model(parser, decoder, 10);
				final Snapshot snapshot=new Snapshot(entry, family, model, targetName, start, stop, strand, state, resolved, number);
				if(scored){
					snapshot.row=codec.create(genome.id, new Locus(target.name, family.name, (int)start, (int)stop, strand),
						state==MASK, state==MASK ? decoded(11) : null, model.name.equals(".") ? null : model.name, path+":"+number);
					genome.scored[side].add(snapshot);
				}
				if(resolved && valid){genome.physical[side].add(snapshot);}
				entry.rows[side]=snapshot;
				family.counts[side][state]++;
				model.counts[side][state]++;
				rowCounts[side]++;
				buffer.clear().append(SIDES[side]).tab().append(number).tab().append(line).nl();
				retained.writer.print(buffer);
			}
		}
		for(final Genome item:genomes){
			Collections.sort(item.physical[side], PHYSICAL);
			Collections.sort(item.scored[side], PHYSICAL);
			Snapshot previous=null;
			for(final Snapshot current:item.physical[side]){
				require(previous==null || PHYSICAL.compare(previous, current)!=0, "Duplicate/contradictory physical dispositions within a scenario/input");
				previous=current;
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Comparison          ----------------*/
	/*--------------------------------------------------------------*/

	private void process(final Timer timer){
		final NcrnaReferenceCalls calls=new NcrnaReferenceCalls(loadFamilies());
		loadInputs();
		loadContigs();
		try{Files.createDirectory(output);}
		catch(IOException failure){throw new IllegalStateException("Could not create diagnostic directory", failure);}
		writeSources();
		try(Output retained=new Output(path("dispositions.tsv.gz"))){
			retained.writer.println("scenario\tsource_line\t"+HEADER);
			readScenario(beforePath, 0, retained);
			readScenario(afterPath, 1, retained);
		}
		try(Output changes=new Output(path("changes.tsv")); Output outcomes=new Output(path("calls.tsv.gz")); Output inventory=new Output(path("inputs.tsv"))){
			changes.writer.println(CHANGE_HEADER);
			outcomes.writer.println(CALL_HEADER);
			inventory.writer.println("input_id\tfasta\tquery_gff\tfeatures\tCDS\ttRNA\trRNA\tincluded_calls\texcluded_calls\tbefore_decisions\tafter_decisions\tqualification");
			for(final Genome genome:genomes){
				writeChanges(genome, changes);
				compareCalls(genome, calls.read(genome.gff, genome.id), outcomes, inventory);
			}
		}
		writeCounts();
		writeSummary();
		timer.stop();
		outstream.println("Compared two explicit hypothetical references for "+genomes.size()+" inputs. Time: "+timer);
		outstream.println(QUALIFICATION+": metric changes follow changed reference assumptions, not measured biological improvement.");
	}

	private void writeChanges(final Genome genome, final Output out){
		Collections.sort(genome.order, (a, b)->a.key.compareTo(b.key));
		for(final Entry entry:genome.order){
			final Snapshot a=entry.rows[0], b=entry.rows[1];
			final boolean paired=a!=null && b!=null;
			final boolean coordinate=paired && (a.start!=b.start || a.stop!=b.stop || a.strand!=b.strand || !a.target.equals(b.target));
			final boolean family=paired && (!a.family.name.equals(b.family.name) || a.resolved!=b.resolved);
			final boolean model=paired && !a.model.name.equals(b.model.name);
			final boolean disposition=paired && a.state!=b.state;
			buffer.clear();
			field(genome.id).tab();
			field(entry.key).tab().append(a==null ? "ADDED" : b==null ? "REMOVED" : "PAIRED").tab();
			for(final Snapshot row:entry.rows){
				if(row==null){buffer.append(".\tABSENT\t.\t.\t.\t.\t.\t.\t.\t");}
				else{
					buffer.append(row.line).tab().append(STATES[row.state]).tab();
					field(row.family.name).tab();
					field(row.model.name).tab();
					field(row.target).tab();
					coordinate(row.start).tab();
					coordinate(row.stop).tab().append(row.strand).tab().append(row.resolved).tab();
				}
			}
			buffer.append(coordinate).tab().append(family).tab().append(model).tab().append(disposition).tab();
			out.writer.print(buffer.append(QUALIFICATION).nl());
		}
	}

	private void compareCalls(final Genome genome, final NcrnaReferenceCalls.Result calls, final Output out, final Output inventory){
		final ArrayList<Locus> selected=new ArrayList<Locus>();
		for(final Call call:calls.calls){
			final Target target=genome.targets.get(call.locus.seqid);
			require(target!=null && call.locus.stop<=target.length, "Saved query is outside the declared observed contig");
			if(call.included){selected.add(call.locus);}
		}
		final Result[] matches=new Result[2];
		for(int side=0; side<2; side++){
			final ArrayList<Snapshot> refs=genome.scored[side];
			final ArrayList<Locus> loci=new ArrayList<Locus>(refs.size());
			final boolean[] masks=new boolean[refs.size()];
			for(int i=0; i<refs.size(); i++){
				final Snapshot ref=refs.get(i);
				loci.add(ref.row.locus);
				masks[i]=ref.state==MASK;
			}
			final Result result=NcrnaReferenceMatcher.match(selected, ones(selected.size()), loci, ones(loci.size()), masks, overlap);
			matches[side]=result;
			for(int i=0; i<refs.size(); i++){
				final Snapshot ref=refs.get(i);
				final long[] counts=ref.family.metrics[side];
				if(ref.state==MASK){counts[MASK_REFS]++;}
				else{
					counts[SCORED_REFS]++;
					if(result.referenceToCall[i]<0){counts[FN]++;}
				}
			}
		}
		int index=0;
		for(final Call call:calls.calls){
			buffer.clear();
			field(genome.id).tab().append(call.line).tab();
			field(call.locus.family).tab();
			field(call.model).tab();
			final int ci=call.included ? index++ : -1;
			Snapshot before=null, after=null;
			int first=-1, last=-1;
			for(int side=0; side<2; side++){
				final Result result=matches[side];
				final int pair=ci<0 ? -1 : result.callToReference[ci], mask=ci<0 ? -1 : result.callMaskReference[ci];
				final int outcome=ci<0 ? EXCLUDED : pair>=0 ? TP : mask>=0 ? MASK_CALLS : FP;
				final Snapshot witness=pair>=0 ? genome.scored[side].get(pair) : mask>=0 ? genome.scored[side].get(mask) : null;
				familyByName.get(call.locus.family).metrics[side][outcome]++;
				if(ci>=0 && result.duplicateCalls[ci]){familyByName.get(call.locus.family).metrics[side][DUPLICATES]++;}
				buffer.append(outcome==TP ? "TP" : outcome==FP ? "FP" : outcome==MASK_CALLS ? "MASKED" : "EXCLUDED").tab();
				field(witness==null ? "." : witness.entry.key).tab();
				field(witness==null ? "." : witness.row.id).tab();
				if(witness==null){buffer.append('.');}else{buffer.append(witness.line);}
				buffer.tab().append(ci>=0 && result.duplicateCalls[ci]).tab().append(ci>=0 && result.scoredEligible[ci]).tab();
				if(side==0){
					before=witness;
					first=outcome;
				}else{
					after=witness;
					last=outcome;
				}
			}
			final boolean changed=first!=last || !sameWitness(before, after)
				|| (ci>=0 && matches[0].scoredEligible[ci]!=matches[1].scoredEligible[ci]);
			out.writer.print(buffer.append(changed).tab().append(QUALIFICATION).nl());
		}
		assert(index==selected.size()) : "Included call indices must match both matcher result arrays";
		buffer.clear();
		field(genome.id).tab();
		field(genome.fasta).tab();
		field(genome.gff).tab().append(calls.features).tab().append(calls.cds).tab().append(calls.trna).tab().append(calls.rrna)
			.tab().append(calls.included).tab().append(calls.excluded);
		for(int side=0; side<2; side++){
			long count=0;
			for(final Entry entry:genome.order){if(entry.rows[side]!=null){count++;}}
			buffer.tab().append(count);
		}
		inventory.writer.print(buffer.tab().append(QUALIFICATION).nl());
	}

	/*--------------------------------------------------------------*/
	/*----------------            Reports           ----------------*/
	/*--------------------------------------------------------------*/

	/** Small provenance map for the preserved source-line witnesses; not a reference manifest. */
	private void writeSources(){
		try(Output out=new Output(path("sources.tsv"))){
			out.writer.println("role\tpath\tqualification");
			final String[] roles={"before", "after", "calls", "contigs", "familymap"};
			final String[] paths={beforePath, afterPath, callsPath, contigsPath, familyPath};
			for(int i=0; i<paths.length; i++){
				buffer.clear().append(roles[i]).tab();
				out.writer.print(field(new File(paths[i]).getAbsolutePath()).tab().append(QUALIFICATION).nl());
			}
		}
	}

	private void writeCounts(){
		final long[] totals=new long[2];
		try(Output out=new Output(path("counts.tsv"))){
			out.writer.println("scope\tfamily\tmodel\tscenario\tscored\tdont_care\tunresolved\tremoved\tcontext\tqualification");
			for(final Family family:families){
				writeCount(out, "FAMILY", family.name, ".", family.counts);
				for(final Model model:family.models){writeCount(out, "MODEL", family.name, model.name, model.counts);}
				for(int side=0; side<2; side++){for(final long count:family.counts[side]){totals[side]+=count;}}
			}
		}
		require(Arrays.equals(totals, rowCounts), "Disposition accounting must preserve every scenario row including context and unresolved evidence");
	}
	private void writeCount(final Output out, final String scope, final String family, final String model, final long[][] counts){
		for(int side=0; side<2; side++){
			buffer.clear().append(scope).tab();
			field(family).tab();
			field(model).tab().append(SIDES[side]);
			for(final long count:counts[side]){buffer.tab().append(count);}
			out.writer.print(buffer.tab().append(QUALIFICATION).nl());
		}
	}
	/** Written last; it is a diagnostic completion table, never a reference manifest. */
	private void writeSummary(){
		try(Output out=new Output(path("summary.tsv"))){
			out.writer.println(SUMMARY_HEADER);
			for(final Family family:families){
				if(!family.included && !family.hasQueries()){continue;}
				buffer.clear();
				field(family.name).tab().append(genomes.size());
				for(int side=0; side<2; side++){
					final long[] counts=family.metrics[side];
					require(counts[TP]+counts[FN]==counts[SCORED_REFS], "Hypothetical TP/FN must conserve that scenario's scored denominator");
					for(final long count:counts){buffer.tab().append(count);}
				}
				buffer.tab().append(family.metrics[1][TP]-family.metrics[0][TP]).tab().append(family.metrics[1][FP]-family.metrics[0][FP])
					.tab().append(family.metrics[1][FN]-family.metrics[0][FN]).tab().append(family.metrics[1][MASK_CALLS]-family.metrics[0][MASK_CALLS]);
				out.writer.print(buffer.tab().append(overlap, 6).tab().append(QUALIFICATION).nl());
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Helpers           ----------------*/
	/*--------------------------------------------------------------*/

	private Family family(final String name){
		Family family=familyByName.get(name);
		if(family==null){
			require(!name.isEmpty(), "Scenario family must be named or explicitly unknown as a dot");
			family=new Family(name);
			familyByName.put(name, family);
			families.add(family);
		}
		return family;
	}
	private String path(final String file){return output.resolve(file).toString();}
	private String decoded(final int column){
		parser.setBounds(column);
		return decoder.decode(parser.line(), parser.a(), parser.b());
	}
	private String literal(final int column){
		parser.setBounds(column);
		final String value=new String(parser.line(), parser.a(), parser.b()-parser.a(), StandardCharsets.UTF_8);
		require(value.indexOf('\uFFFD')<0, "Invalid UTF-8 changes exact input identity");
		return value;
	}
	private boolean nonempty(final int column){
		parser.setBounds(column);
		if(parser.a()==parser.b() || parser.termEquals(".", column)){return false;}
		for(int i=parser.a(); i<parser.b(); i++){if(parser.line()[i]=='%' || parser.line()[i]<0){return present(decoded(column));}}
		return true;
	}
	private int state(final int column){
		for(int i=0; i<STATES.length; i++){if(parser.termEquals(STATES[i], column)){return i;}}
		throw new IllegalArgumentException("Unknown explicit scenario disposition");
	}
	private long coordinate(final int column){return parser.termEquals(".", column) ? Long.MIN_VALUE : parser.parseLong(column);}
	private ByteBuilder coordinate(final long value){return value==Long.MIN_VALUE ? buffer.append('.') : buffer.append(value);}
	private ByteBuilder field(final String value){return BuildNcrnaReference.field(buffer, value);}
	private boolean same(final int column, final byte[] encoded){return same(parser, column, encoded);}
	private static boolean same(final LineParser1 parser, final int column, final byte[] encoded){
		parser.setBounds(column);
		if(parser.b()-parser.a()!=encoded.length){return false;}
		for(int i=0; i<encoded.length; i++){if(parser.line()[parser.a()+i]!=encoded[i]){return false;}}
		return true;
	}
	private static byte[] encoded(final String text){return BuildNcrnaReference.field(new ByteBuilder(), text).toBytes();}
	private static boolean present(final String value){return !value.isEmpty() && !value.equals(".");}
	private static boolean sameWitness(final Snapshot a, final Snapshot b){
		return a==null ? b==null : b!=null && a.entry.key.equals(b.entry.key) && a.row.id.equals(b.row.id);
	}
	private static int[] ones(final int size){
		final int[] result=new int[size];
		Arrays.fill(result, 1);
		return result;
	}
	private static void require(final boolean valid, final String reason){if(!valid){throw new IllegalArgumentException(reason);}}

	/*--------------------------------------------------------------*/
	/*----------------             Types            ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Genome{
		@SuppressWarnings("unchecked")
		Genome(final String id_, final String fasta_, final String gff_){
			id=id_;
			fasta=fasta_;
			gff=gff_;
			encoded=encoded(id);
			fastaEncoded=encoded(fasta);
			for(int i=0; i<2; i++){
				scored[i]=new ArrayList<Snapshot>();
				physical[i]=new ArrayList<Snapshot>();
			}
		}
		final String id, fasta, gff;
		final byte[] encoded, fastaEncoded;
		final ObjectMap<String,Target> targets=new ObjectMap<String,Target>(16, String.class, Target.class);
		final ObjectMap<String,Entry> entries=new ObjectMap<String,Entry>(16, String.class, Entry.class);
		final ArrayList<Entry> order=new ArrayList<Entry>();
		@SuppressWarnings("unchecked")
		final ArrayList<Snapshot>[] scored=new ArrayList[2], physical=new ArrayList[2];
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
	private static final class Entry{
		Entry(final String key_){key=key_;}
		final String key;
		final Snapshot[] rows=new Snapshot[2];
	}
	private static final class Family{
		Family(final String name_){
			name=name_;
			encoded=encoded(name);
		}
		Model model(final LineParser1 parser, final NcrnaHistoricalEvidence.Attributes decoder, final int column){
			if(lastModel!=null && same(parser, column, lastModel.encoded)){return lastModel;}
			parser.setBounds(column);
			final String name=decoder.decode(parser.line(), parser.a(), parser.b());
			require(!name.isEmpty(), "Model must be named or explicitly absent as a dot");
			lastModel=byModel.get(name);
			if(lastModel==null){
				lastModel=new Model(name);
				byModel.put(name, lastModel);
				models.add(lastModel);
			}
			return lastModel;
		}
		boolean hasQueries(){return metrics[0][EXCLUDED]>0 || metrics[1][EXCLUDED]>0;}
		final String name;
		final byte[] encoded;
		boolean included;
		final long[][] counts=new long[2][5], metrics=new long[2][8];
		final ObjectMap<String,Model> byModel=new ObjectMap<String,Model>(4, String.class, Model.class);
		final ArrayList<Model> models=new ArrayList<Model>();
		Model lastModel;
	}
	private static final class Model{
		Model(final String name_){
			name=name_;
			encoded=encoded(name);
		}
		final String name;
		final byte[] encoded;
		final long[][] counts=new long[2][5];
	}
	private static final class Snapshot{
		Snapshot(final Entry entry_, final Family family_, final Model model_, final String target_, final long start_, final long stop_,
			final char strand_, final int state_, final boolean resolved_, final long line_){
			entry=entry_;
			family=family_;
			model=model_;
			target=target_;
			start=start_;
			stop=stop_;
			strand=strand_;
			state=state_;
			resolved=resolved_;
			line=line_;
		}
		final Entry entry;
		final Family family;
		final Model model;
		final String target;
		final long start, stop, line;
		final char strand;
		final int state;
		final boolean resolved;
		NcrnaReferenceGff.Row row;
	}
	private static final Comparator<Snapshot> PHYSICAL=(a, b)->{
		int c=a.target.compareTo(b.target);
		if(c==0){c=a.family.name.compareTo(b.family.name);}
		if(c==0){c=Long.compare(a.start, b.start);}
		if(c==0){c=Long.compare(a.stop, b.stop);}
		return c==0 ? Character.compare(a.strand, b.strand) : c;
	};

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String beforePath, afterPath, callsPath, contigsPath, familyPath;
	private float overlap=.5f;
	private final Path output;
	private final PrintStream outstream;
	private final LineParser1 parser=new LineParser1('\t');
	private final NcrnaHistoricalEvidence.Attributes decoder=new NcrnaHistoricalEvidence.Attributes();
	private final NcrnaReferenceGff codec=new NcrnaReferenceGff();
	private final ByteBuilder buffer=new ByteBuilder();
	private final ObjectMap<String,Genome> inputs=new ObjectMap<String,Genome>(16, String.class, Genome.class);
	private final ArrayList<Genome> genomes=new ArrayList<Genome>();
	private final ObjectMap<String,Family> familyByName=new ObjectMap<String,Family>(16, String.class, Family.class);
	private final ArrayList<Family> families=new ArrayList<Family>();
	private final long[] rowCounts=new long[2];
	private static final int SCORED=0, MASK=1;
	private static final int TP=0, FP=1, FN=2, SCORED_REFS=3, MASK_REFS=4, MASK_CALLS=5, DUPLICATES=6, EXCLUDED=7;
	private static final String[] SIDES={"before", "after"};
	private static final String[] STATES={"scored", "dont_care", "unresolved", "removed", "context"};
	public static final String HEADER="decision_id\t"+WriteNcrnaReference.DECISION_HEADER;
	public static final String QUALIFICATION="HYPOTHETICAL_REFERENCE_COMPARISON_NOT_BIOLOGICAL_ACCURACY";
	public static final String CHANGE_HEADER="input_id\tdecision_id\trelation\tbefore_line\tbefore_state\tbefore_family\tbefore_model\tbefore_target\tbefore_start\tbefore_stop\tbefore_strand\tbefore_resolved\tafter_line\tafter_state\tafter_family\tafter_model\tafter_target\tafter_start\tafter_stop\tafter_strand\tafter_resolved\tcoordinates_changed\tfamily_changed\tmodel_changed\tdisposition_changed\tqualification";
	public static final String CALL_HEADER="input_id\tquery_line\tfamily\tquery_model\tbefore_outcome\tbefore_decision\tbefore_physical_id\tbefore_source_line\tbefore_duplicate\tbefore_scored_eligible\tafter_outcome\tafter_decision\tafter_physical_id\tafter_source_line\tafter_duplicate\tafter_scored_eligible\tchanged\tqualification";
	public static final String SUMMARY_HEADER="family\tinputs\tbefore_TP\tbefore_FP\tbefore_FN\tbefore_scored\tbefore_masks\tbefore_masked_calls\tbefore_duplicates\tbefore_excluded\tafter_TP\tafter_FP\tafter_FN\tafter_scored\tafter_masks\tafter_masked_calls\tafter_duplicates\tafter_excluded\tdelta_TP\tdelta_FP\tdelta_FN\tdelta_masked_calls\toverlap\tqualification";
}
