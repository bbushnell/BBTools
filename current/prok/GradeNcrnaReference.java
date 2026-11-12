package prok;

import java.io.File;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import map.ObjectMap;
import map.ObjectIntMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaReferenceCalls.Call;
import prok.NcrnaReferenceCalls.Rule;
import prok.NcrnaReferenceGff.Row;
import prok.NcrnaReferenceMatcher.Result;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntList;

/** Grades per-input CallGenes queries against explicitly scored/masked references.
 * Adapted from template.A_SampleByteFile: PreParser/Parser, native byte streams,
 * reusable output buffers and checked closure. Computation processes one input
 * at a time; no sequence payloads, calling, admission or threshold selection.
 *
 * By default query/reference manifests must have equal exact input IDs and FASTA
 * paths. Explicit restrictreference=t selects the complete query input population
 * and included family projection from a larger reference. Excluded input counts
 * are separate from reference-row counts in inspected on-panel files.
 * These are checked declarations, not independent proof of producer provenance.
 * A combined saved GFF needs its original input mapping before use here.
 * Reference GFF status controls masking, with scored eligibility taking precedence.
 * Summary is written only after all detailed outputs close successfully.
 * @author Brian Bushnell, Raiden
 */
public final class GradeNcrnaReference{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final GradeNcrnaReference tool=new GradeNcrnaReference(args);
		try{tool.process(timer);}
		finally{Shared.closeStream(tool.outstream);}
	}

	private GradeNcrnaReference(final String[] arguments){
		final PreParser pp=new PreParser(arguments, GradeNcrnaReference.class, false);
		outstream=pp.outstream;
		Shared.setThreads(1);
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:pp.args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("reference") || a.equals("ref")){
				reference=b;
			}else if(a.equals("familymap")){
				familymap=b;
			}else if(a.equals("overlap")){
				overlap=Float.parseFloat(b);
			}else if(a.equals("restrictreference")){
				restrictReference=Parse.parseBoolean(b);
			}else if(a.equals("subtypes")){
				subtypes=b;
			}else if(a.equals("witnesses")){
				witnesses=b;
			}else if(a.equals("events")){
				events=b;
			}else if(!parser.parse(arg, a, b)){
				throw new IllegalArgumentException("Unknown parameter: "+arg);
			}
		}
		input=parser.in1;
		prefix=parser.out1;
		require(input!=null && reference!=null && familymap!=null && prefix!=null,
			"Supply in=query_manifest reference=reference_manifest familymap=rules out=prefix");
		require(!parser.append && !parser.overwrite, "Grading outputs must be fresh; append/overwrite are unsupported");
		require(overlap>0 && overlap<=1, "Scored reciprocal overlap must be finite and in (0,1]");
		require(Tools.testInputFiles(false, true, input, reference, familymap), "Grading manifests and family map must be readable");
		require(subtypes==null ? witnesses==null && events==null : witnesses!=null && events!=null,
			"Subtype reporting requires subtypes= mapping, witnesses= and events= together");
		if(subtypes!=null){
			require(Tools.testInputFiles(false, true, subtypes, witnesses, events), "Subtype mapping and evidence sidecars must be readable");
		}
		final int baseOutputs=restrictReference ? 4 : 3;
		outputPaths=new String[baseOutputs+(subtypes==null ? 0 : NcrnaSubtypeGrade.SUFFIXES.length)];
		for(int i=0; i<baseOutputs; i++){outputPaths[i]=new File(prefix+SUFFIXES[i]).getAbsolutePath();}
		for(int i=baseOutputs; i<outputPaths.length; i++){
			outputPaths[i]=new File(prefix+NcrnaSubtypeGrade.SUFFIXES[i-baseOutputs]).getAbsolutePath();
		}
		require(Tools.testOutputFiles(false, false, false, outputPaths),
			"Grading outputs must be fresh and writable");
		for(final String path:new String[]{input, reference, familymap, subtypes, witnesses, events}){
			if(path!=null){checkOutputCollision(new File(path).getAbsolutePath());}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Processing          ----------------*/
	/*--------------------------------------------------------------*/

	private void process(final Timer timer){
		final Rule[] rules=loadRules(familymap);
		final NcrnaReferenceCalls reader=new NcrnaReferenceCalls(rules);
		final ArrayList<Binding> referencePopulation=loadManifest(reference);
		final ArrayList<Binding> queries=loadManifest(input);
		for(final Binding binding:referencePopulation){checkBindingOutputs(binding);}
		for(final Binding binding:queries){checkBindingOutputs(binding);}
		require(restrictReference || referencePopulation.size()==queries.size(), "Query and reference manifests must have equal input populations, including zero-call files");
		final ArrayList<Binding> references=new ArrayList<Binding>(queries.size());
		final ObjectMap<String, Binding> queryById=new ObjectMap<String, Binding>(16, String.class, Binding.class);
		for(final Binding query:queries){queryById.put(query.id, query);}
		for(final Binding ref:referencePopulation){
			final Binding query=queryById.get(ref.id);
			if(query==null && restrictReference){continue;}
			require(query!=null, "Missing query input: "+ref.id);
			require(ref.fasta.equals(query.fasta), "Original FASTA paths disagree for exact input: "+ref.id);
			query.matched=true;
			references.add(ref);
		}
		for(final Binding query:queries){require(query.matched, "Missing reference input: "+query.id);}
		final NcrnaSubtypeGrade subtypeGrade;
		if(subtypes==null){subtypeGrade=null;}
		else{
			final ArrayList<String> inputIds=new ArrayList<String>(queries.size());
			for(final Binding query:queries){inputIds.add(query.id);}
			subtypeGrade=new NcrnaSubtypeGrade(subtypes, witnesses, events, rules, inputIds);
		}
		for(final Rule rule:rules){
			Stats stats=byFamily.get(rule.family);
			if(stats==null){
				stats=new Stats(rule.family);
				byFamily.put(rule.family, stats);
				families.add(stats);
			}
			stats.included|=rule.include;
		}
		try(Output raw=new Output(prefix+".raw.tsv"); Output inputs=new Output(prefix+".inputs.tsv");
			Output selection=restrictReference ? new Output(prefix+".selection.tsv") : null;
			Output subtypeRaw=subtypeGrade==null ? null : new Output(prefix+NcrnaSubtypeGrade.SUFFIXES[0]);
			Output subtypeSupport=subtypeGrade==null ? null : new Output(prefix+NcrnaSubtypeGrade.SUFFIXES[1])){
			raw.writer.println(RAW_HEADER);
			inputs.writer.println(INPUT_HEADER);
			if(subtypeGrade!=null){
				subtypeRaw.writer.println(NcrnaSubtypeGrade.RAW_HEADER);
				subtypeSupport.writer.println(NcrnaSubtypeGrade.SUPPORT_HEADER);
			}
			if(selection!=null){
				selection.writer.println(SELECTION_HEADER);
				buffer.clear().append("INPUT_POPULATION\t.\t.\t").append(references.size()).tab()
					.append(referencePopulation.size()-references.size()).append("\t.\t.\t.\t.\tEXACT_QUERY_INPUT_SET_OFF_PANEL_GFFS_NOT_READ\n");
				flush(selection);
			}
			for(final Binding ref:references){
				grade(ref, queryById.get(ref.id), reader, raw, inputs, selection, subtypeGrade, subtypeRaw, subtypeSupport);
			}
		}
		if(subtypeGrade!=null){
			subtypeGrade.writeSummary(prefix, references.size(), overlap);
			for(final String qualification:NcrnaSubtypeGrade.QUALIFICATIONS){outstream.println(qualification);}
		}
		writeSummary(references.size());
		timer.stop();
		outstream.println("Graded "+references.size()+" declared inputs against explicit reference dispositions. Time: "+timer);
		outstream.println("Manifest agreement is not independent producer provenance; no biological admission or confidence calibration.");
	}

	/** Retains only one input's loci while accumulating primitive per-family counts. */
	private void grade(final Binding ref, final Binding query, final NcrnaReferenceCalls reader,
		final Output raw, final Output inputOutput, final Output selection, final NcrnaSubtypeGrade subtypeGrade,
		final Output subtypeRaw, final Output subtypeSupport){
		final ArrayList<Row> rows=codec.read(ref.gff, ref.id);
		if(restrictReference){restrict(rows, ref.id, selection);}
		final NcrnaReferenceCalls.Result queries=reader.read(query.gff, query.id);
		final ArrayList<Locus> truth=new ArrayList<Locus>(rows.size());
		final boolean[] masks=new boolean[rows.size()];
		for(int i=0; i<rows.size(); i++){
			final Row row=rows.get(i);
			final Stats stats=byFamily.get(row.locus.family);
			if(stats==null || !stats.included){throw new IllegalArgumentException("Reference family lacks an included query projection: "+row.locus.family);}
			truth.add(row.locus);
			masks[i]=row.mask;
			if(row.mask){stats.maskReferences++;}
			else{stats.scoredReferences++;}
		}
		final ArrayList<Locus> selected=new ArrayList<Locus>();
		final IntList original=new IntList();
		for(int i=0; i<queries.calls.size(); i++){
			final Call call=queries.calls.get(i);
			if(call.included){
				selected.add(call.locus);
				original.add(i);
			}
		}
		final Result match=NcrnaReferenceMatcher.match(selected, ones(selected.size()), truth, ones(truth.size()), masks, overlap);
		int selectedIndex=0;
		for(final Call call:queries.calls){
			final Stats stats=byFamily.get(call.locus.family);
			if(!call.included){
				stats.excluded++;
				writeCall(raw, ref.id, call, -1, rows, match, stats);
			}else{
				writeCall(raw, ref.id, call, selectedIndex++, rows, match, stats);
			}
		}
		assert(selectedIndex==selected.size()) : "Every included query must consume exactly one original-index matcher result";
		for(int i=0; i<rows.size(); i++){
			final Row row=rows.get(i);
			final int matched=match.referenceToCall[i];
			final String outcome=row.mask ? "DONT_CARE" : matched>=0 ? "TP" : "FN";
			if(!row.mask && matched<0){byFamily.get(row.locus.family).fn++;}
			startLine(ref.id, "REFERENCE", row.locus, outcome);
			field(row.id).append("\t.\t.\t");
			if(matched<0){buffer.append('.');}
			else{buffer.append(queries.calls.get(original.get(matched)).line);}
			buffer.append("\t.\t.\t.\t");
			field(row.model==null ? "." : row.model).append("\t.\t.\n");
			flush(raw);
		}
		field(ref.id).tab();
		field(ref.fasta).tab();
		field(query.gff).tab();
		field(ref.gff).tab().append(queries.features).tab().append(queries.cds).tab().append(queries.trna)
			.tab().append(queries.rrna).tab().append(queries.included).tab().append(queries.excluded).nl();
		flush(inputOutput);
		if(subtypeGrade!=null){subtypeGrade.grade(ref.id, rows, queries.calls, overlap, subtypeRaw, subtypeSupport);}
	}

	/** Compacts the one-input list in place. Family enablement comes only from the
	 * declared projection, never observed calls; zero-call inputs therefore keep FNs.
	 * Counters are primitive, reused across inputs, and count both retained statuses.
	 */
	private void restrict(final ArrayList<Row> rows, final String inputId, final Output out){
		assert(out!=null) : "Explicit reference restriction always publishes its separate selection accounting";
		selectionIndex.clear();
		selectionFamilies.clear();
		for(final IntList counts:selectionCounts){counts.clear();}
		final int original=rows.size();
		int kept=0;
		for(int i=0; i<original; i++){
			final Row row=rows.get(i);
			final String family=row.locus.family;
			int index=selectionIndex.get(family);
			if(index<0){
				index=selectionFamilies.size();
				selectionIndex.put(family, index);
				selectionFamilies.add(family);
				for(final IntList counts:selectionCounts){counts.add(0);}
			}
			final Stats stats=byFamily.get(family);
			final boolean keep=stats!=null && stats.included;
			final IntList counts=selectionCounts[(keep ? 0 : 2)+(row.mask ? 1 : 0)];
			counts.set(index, counts.get(index)+1);
			if(keep){rows.set(kept++, row);}
		}
		while(rows.size()>kept){rows.remove(rows.size()-1);}
		long accounted=0;
		for(int i=0; i<selectionFamilies.size(); i++){
			final String family=selectionFamilies.get(i);
			final Stats stats=byFamily.get(family);
			buffer.clear().append("REFERENCE_FAMILY\t");
			field(inputId).tab();
			field(family).append("\t.\t.");
			for(final IntList counts:selectionCounts){
				final int n=counts.get(i);
				accounted+=n;
				buffer.tab().append(n);
			}
			buffer.tab().append(stats==null ? "UNMAPPED_REFERENCE_FAMILY" : stats.included ? "INCLUDED_PROJECTION" : "EXCLUDED_QUERY_FAMILY").nl();
			flush(out);
		}
		require(accounted==original && rows.size()==kept, "Reference selection must conserve all inspected on-panel rows");
	}

	/** Writes a query outcome; TP counts and endpoint offsets are computed here once. */
	private void writeCall(final Output out, final String inputId, final Call call, final int index,
		final ArrayList<Row> references, final Result match, final Stats stats){
		final int paired=index<0 ? -1 : match.callToReference[index];
		final int masked=index<0 ? -1 : match.callMaskReference[index];
		final String outcome=index<0 ? "EXCLUDED" : paired>=0 ? "TP" : masked>=0 ? "MASKED" : "FP";
		if(index>=0){
			if(paired>=0){stats.tp++;}
			else if(masked>=0){stats.maskCalls++;}
			else{stats.fp++;}
			if(match.duplicateCalls[index]){stats.duplicates++;}
		}
		startLine(inputId, "QUERY", call.locus, outcome);
		buffer.append(".\t").append(call.line).tab();
		field(paired<0 ? "." : references.get(paired).id).append("\t.\t");
		field(masked<0 ? "." : references.get(masked).id).tab();
		if(index<0){buffer.append(".\t.\t");}
		else{buffer.append(match.duplicateCalls[index]).tab().append(match.scoredEligible[index]).tab();}
		field(call.model).tab();
		if(paired>=0){
			final Locus truth=references.get(paired).locus;
			final boolean minus=call.locus.strand=='-';
			final int d5=minus ? truth.stop-call.locus.stop : call.locus.start-truth.start;
			final int d3=minus ? truth.start-call.locus.start : call.locus.stop-truth.stop;
			buffer.append(d5).tab().append(d3).nl();
			if(d5==0){stats.exact5++;}
			if(d3==0){stats.exact3++;}
			if(d5==0 && d3==0){stats.exactBoth++;}
		}else{buffer.append(".\t.\n");}
		flush(out);
	}

	/** Common raw-row prefix; every column after this is supplied by its role writer. */
	private void startLine(final String inputId, final String role, final Locus locus, final String outcome){
		assert(buffer.length()==0) : "A completed report row is flushed before the next original row";
		field(inputId).tab().append(role).tab();
		field(locus.family).tab();
		field(locus.seqid).tab().append(locus.start).tab().append(locus.stop).tab()
			.append(locus.strand).tab().append(outcome).tab();
	}
	private ByteBuilder field(final String text){return BuildNcrnaReference.field(buffer, text);}
	private void flush(final Output out){
		assert(buffer.length()>0) : "A report write must contain the completed row, not an empty placeholder";
		out.writer.print(buffer);
		buffer.clear();
	}
	private static int[] ones(final int count){
		final int[] ids=new int[count];
		Arrays.fill(ids, 1);
		return ids;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Metadata           ----------------*/
	/*--------------------------------------------------------------*/

	public static Rule[] loadRules(final String familymap){
		final ArrayList<Rule> rules=new ArrayList<Rule>();
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(familymap)){
			RrnaResourceIO.header(in.reader, RULE_HEADER, familymap);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==4, "Family map requires exactly four fields");
				final boolean include=p.termEquals("true", 3) || p.termEquals("t", 3);
				require(include || p.termEquals("false", 3) || p.termEquals("f", 3), "Family-map include must be true/t or false/f");
				rules.add(new Rule(text(p, 0), text(p, 1), text(p, 2), include));
			}
		}
		return rules.toArray(new Rule[rules.size()]);
	}

	/** Manifests use literal UTF-8 fields; no path aliasing or percent decoding. */
	public static ArrayList<Binding> loadManifest(final String path){
		final ArrayList<Binding> result=new ArrayList<Binding>();
		final ObjectMap<String, Binding> ids=new ObjectMap<String, Binding>(16, String.class, Binding.class);
		final ObjectMap<String, Binding> fastas=new ObjectMap<String, Binding>(16, String.class, Binding.class);
		final ObjectMap<String, Binding> gffs=new ObjectMap<String, Binding>(16, String.class, Binding.class);
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, MANIFEST_HEADER, path);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==3, "Each grading manifest requires exactly input_id, fasta and gff");
				final Binding binding=new Binding(text(p, 0), text(p, 1), text(p, 2));
				require(!binding.id.isEmpty() && new File(binding.fasta).isAbsolute() && new File(binding.gff).isAbsolute(),
					"Manifest rows require nonempty exact IDs and absolute FASTA/GFF paths");
				require(!ids.contains(binding.id), "Duplicate input_id in grading manifest: "+binding.id);
				require(!fastas.contains(binding.fasta), "Duplicate exact FASTA path in grading manifest: "+binding.fasta);
				require(!gffs.contains(binding.gff), "One GFF cannot stand for multiple declared inputs: "+binding.gff);
				ids.put(binding.id, binding);
				fastas.put(binding.fasta, binding);
				gffs.put(binding.gff, binding);
				result.add(binding);
			}
		}
		return result;
	}
	private void checkBindingOutputs(final Binding binding){
		checkOutputCollision(binding.fasta);
		checkOutputCollision(binding.gff);
	}
	private void checkOutputCollision(final String path){
		for(final String output:outputPaths){
			require(!path.equals(output), "An output cannot overwrite a grading input");
		}
	}
	private static String text(final LineParser1 p, final int column){
		p.setBounds(column);
		final String value=new String(p.line(), p.a(), p.b()-p.a(), StandardCharsets.UTF_8);
		require(value.indexOf('\uFFFD')<0, "Invalid UTF-8 would change grading metadata identity");
		return value;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Reports           ----------------*/
	/*--------------------------------------------------------------*/

	private void writeSummary(final int inputs){
		try(Output out=new Output(prefix+".summary.tsv")){
			out.writer.println(SUMMARY_HEADER);
			for(final Stats stats:families){
				require(stats.tp+stats.fn==stats.scoredReferences, "Matched and missed scored references must conserve the family denominator");
				field(stats.family).tab().append(inputs).tab().append(stats.tp).tab().append(stats.fp).tab().append(stats.fn)
					.tab().append(stats.scoredReferences).tab().append(stats.maskReferences).tab().append(stats.maskCalls)
					.tab().append(stats.duplicates).tab().append(stats.excluded).tab().append(stats.exact5)
					.tab().append(stats.exact3).tab().append(stats.exactBoth).tab();
				ratio(stats.tp, stats.tp+stats.fp).tab();
				ratio(stats.tp, stats.tp+stats.fn).tab().append(overlap, 6).append("\tDECLARED_INPUT_BINDINGS_EXPLICIT_REFERENCE_STATUS\n");
				flush(out);
			}
		}
	}
	private ByteBuilder ratio(final long numerator, final long denominator){
		return denominator==0 ? buffer.append("NA") : buffer.append(numerator/(double)denominator, 6);
	}
	private static void require(final boolean valid, final String reason){
		if(!valid){throw new IllegalArgumentException(reason);}
	}

	/*--------------------------------------------------------------*/
	/*----------------             Types            ----------------*/
	/*--------------------------------------------------------------*/

	public static final class Binding{
		Binding(final String id_, final String fasta_, final String gff_){
			id=id_;
			fasta=fasta_;
			gff=gff_;
		}
		public final String id, fasta, gff;
		boolean matched;
	}
	private static final class Stats{
		Stats(final String family_){family=family_;}
		final String family;
		boolean included;
		long tp, fp, fn, scoredReferences, maskReferences, maskCalls, duplicates, excluded;
		long exact5, exact3, exactBoth;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String input, reference, familymap, prefix;
	private String subtypes, witnesses, events;
	private float overlap=.5f;
	private boolean restrictReference;
	private final PrintStream outstream;
	private final NcrnaReferenceGff codec=new NcrnaReferenceGff();
	private final ByteBuilder buffer=new ByteBuilder();
	private final ObjectMap<String, Stats> byFamily=new ObjectMap<String, Stats>(16, String.class, Stats.class);
	private final ArrayList<Stats> families=new ArrayList<Stats>();
	private final String[] outputPaths;
	private final ObjectIntMap<String> selectionIndex=new ObjectIntMap<String>(String.class);
	private final ArrayList<String> selectionFamilies=new ArrayList<String>();
	private final IntList[] selectionCounts={new IntList(), new IntList(), new IntList(), new IntList()};
	public static final String MANIFEST_HEADER="input_id\tfasta\tgff";
	public static final String RULE_HEADER="family\tmodel_prefix\tdisplay_label\tinclude";
	public static final String RAW_HEADER="input_id\trole\tfamily\tseqid\tstart\tstop\tstrand\toutcome\treference_id\tquery_line\tmatched_reference_id\tmatched_query_line\tmask_reference_id\tduplicate\tscored_eligible\tmodel\td5\td3";
	public static final String INPUT_HEADER="input_id\tfasta\tquery_gff\treference_gff\tfeatures\tcds\ttrna\trrna\tincluded_queries\texcluded_queries";
	public static final String SUMMARY_HEADER="family\tinputs\tTP\tFP\tFN\tscored_references\tmasked_references\tmasked_calls\tduplicate_queries\texcluded_queries\texact_5prime\texact_3prime\texact_both\tprecision\trecall\toverlap\tqualification";
	public static final String SELECTION_HEADER="record\tinput_id\tfamily\tretained_inputs\tomitted_inputs\tscored_kept\tmasked_kept\tscored_omitted\tmasked_omitted\tqualification";
	private static final String[] SUFFIXES={".raw.tsv", ".inputs.tsv", ".summary.tsv", ".selection.tsv"};
}
