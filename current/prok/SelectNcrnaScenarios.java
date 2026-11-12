package prok;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Locale;

import map.ObjectIntMap;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaReferenceCalls.Call;
import prok.NcrnaReferenceCalls.Rule;
import prok.NcrnaReferenceGeometry.Target;
import prok.NcrnaTextIO.ConfigInput;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/** Selects identical input/family scopes from paired hypothetical scenarios.
 * Derived from A_SampleByteFile. Does not alter a disposition or infer a biological
 * family. An explicit original caller config and complete adapter bindings protect
 * zero-call inputs. Gate interpretation is deliberately limited to unambiguous,
 * explicitly named canonical flags, not current caller defaults or arbitrary
 * configuration side effects. The admitted profile requires euk=f.
 *
 * Paired rows must have the same input-local decision key and family in the same
 * order. All rows remain in global companions; subsets preserve included row bytes.
 * A missing/disabled/unsupported population is counted separately from a complete
 * input with zero calls. One thread owns metadata, parsers and primitive counters;
 * only one input's parsed query calls are retained at a time.
 * @author Brian Bushnell, Raiden
 */
public final class SelectNcrnaScenarios{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final PreParser pp=new PreParser(args, SelectNcrnaScenarios.class, false);
		try{new SelectNcrnaScenarios(pp.args, pp.outstream).process(timer);}
		finally{Shared.closeStream(pp.outstream);}
	}
	private SelectNcrnaScenarios(final String[] args, final PrintStream stream){
		outstream=stream;
		Shared.setThreads(1);
		final Parser common=new Parser();
		common.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("calls")){callsPath=b;}
			else if(a.equals("bindings")){bindingsPath=b;}
			else if(a.equals("models")){modelsPath=b;}
			else if(a.equals("scope")){scopePath=b;}
			else if(a.equals("callerconfig")){configPath=b;}
			else if(!common.parse(arg, a, b)){throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		first=common.in1;
		second=common.in2;
		require(first!=null && second!=null && callsPath!=null && bindingsPath!=null && modelsPath!=null && scopePath!=null && configPath!=null && common.out1!=null,
			"Supply in=GA in2=margin calls=manifest bindings=adapter_bindings models=prefix_map scope=flags callerconfig=original out=new_directory");
		require(!common.overwrite && !common.append, "Scenario subsets require fresh outputs");
		require(Tools.testInputFiles(false, true, first, second, callsPath, bindingsPath, modelsPath, scopePath, configPath), "Selection inputs must be readable");
		output=new File(common.out1).getAbsoluteFile().toPath();
		require(!Files.exists(output), "Selection output directory must not exist");
	}

	private void process(final Timer timer){
		assert(families.isEmpty() && inputs.isEmpty()) : "A selector instance owns exactly one fresh output population";
		readScope();
		readConfig();
		readModels();
		readInputs();
		validateCalls();
		try{Files.createDirectory(output);}
		catch(IOException failure){throw new IllegalStateException("Could not create scenario subset directory", failure);}
		writeMetadata();
		selectRows();
		writeCounts();
		writeSummary();
		timer.stop();
		outstream.println("Selected "+selected+" paired decisions for "+inputs.size()+" complete declared inputs; time: "+timer);
		outstream.println(QUALIFICATION);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Explicit scope       ----------------*/
	/*--------------------------------------------------------------*/

	/** Loads retained family/gate metadata; no per-event text is allocated here. */
	private void readScope(){
		try(Input in=new Input(scopePath)){
			RrnaResourceIO.header(in.reader, SCOPE_HEADER, scopePath);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				left.set(row);
				require(left.terms()==3, "Scope rows require family, canonical required flags and supported status");
				final String name=literal(left, 0);
				require(!name.isEmpty() && !byName.contains(name), "Scope family must be named exactly once");
				final String[] flags=literal(left, 1).split(",", -1);
				final ObjectSet<String> seen=new ObjectSet<String>(String.class);
				for(final String flag:flags){
					require(flag.matches("[a-z][a-z0-9]*") && seen.add(flag), "Scope flags must be nonempty unique canonical lowercase names");
					required.put(flag, 0);
				}
				final Family family=new Family(name, flags, bool(left, 2), true);
				register(family);
			}
		}
		require(!families.isEmpty(), "At least one explicitly supported or unsupported family is required");
		required.put("euk", 0);
	}

	/** This is an explicit gate declaration, never a replay of arbitrary CallGenes args. */
	private void readConfig(){
		try(ConfigInput in=new ConfigInput(configPath)){
			while(in.next()){
				final String key=in.key.toLowerCase(Locale.ROOT);
				require(!key.equals("config") && !key.equals("in2"), "Nested or paired caller configs need a separately resolved binding");
				if(key.equals("in") || key.equals("in1")){
					require(new File(in.value).isAbsolute(), "Original config inputs must be absolute single-FASTA paths");
					originalInputs.add(in.value);
					continue;
				}
				require(!key.endsWith("only") && !unsupportedGate(key), "Unsupported caller alias/only-mode can change family enablement");
				if(!required.contains(key)){continue;}
				require(!seenFlags.contains(key), "Missing/duplicate/conflicting canonical gate flags cannot be inferred");
				seenFlags.add(key);
				final boolean value=parseBoolean(in.value);
				required.put(key, value ? 1 : -1);
			}
		}
		require(seenFlags.size()==required.size() && required.get("euk")==-1 && !originalInputs.isEmpty(),
			"Every canonical gate must be explicit once, euk=f, and the complete original input list present");
		for(final Family family:families){
			family.enabled=true;
			for(final String flag:family.flags){family.enabled&=required.get(flag)==1;}
		}
	}
	private static boolean unsupportedGate(final String key){
		return key.equals("generalncrna") || key.equals("ssrs") || key.equals("6s") || key.equals("b12")
			|| key.equals("r58lsu") || key.equals("s18") || key.equals("euk18s") || key.equals("euk5s")
			|| key.equals("rrna17") || key.equals("rrna17mer");
	}

	/** Binds explicit query prefixes to the fixed family selection, including exclusions. */
	private void readModels(){
		try(Input in=new Input(modelsPath)){
			RrnaResourceIO.header(in.reader, MODEL_HEADER, modelsPath);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				left.set(row);
				require(left.terms()==3, "Model mapping needs family, model prefix and exact display label");
				final Family family=byName.get(literal(left, 0));
				require(family!=null, "Every caller model needs an explicit baseline family scope");
				final Rule rule=new Rule(family.name, literal(left, 1), literal(left, 2), family.enabled && family.supported);
				rules.add(rule);
				family.modelRules++;
			}
		}
		for(final Family family:families){require(family.modelRules>0 || !family.supported, "A supported scoring family needs an explicit caller-model mapping");}
		reader=new NcrnaReferenceCalls(rules.toArray(new Rule[rules.size()]));
	}

	/*--------------------------------------------------------------*/
	/*----------------      Complete input binding  ----------------*/
	/*--------------------------------------------------------------*/

	/** Requires an ordered one-to-one projection of the complete original input list. */
	private void readInputs(){
		final ObjectSet<String> fastas=new ObjectSet<String>(String.class), gffs=new ObjectSet<String>(String.class);
		try(Input calls=new Input(callsPath); Input bindings=new Input(bindingsPath)){
			RrnaResourceIO.header(calls.reader, GradeNcrnaReference.MANIFEST_HEADER, callsPath);
			RrnaResourceIO.header(bindings.reader, BINDING_HEADER, bindingsPath);
			int index=0;
			for(byte[] a=calls.reader.nextLine(), b=bindings.reader.nextLine(); a!=null || b!=null; a=calls.reader.nextLine(), b=bindings.reader.nextLine()){
				require(a!=null && b!=null && index<originalInputs.size(), "Call manifest must retain the complete adapter/config population, including empty queries");
				left.set(a);
				right.set(b);
				require(left.terms()==3 && right.terms()==8, "Query and adapter binding schemas changed");
				final String id=literal(left, 0), fasta=literal(left, 1), gff=literal(left, 2);
				final String original=decoded(right, 3);
				require(!id.isEmpty() && id.equals(decoded(right, 0)) && fasta.equals(decoded(right, 4)) && gff.equals(decoded(right, 6))
					&& original.equals(originalInputs.get(index)) && right.parseLong(5)==index+2L
					&& right.termEquals("GENERIC_RNA_ONLY_ORIGINAL_OTHER_FEATURES_RETAINED", 7), "Original/canonical FASTA aliases and ordered query bindings must agree exactly");
				require(new File(fasta).isAbsolute() && new File(gff).isAbsolute() && fastas.add(fasta) && gffs.add(gff), "Each canonical FASTA/query path must be absolute and unique");
				final Binding binding=new Binding(id, original, fasta, gff);
				require(!byInput.contains(binding.key), "Duplicate query input identity");
				byInput.put(binding.key, binding);
				inputs.add(binding);
				index++;
			}
			// The original config, not an observed nonzero call, defines the denominator.
			require(index==originalInputs.size(), "Shortened manifests cannot remove zero-call inputs from the original config");
		}
	}

	/** Checks every projected RNA row; unknown attribution or disabled output is fatal. */
	private void validateCalls(){
		assert(reader!=null) : "Query validation must use the loaded historical prefix mapping";
		for(final Binding input:inputs){
			final NcrnaReferenceCalls.Result calls=reader.read(input.gff, input.id);
			require(calls.features==calls.calls.size(), "Adapter query files must contain generic RNA only; original other types stay in the source inventory");
			input.calls=calls.calls.size();
			if(input.calls==0){emptyInputs++;}
			for(final Call call:calls.calls){
				final Family family=byName.get(call.locus.family);
				require(family!=null && family.enabled, "Saved query contradicts a disabled family in the original caller configuration");
				family.calls++;
				if(call.included){input.included++;}
				else{input.excluded++;}
			}
			require(input.calls==input.included+input.excluded, "Every generic query must retain exactly one included/excluded attribution");
			queryCalls+=input.calls;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Paired selection      ----------------*/
	/*--------------------------------------------------------------*/

	/** Streams both sides together and retains one byte key per global decision.
	 * Canonical TSV escaping makes raw-key equality equivalent to decoded equality.
	 * Excluded rows also participate in duplicate detection and conservation.
	 */
	private void selectRows(){
		final ObjectSet<Target> decisions=new ObjectSet<Target>(Target.class);
		try(Input ga=new Input(first); Input margin=new Input(second); Output low=new Output(path("ga.tsv.gz")); Output high=new Output(path("margin.tsv.gz"))){
			RrnaResourceIO.header(ga.reader, CompareNcrnaScenarios.HEADER, first);
			RrnaResourceIO.header(margin.reader, CompareNcrnaScenarios.HEADER, second);
			low.writer.println(CompareNcrnaScenarios.HEADER);
			high.writer.println(CompareNcrnaScenarios.HEADER);
			for(byte[] a=ga.reader.nextLine(), b=margin.reader.nextLine(); a!=null || b!=null; a=ga.reader.nextLine(), b=margin.reader.nextLine()){
				require(a!=null && b!=null, "Asymmetric scenario populations cannot define a common baseline subset");
				left.set(a);
				right.set(b);
				require(left.terms()==14 && right.terms()==14 && same(left, 0, right, 0) && same(left, 1, right, 1) && same(left, 6, right, 6),
					"Paired scenarios require identical ordered input_id/decision_id/family keys");
				canonicalKey(left, 0);
				canonicalKey(left, 1);
				canonicalKey(left, 6);
				left.setBounds(1);
				final int end=left.b();
				require(left.a()<end, "Every scenario decision needs an input identity");
				left.setBounds(0);
				require(left.a()<left.b(), "Every scenario row needs a decision identity");
				// Retain one framed key per pair, including exclusions; never a field String.
				probe.set(a, left.a(), end);
				require(decisions.add(probe.copy()), "Duplicate input-local decision in global scenarios");
				probe.set(left, 1);
				final Binding input=byInput.get(probe);
				probe.set(left, 6);
				Family family=byFamily.get(probe);
				if(family==null){
					family=new Family(decoded(left, 6), new String[0], false, false);
					register(family);
				}
				final int before=state(left), after=state(right);
				require(family.defined || before>=UNRESOLVED && after>=UNRESOLVED,
					"A scored/masked global family needs explicit baseline scope, even when unsupported");
				final int selection=input==null ? OUTSIDE_INPUT : !family.defined || !family.supported ? UNSUPPORTED
					: !family.enabled ? DISABLED : INCLUDED;
				family.decisions[selection][before][after]++;
				globalRows++;
				if(selection!=INCLUDED){continue;}
				low.writer.println(a);
				high.writer.println(b);
				input.decisions++;
				selected++;
			}
		}
		assert(decisions.size()==globalRows) : "Unique paired keys must conserve all rows before subset filtering";
	}

	/*--------------------------------------------------------------*/
	/*----------------             Reports          ----------------*/
	/*--------------------------------------------------------------*/

	/** Records original companions and emits the comparator's explicit family map. */
	private void writeMetadata(){
		try(Output out=new Output(path("sources.tsv"))){
			out.writer.println("role\tpath\tqualification");
			final String[] roles={"global_ga", "global_margin", "complete_query_manifest", "adapter_bindings", "model_map", "family_scope", "original_caller_config"};
			final String[] paths={first, second, callsPath, bindingsPath, modelsPath, scopePath, configPath};
			for(int i=0; i<roles.length; i++){
				buffer.clear().append(roles[i]).tab();
				field(new File(paths[i]).getAbsolutePath()).tab().append(QUALIFICATION).nl();
				out.writer.print(buffer);
			}
		}
		try(Output out=new Output(path("familymap.tsv"))){
			out.writer.println(GradeNcrnaReference.RULE_HEADER);
			for(final Rule rule:rules){
				buffer.clear();
				field(rule.family).tab();
				field(rule.modelPrefix).tab();
				field(rule.displayLabel).tab().append(rule.include).nl();
				out.writer.print(buffer);
			}
		}
	}
	/** Publishes reason/state cross-counts and every input, including empty queries. */
	private void writeCounts(){
		long counted=0, countedIncluded=0;
		try(Output out=new Output(path("selection_counts.tsv"))){
			out.writer.println("family\tselection\tbefore_state\tafter_state\tpaired_rows\tqualification");
			for(final Family family:families){
				for(int reason=0; reason<SELECTIONS.length; reason++){
					for(int before=0; before<STATES.length; before++){
						for(int after=0; after<STATES.length; after++){
							final long count=family.decisions[reason][before][after];
							if(count==0){continue;}
							counted+=count;
							if(reason==INCLUDED){countedIncluded+=count;}
							buffer.clear();
							field(family.name).tab().append(SELECTIONS[reason]).tab().append(STATES[before]).tab().append(STATES[after])
								.tab().append(count).tab().append(QUALIFICATION).nl();
							out.writer.print(buffer);
						}
					}
				}
			}
		}
		require(counted==globalRows && countedIncluded==selected, "Inclusion/exclusion cells must conserve every paired global decision");
		try(Output out=new Output(path("inputs.tsv"))){
			out.writer.println("input_id\toriginal_fasta\tcanonical_fasta\tquery_gff\tgeneric_RNA\tincluded_calls\texcluded_calls\tselected_decisions\tqualification");
			for(final Binding input:inputs){
				buffer.clear();
				field(input.id).tab();
				field(input.original).tab();
				field(input.fasta).tab();
				field(input.gff).tab().append(input.calls).tab().append(input.included).tab().append(input.excluded)
					.tab().append(input.decisions).tab().append(QUALIFICATION).nl();
				out.writer.print(buffer);
			}
		}
		try(Output out=new Output(path("scope.tsv"))){
			out.writer.println("family\tdeclared_scope\tenabled\tsupported\tquery_calls\tqualification");
			for(final Family family:families){
				buffer.clear();
				field(family.name).tab().append(family.defined).tab().append(family.enabled).tab().append(family.supported)
					.tab().append(family.calls).tab().append(QUALIFICATION).nl();
				out.writer.print(buffer);
			}
		}
	}
	/** Written only after successful row selection and closure of all prior outputs. */
	private void writeSummary(){
		long decisions=0;
		for(final Binding input:inputs){decisions+=input.decisions;}
		require(decisions==selected && inputs.size()==originalInputs.size(), "Complete original inputs and selected decisions must conserve the baseline denominator");
		try(Output out=new Output(path("summary.tsv"))){
			out.writer.println("inputs\tempty_query_inputs\tgeneric_RNA\tglobal_paired_decisions\tselected_paired_decisions\texcluded_paired_decisions\tqualification");
			buffer.clear().append(inputs.size()).tab().append(emptyInputs).tab().append(queryCalls).tab().append(globalRows)
				.tab().append(selected).tab().append(globalRows-selected).tab().append(QUALIFICATION).nl();
			out.writer.print(buffer);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------             Helpers          ----------------*/
	/*--------------------------------------------------------------*/

	private void register(final Family family){
		assert(family!=null) : "Both declared and context-only family counters require retained metadata";
		require(!byFamily.contains(family.key), "Family names must remain unique after evidence TSV encoding");
		byFamily.put(family.key, family);
		byName.put(family.name, family);
		families.add(family);
	}
	private static int state(final LineParser1 parser){
		for(int i=0; i<STATES.length; i++){if(parser.termEquals(STATES[i], 8)){return i;}}
		throw new IllegalArgumentException("Unknown scenario disposition");
	}
	private static boolean same(final LineParser1 a, final int ac, final LineParser1 b, final int bc){
		assert(ac>=0 && ac<a.terms() && bc>=0 && bc<b.terms()) : "Pair comparison requires validated scenario columns";
		a.setBounds(ac);
		b.setBounds(bc);
		if(a.b()-a.a()!=b.b()-b.a()){return false;}
		for(int i=0; i<a.b()-a.a(); i++){if(a.line()[a.a()+i]!=b.line()[b.a()+i]){return false;}}
		return true;
	}
	/** BuildNcrnaReference.field escapes only percent, tab, CR and LF, in uppercase. */
	private static void canonicalKey(final LineParser1 parser, final int column){
		parser.setBounds(column);
		final byte[] data=parser.line();
		for(int i=parser.a(); i<parser.b(); i++){
			if(data[i]!='%'){continue;}
			require(i+2<parser.b() && ((data[i+1]=='2' && data[i+2]=='5')
				|| (data[i+1]=='0' && (data[i+2]=='9' || data[i+2]=='A' || data[i+2]=='D'))),
				"Scenario identity requires canonical evidence TSV escaping");
			i+=2;
		}
	}
	private String decoded(final LineParser1 parser, final int column){
		parser.setBounds(column);
		return decoder.decode(parser.line(), parser.a(), parser.b());
	}
	private static String literal(final LineParser1 parser, final int column){
		assert(column>=0 && column<parser.terms()) : "Metadata text must remain inside its validated TSV schema";
		parser.setBounds(column);
		final String text=new String(parser.line(), parser.a(), parser.b()-parser.a(), StandardCharsets.UTF_8);
		require(text.indexOf('\uFFFD')<0, "Malformed UTF-8 changes declared input/family identity");
		return text;
	}
	private static Target key(final String text){
		final byte[] data=BuildNcrnaReference.field(new ByteBuilder(), text).toBytes();
		final Target probe=new Target();
		probe.set(data, 0, data.length);
		return probe.copy();
	}
	private static boolean bool(final LineParser1 parser, final int column){return parseBoolean(literal(parser, column));}
	private static boolean parseBoolean(final String text){
		if(text.equalsIgnoreCase("true") || text.equalsIgnoreCase("t")){return true;}
		require(text.equalsIgnoreCase("false") || text.equalsIgnoreCase("f"), "Declared boolean flags require true/t or false/f");
		return false;
	}
	private ByteBuilder field(final String text){return BuildNcrnaReference.field(buffer, text);}
	private String path(final String file){return output.resolve(file).toString();}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	private static final class Family{
		Family(final String name_, final String[] flags_, final boolean supported_, final boolean defined_){
			name=name_;
			flags=flags_;
			supported=supported_;
			defined=defined_;
			key=key(name);
		}
		final String name;
		final String[] flags;
		final boolean supported, defined;
		final Target key;
		boolean enabled;
		int modelRules;
		long calls;
		final long[][][] decisions=new long[SELECTIONS.length][STATES.length][STATES.length];
	}
	private static final class Binding{
		Binding(final String id_, final String original_, final String fasta_, final String gff_){
			id=id_;
			original=original_;
			fasta=fasta_;
			gff=gff_;
			key=key(id);
		}
		final String id, original, fasta, gff;
		final Target key;
		long calls, included, excluded, decisions;
	}
	private String first, second, callsPath, bindingsPath, modelsPath, scopePath, configPath;
	private final Path output;
	private final PrintStream outstream;
	private NcrnaReferenceCalls reader;
	private final LineParser1 left=new LineParser1('\t'), right=new LineParser1('\t');
	private final NcrnaHistoricalEvidence.Attributes decoder=new NcrnaHistoricalEvidence.Attributes();
	private final Target probe=new Target();
	private final ByteBuilder buffer=new ByteBuilder();
	private final ArrayList<Family> families=new ArrayList<Family>();
	private final ArrayList<Binding> inputs=new ArrayList<Binding>();
	private final ArrayList<Rule> rules=new ArrayList<Rule>();
	private final ArrayList<String> originalInputs=new ArrayList<String>();
	private final ObjectIntMap<String> required=new ObjectIntMap<String>(String.class);
	private final ObjectSet<String> seenFlags=new ObjectSet<String>(String.class);
	private final ObjectMap<String,Family> byName=new ObjectMap<String,Family>(String.class, Family.class);
	private final ObjectMap<Target,Family> byFamily=new ObjectMap<Target,Family>(Target.class, Family.class);
	private final ObjectMap<Target,Binding> byInput=new ObjectMap<Target,Binding>(Target.class, Binding.class);
	private long globalRows, selected, emptyInputs, queryCalls;
	private static final int INCLUDED=0, OUTSIDE_INPUT=1, UNSUPPORTED=2, DISABLED=3, UNRESOLVED=2;
	private static final String[] SELECTIONS={"INCLUDED", "OUTSIDE_BASELINE_INPUT", "UNSUPPORTED_OR_UNDECLARED_FAMILY", "DISABLED_FAMILY"};
	private static final String[] STATES={"scored", "dont_care", "unresolved", "removed", "context"};
	public static final String MODEL_HEADER="family\tmodel_prefix\tdisplay_label";
	public static final String SCOPE_HEADER="family\trequired_flags\tsupported";
	public static final String BINDING_HEADER="input_id\tpanel_genome\tassembly\toriginal_fasta\tcatalog_fasta\tpanel_line\tquery_gff\tprojection";
	public static final String QUALIFICATION="DECLARED_BASELINE_SCOPE_HYPOTHETICAL_REFERENCE_NOT_BIOLOGICAL_ACCURACY";
}
