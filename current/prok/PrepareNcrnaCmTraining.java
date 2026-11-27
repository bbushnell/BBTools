package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Locale;

import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaReferenceGeometry.Target;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Inventories saved ordinary accepted calls, preserving unmatched calls too.
 * Starts from the full formal roster, never an annotation-positive subset.
 * Reuses completed member-wide assembly/alias/duplicate exclusion audits and
 * the older evaluation-assembly exclusion ledger. No genome is called or aligned.
 * CLI/lifecycle follow template.A_SampleByteFile; native stream owners support
 * multiple complete inputs. Byte-range keys avoid per-query field allocation.
 * @author Brian Bushnell, Raiden
 */
public final class PrepareNcrnaCmTraining{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, PrepareNcrnaCmTraining.class, false);
		try{new PrepareNcrnaCmTraining(pp.args).process();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private PrepareNcrnaCmTraining(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("bindings")){bindings=b;}
			else if(a.equals("audits")){audits=b.split(",", -1);}
			else if(a.equals("exclude")){exclude=b;}
			else if(a.equals("refs")){refs=b.split(",", -1);}
			else if(a.equals("expectedmodels")){expectedModels=Parse.parseIntKMG(b);}
			else if(a.equals("expectedinputs")){expectedInputs=Parse.parseIntKMG(b);}
			else if(a.equals("expectedprefixes")){expectedPrefixes=Parse.parseIntKMG(b);}
			else if(a.equals("mincalls")){minCalls=Parse.parseIntKMG(b);}
			else if(a.equals("minassemblies")){minAssemblies=Parse.parseIntKMG(b);}
			else{require(parser.parse(arg, a, b), "Unknown CM training-input option: "+arg);}
		}
		prefixList=parser.in1;
		out=parser.out1;
		require(prefixList!=null && bindings!=null && audits!=null && audits.length==3
			&& refs!=null && refs.length==3 && exclude!=null && out!=null
			&& expectedModels>0 && expectedInputs>0 && expectedPrefixes>0 && minCalls>0 && minAssemblies>0
			&& !parser.overwrite && !parser.append && parser.in2==null && parser.out2==null
			&& parser.maxReads<0 && parser.samplerate==1,
			"Bind complete saved prefix list, full bindings, three member audits/refs, evaluation exclude ledger and fresh out prefix");
		RrnaResourceIO.fresh(out, ".calls.tsv.gz", ".sources.tsv", ".models.tsv", ".summary.tsv");
	}

	/*--------------------------------------------------------------*/
	/*----------------          Processing          ----------------*/
	/*--------------------------------------------------------------*/

	private void process(){
		readSources();
		for(int member=0; member<3; member++){
			readAudit(audits[member], member==0);
			readModels(refs[member], member);
		}
		require(models.size()==expectedModels, "Consensus model count differs from the declared frozen model axis");
		readExclusions();
		writeSources();
		readCalls();
		writeCensus();
		System.err.println("NCRNA_CM_TRAINING_INPUTS_PASS_NO_CALLER_NO_CM_SEARCH");
	}

	private void readSources(){
		try(Input input=new Input(bindings)){
			RrnaResourceIO.header(input.reader, NcrnaDdlHoldout.PANEL_HEADER, bindings);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==5 && p.parseInt(0)==sources.size(), "Full source roster needs dense ordinals and five fields");
				final Source source=new Source(sources.size(), p.parseString(1), p.parseString(2), p.parseString(3), p.parseString(4));
				probe.set(p, 1);
				require(!byAssembly.contains(probe), "Formal assemblies must be distinct");
				byAssembly.put(probe.copy(), source);
				probe.set(p, 2);
				require(!byInput.contains(probe), "Saved input IDs must bind exactly one formal assembly");
				byInput.put(probe.copy(), source);
				sources.add(source);
			}
		}
		require(sources.size()==expectedInputs, "Full source roster differs from expectedinputs; never substitute a positive subset");
		blocked=new boolean[sources.size()];
	}

	private void readAudit(final String path, final boolean first){
		int index=0;
		final boolean[] seen=new boolean[sources.size()];
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, AUDIT_HEADER, path);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(index<sources.size() && p.terms()==8, "Holdout audit must conserve the full formal source axis");
				probe.set(p, 0);
				final Source source=byAssembly.get(probe);
				require(source!=null && !seen[source.ordinal] && p.termEquals(source.input, 1),
					"Every audit assembly/input must bind once; formal and shuffled caller orders may differ");
				seen[source.ordinal]=true;
				index++;
				final int group=p.parseInt(4);
				require(group>=0 && group<blocked.length && (first || source.group==group), "Member audits must share the same resolved assembly/alias groups");
				source.group=group;
				blocked[group]|=booleanField(5) | booleanField(6);
			}
		}
		require(index==sources.size(), "Truncated holdout audit cannot establish training eligibility");
	}

	private void readModels(final String path, final int member){
		final NcrnaSupport.Library library=NcrnaSupport.loadLibrary(path);
		for(int local=0; local<library.names.length; local++){
			final String name=library.names[local];
			require(name.startsWith("cobalamin_"+MEMBERS[member]+"_"), "Consensus order/member identity must match the frozen caller libraries");
			final Model model=new Model(models.size(), member, name, sources.size());
			final byte[] bytes=name.getBytes(StandardCharsets.UTF_8);
			probe.set(bytes, 0, bytes.length);
			require(!byModel.contains(probe), "Consensus identifiers must be unique across the clan");
			byModel.put(probe.copy(), model);
			models.add(model);
		}
	}

	private void readExclusions(){
		try(Input input=new Input(exclude)){
			RrnaResourceIO.header(input.reader, EXCLUSION_LABEL, exclude);
			RrnaResourceIO.header(input.reader, EXCLUSION_HEADER, exclude);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==5, "Existing evaluation exclusion ledger requires five fields");
				if(p.termEquals("UNMATCHED_MANIFEST_PATH_ALIAS", 4)){
					require(p.termEquals(".", 3), "Unmatched historical evaluation source must not claim a formal assembly");
					unmatchedExclusions++;
					continue;
				}
				require(p.termEquals("EXCLUDED_ASSEMBLY", 4), "Unknown existing evaluation exclusion disposition");
				probe.set(p, 3);
				final Source source=byAssembly.get(probe);
				require(source!=null, "Evaluation assembly must exist in the full bound roster");
				blocked[source.group]=true;
			}
		}
	}

	private void writeSources(){
		try(Output output=new Output(out+".sources.tsv")){
			output.writer.println(NcrnaDdlHoldout.PANEL_HEADER);
			for(final Source source:sources){
				source.allowed=!blocked[source.group];
				if(!source.allowed){continue;}
				output.writer.print(row.clear().append(trainingInputs++).tab().append(source.assembly).tab().append(source.input)
					.tab().append(source.fasta).tab().append(source.gff).nl());
			}
		}
		require(trainingInputs>0, "No nonexcluded training assembly remains");
	}

	private void readCalls(){
		int prefixes=0;
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		try(Input list=new Input(prefixList); Output output=new Output(out+".calls.tsv.gz")){
			RrnaResourceIO.header(list.reader, "prefix", prefixList);
			output.writer.println(CALL_HEADER);
			for(byte[] line; (line=list.reader.nextLine())!=null;){
				final String prefix=new String(line, StandardCharsets.UTF_8);
				require(!prefix.isEmpty() && seen.add(prefix), "Saved ordinary-call prefixes must be nonempty and distinct");
				readQueries(prefix, prefixes++, output);
			}
		}
		require(prefixes==expectedPrefixes && selectedCalls>0, "Saved prefix census must be complete and retain actual accepted calls");
	}

	private void readQueries(final String prefix, final int owner, final Output output){
		try(Input input=new Input(prefix+".queries.tsv")){
			RrnaResourceIO.header(input.reader, QUERY_HEADER, prefix);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==10 && p.parseInt(2)>0 && p.parseInt(3)>=p.parseInt(2)
					&& (p.termEquals("+", 4) || p.termEquals("-", 4)), "Saved accepted-call coordinates and strand must be intact");
				probe.set(p, 0);
				final Source source=byInput.get(probe);
				require(source!=null && (source.owner<0 || source.owner==owner), "Every saved caller source must belong to exactly one disjoint prefix");
				source.owner=owner;
				p.setBounds(6);
				int stop=p.a();
				while(stop<p.b() && line[stop]>32){stop++;}
				probe.set(line, p.a(), stop);
				final Model model=byModel.get(probe);
				require(model!=null && p.termEquals(MEMBERS[model.member], 5), "Saved winner must bind a real frozen consensus model and matching member");
				allCalls++;
				if(!source.allowed){excludedCalls++; continue;}
				model.calls++;
				if(model.addAssembly(source.group)){model.assemblies++;}
				if(p.termEquals("TP", 7)){matchedCalls++;}
				output.writer.print(row.clear().append(selectedCalls++).tab().append(source.assembly).tab().append(model.ordinal)
					.tab().append(line).nl());
			}
		}
	}

	private void writeCensus(){
		int supported=0, below=0;
		long pathCalls=0;
		try(Output output=new Output(out+".models.tsv")){
			output.writer.println("model_ordinal0\tmember\tmodel\taccepted_calls\tresolved_assembly_groups\tlearned_corridor_supported");
			for(final Model model:models){
				if(model.member<2){pathCalls+=model.calls;}
				final boolean enough=model.member<2 && model.calls>=minCalls && model.assemblies>=minAssemblies;
				if(model.member<2){if(enough){supported++;}else{below++;}}
				output.writer.print(row.clear().append(model.ordinal).tab().append(MEMBERS[model.member]).tab().append(model.name)
					.tab().append(model.calls).tab().append(model.assemblies).tab().append(enough).nl());
			}
		}
		require(allCalls==selectedCalls+excludedCalls, "Training and excluded accepted calls must conserve the complete saved population");
		try(Output output=new Output(out+".summary.tsv")){
			output.writer.println("formal_inputs\ttraining_inputs\tall_accepted_calls\texcluded_calls\ttraining_calls\tmatched_training_calls\tmain_second_supported_models\tmain_second_below_support_models\tunmatched_historical_exclusions\tmain_second_path_calls");
			output.writer.print(row.clear().append(sources.size()).tab().append(trainingInputs).tab().append(allCalls).tab().append(excludedCalls)
				.tab().append(selectedCalls).tab().append(matchedCalls).tab().append(supported).tab().append(below).tab().append(unmatchedExclusions).tab().append(pathCalls).nl());
		}
	}

	private boolean booleanField(final int column){
		require(p.termEquals("true", column) || p.termEquals("false", column), "Assembly exclusion flags must be explicit booleans");
		return p.termEquals("true", column);
	}
	private static void require(final boolean ok, final String why){if(!ok){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Source{
		Source(final int ordinal_, final String assembly_, final String input_, final String fasta_, final String gff_){
			ordinal=ordinal_;
			assembly=assembly_;
			input=input_;
			fasta=fasta_;
			gff=gff_;
		}
		final int ordinal;
		final String assembly, input, fasta, gff;
		int group=-1, owner=-1;
		boolean allowed;
	}
	private static final class Model{
		Model(final int ordinal_, final int member_, final String name_, final int sources){
			ordinal=ordinal_;
			member=member_;
			name=name_;
			assemblyBits=new long[(sources+63)/64];
		}
		boolean addAssembly(final int group){
			assert(group>=0 && (group>>>6)<assemblyBits.length) : "Resolved alias-group ordinal must fit the per-model distinct-assembly bitset before indexing: "+group;
			final int word=group>>>6;
			final long mask=1L<<(group&63);
			final boolean absent=(assemblyBits[word]&mask)==0;
			assemblyBits[word]|=mask;
			return absent;
		}
		final int ordinal, member;
		final String name;
		final long[] assemblyBits;
		long calls;
		int assemblies;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String prefixList, bindings, exclude, out;
	private String[] audits, refs;
	private int expectedModels, expectedInputs, expectedPrefixes, minCalls=100, minAssemblies=10, trainingInputs;
	private long allCalls, excludedCalls, selectedCalls, matchedCalls, unmatchedExclusions;
	private boolean[] blocked;
	private final ArrayList<Source> sources=new ArrayList<Source>();
	private final ArrayList<Model> models=new ArrayList<Model>();
	private final ObjectMap<Target, Source> byInput=new ObjectMap<Target, Source>(Target.class, Source.class);
	private final ObjectMap<Target, Source> byAssembly=new ObjectMap<Target, Source>(Target.class, Source.class);
	private final ObjectMap<Target, Model> byModel=new ObjectMap<Target, Model>(Target.class, Model.class);
	private final Target probe=new Target();
	private final LineParser1 p=new LineParser1('\t');
	private final ByteBuilder row=new ByteBuilder();
	private static final String[] MEMBERS={"RF00174", "RF01482", "RF01689"};
	private static final String AUDIT_HEADER="assembly_id\tinput_id\tcanonical_fasta\thistorical_alias_group\tresolved_group\theld_out\tduplicate_original_group\tsampled";
	private static final String EXCLUSION_LABEL="# HELD OUT (panel genomes excluded from tables and training): exclusion bindings; not a fitted result";
	private static final String EXCLUSION_HEADER="panel\trow1\tfasta\tassembly\tstatus";
	public static final String QUERY_HEADER="input_id\tseqid\tstart\tstop\tstrand\temitted_member\tmodel\toutcome\tscored_eligible\tmatched_reference_id";
	public static final String CALL_HEADER="training_call_ordinal0\tassembly_id\tmodel_ordinal0\t"+QUERY_HEADER;
}
