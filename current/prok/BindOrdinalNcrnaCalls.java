package prok;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Locale;

import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.ConfigInput;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import structures.ByteBuilder;

/** Binds an ordered, input-ordinal-qualified saved caller GFF without reading DNA.
 * Each source ordinal maps to the declared global reference input ID, not a contig
 * name, assembly name or taxID. All feature rows and within-input order survive;
 * normalized per-file headers and a source-line ledger make the projection explicit.
 * No family attribution, biological labeling or matching occurs here.
 * Native streams and lifecycle follow template.A_SampleByteFile. One consumer and
 * one per-input writer suffice because boundaries must be contiguous from zero.
 * @author Raiden
 */
public final class BindOrdinalNcrnaCalls{

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final PreParser pp=new PreParser(args, BindOrdinalNcrnaCalls.class, false);
		try{
			final BindOrdinalNcrnaCalls tool=new BindOrdinalNcrnaCalls(pp.args, pp.outstream);
			tool.process();
			timer.stop();
			pp.outstream.println("ORDINAL_NCRNA_CALLS_BOUND time="+timer);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private BindOrdinalNcrnaCalls(final String[] args, final PrintStream stream){
		outstream=stream;
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("bindings")){
				bindings=b;
			}else if(a.equals("callerconfig")){
				callerConfig=b;
			}else if(a.equals("expectedinputs")){
				expectedInputs=Parse.parseIntKMG(b);
			}else{require(parser.parse(arg, a, b), "Unknown ordinal-binding option: "+arg);}
		}
		input=parser.in1;
		directory=parser.out1;
		require(input!=null && bindings!=null && callerConfig!=null && directory!=null
			&& expectedInputs>0,
			"Supply in=saved_GFF bindings=panel.sources.tsv callerconfig=actual_flat_config expectedinputs=N out=fresh_directory");
		require(!parser.overwrite && !parser.append && parser.in2==null && parser.out2==null
			&& parser.maxReads<0 && parser.samplerate==1, "Binding must retain the entire saved query population without overwrite/append");
		require(new File(input).isAbsolute() && !new File(directory).exists(), "Use an absolute saved GFF and a fresh output directory");
		Shared.setThreads(1);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	private void process(){
		loadBindings();
		checkConfig();
		require(new File(directory).mkdir() && new File(directory, "queries").mkdir(), "Could not create fresh binding outputs");
		try(Input source=new Input(input); Output queries=output("calls.tsv"); Output references=output("references.tsv");
			Output origins=output("source_lines.tsv"); Output counts=output("inputs.tsv")){
			queries.writer.println(GradeNcrnaReference.MANIFEST_HEADER);
			references.writer.println(GradeNcrnaReference.MANIFEST_HEADER);
			origins.writer.println("input_ordinal0\tinput_id\tquery_line\tsource_line");
			counts.writer.println("input_ordinal0\tassembly_id\tinput_id\tfeatures\tquery_gff\tcanonical_new_gff");
			byte[] row=next(source);
			while(row!=null && !boundary(row)){
				require(row.length==0 || row[0]=='#', "Caller feature precedes its first input boundary at line "+sourceLine);
				rejectFasta(row);
				row=next(source);
			}
			for(final Binding binding:order){
				require(row!=null && boundary(row), "Missing input boundary for ordinal "+binding.ordinal);
				parser.set(row);
				require(parser.terms()==2 && parser.termEquals(NcrnaLocusComparison.INPUT_HEADER, 0)
					&& parser.parseInt(1)==binding.ordinal, "Input boundaries must be contiguous in declared order: expected="+binding.ordinal+" source_line="+sourceLine);
				final String query=new File(directory, "queries/"+binding.ordinal+".gff").getAbsolutePath();
				long features=0, queryLine=1;
				try(Output out=new Output(query)){
					out.writer.println("##gff-version 3");
					for(row=next(source); row!=null && !boundary(row); row=next(source)){
						rejectFasta(row);
						queryLine++;
						if(row.length>0 && row[0]!='#'){
							parser.set(row);
							// Format failure context only on failure, never once per saved feature.
							if(parser.terms()!=9 || parser.length(0)==0){
								throw new IllegalArgumentException("Caller feature requires nine fields and a contig: source_line="+sourceLine);
							}
							features++;
							origins.writer.print(buffer.clear().append(binding.ordinal).tab().append(binding.id).tab()
								.append(queryLine).tab().append(sourceLine).nl());
						}
						out.writer.println(row);
					}
				}
				queries.writer.print(buffer.clear().append(binding.id).tab().append(binding.fasta).tab().append(query).nl());
				references.writer.print(buffer.clear().append(binding.id).tab().append(binding.fasta).tab().append(binding.gff).nl());
				counts.writer.print(buffer.clear().append(binding.ordinal).tab().append(binding.assembly).tab().append(binding.id)
					.tab().append(features).tab().append(query).tab().append(binding.gff).nl());
				totalFeatures+=features;
			}
			require(row==null, "Saved GFF contains more input boundaries than its declared panel");
		}
		try(Output provenance=output("complete.tsv")){
			provenance.writer.println("#ORDINAL_NCRNA_CALL_BINDING_V2");
			provenance.writer.println("inputs\tfeatures\tsource_gff\tqualification");
			provenance.writer.print(buffer.clear().append(order.size()).tab().append(totalFeatures).tab().append(input)
				.append("\tDECLARED_CONFIG_AND_ORDINAL_BINDINGS_NOT_INDEPENDENT_PRODUCER_PROOF\n"));
		}
		outstream.println("Bound "+totalFeatures+" unchanged feature rows across "+order.size()+" inputs, including empty inputs.");
	}

	/** Validates the complete panel bijection before any call is attributed.
	 * Dense ordinals and unique global IDs make every declared panel ID reachable
	 * exactly once. Config order independently binds each ordinal to its FASTA path.
	 * This validates the supplied panel, not the biological provenance of its IDs.
	 */
	private void loadBindings(){
		final ObjectSet<String> ids=new ObjectSet<String>(String.class), fastas=new ObjectSet<String>(String.class), gffs=new ObjectSet<String>(String.class);
		try(Input in=new Input(bindings)){
			RrnaResourceIO.header(in.reader, BINDING_HEADER, bindings);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				parser.set(row);
				require(parser.terms()==5, "Panel source binding needs five fields at expected ordinal="+order.size());
				final int ordinal=parser.parseInt(0);
				final String id=parser.parseString(2);
				require(ordinal==order.size(), "Missing or repeated panel ordinal="+ordinal+" input_id="+id+" expected_ordinal="+order.size());
				final Binding binding=new Binding(order.size(), parser.parseString(1), parser.parseString(2), parser.parseString(3), parser.parseString(4));
				require(!id.isEmpty() && ids.add(id), "Empty or duplicate panel target: ordinal="+ordinal+" input_id="+id);
				require(!binding.assembly.isEmpty() && fastas.add(binding.fasta) && gffs.add(binding.gff)
					&& new File(binding.fasta).isAbsolute() && new File(binding.gff).isAbsolute(),
					"Repeated or invalid exact source paths: ordinal="+ordinal+" input_id="+id);
				order.add(binding);
			}
		}
		require(order.size()==expectedInputs, "Incomplete panel bijection: next_ordinal="+order.size()+" input_id=UNMAPPED expected_inputs="+expectedInputs);
	}

	/** Confirms the original closed config's complete input order and saved output path. */
	private void checkConfig(){
		int inputs=0, outputs=0;
		try(ConfigInput config=new ConfigInput(callerConfig)){
			while(config.next()){
				final String a=config.key.toLowerCase(Locale.ROOT), b=config.value;
				if(a.equals("in") || a.equals("in1")){
					require(inputs<order.size() && order.get(inputs).fasta.equals(b), "Original caller input order differs at ordinal "+inputs);
					inputs++;
				}else if(a.equals("out") || a.equals("out1")){
					require(++outputs==1 && input.equals(b), "Caller output must be this exact saved GFF");
				}else{
					require(!a.equals("config") && !a.equals("in2"), "Nested configs/paired inputs do not provide this flat source binding");
				}
			}
		}
		require(inputs==order.size() && outputs==1, "Caller config must declare every input and exactly one saved output");
	}

	private byte[] next(final Input input){
		final byte[] row=input.reader.nextLine();
		if(row!=null){sourceLine++;}
		return row;
	}
	private static boolean boundary(final byte[] row){
		if(row.length<BOUNDARY.length){return false;}
		for(int i=0; i<BOUNDARY.length; i++){if(row[i]!=BOUNDARY[i]){return false;}}
		return true;
	}
	private static void rejectFasta(final byte[] row){
		require(row.length<7 || row[0]!='#' || !shared.Tools.startsWith(row, "##FASTA"), "Saved caller GFF must not embed sequence payloads");
	}
	private Output output(final String name){return new Output(new File(directory, name).getPath());}
	private static void require(final boolean condition, final String why){if(!condition){throw new IllegalArgumentException(why);}}

	private static final class Binding{
		Binding(final int ordinal_, final String assembly_, final String id_, final String fasta_, final String gff_){
			ordinal=ordinal_;
			assembly=assembly_;
			id=id_;
			fasta=fasta_;
			gff=gff_;
		}
		final int ordinal;
		final String assembly, id, fasta, gff;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String input, bindings, callerConfig, directory;
	private int expectedInputs=-1;
	private long sourceLine, totalFeatures;
	private final PrintStream outstream;
	private final ArrayList<Binding> order=new ArrayList<Binding>();
	private final LineParser1 parser=new LineParser1('\t');
	private final ByteBuilder buffer=new ByteBuilder();
	private static final String BINDING_HEADER="input_ordinal0\tassembly_id\tinput_id\tfasta\tcanonical_new_gff";
	private static final byte[] BOUNDARY=NcrnaLocusComparison.INPUT_HEADER.getBytes(java.nio.charset.StandardCharsets.US_ASCII);
}
